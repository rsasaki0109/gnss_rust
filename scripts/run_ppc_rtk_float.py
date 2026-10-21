"""Frozen GPS L1 static/kinematic FLOAT comparison on all common PPC epochs."""
import argparse
import collections
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import subprocess


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def save(path, data):
    path.write_text(json.dumps(data, indent=2) + "\n")


def ecef(lat, lon, height):
    lat, lon = math.radians(lat), math.radians(lon)
    e2 = (1 / 298.257223563) * (2 - 1 / 298.257223563)
    radius = 6378137 / math.sqrt(1 - e2 * math.sin(lat) ** 2)
    return [(radius + height) * math.cos(lat) * math.cos(lon),
            (radius + height) * math.cos(lat) * math.sin(lon),
            (radius * (1 - e2) + height) * math.sin(lat)]


def reference_targets(row, city):
    xyz = [float(row[k]) for k in ["ECEF X (m)", "ECEF Y (m)", "ECEF Z (m)"]]
    lat, lon = map(math.radians, [float(row["Latitude (deg)"]), float(row["Longitude (deg)"])])
    targets = {"reference_is_antenna": xyz}
    arms = {"reference_is_IMU_nagoya": [0.593, -0.670, -1.216]} if city == "nagoya" else {
        "reference_is_IMU_tokyo_top_README": [0.31, 0.0, -0.55],
        "reference_is_IMU_tokyo_city_README": [0.593, -0.670, -1.216],
    }
    r, p, h = [math.radians(float(row[k])) for k in ["Roll (deg)", "Pitch (deg)", "Heading (deg)"]]
    for name, (x, y, z) in arms.items():
        y, z = math.cos(r)*y-math.sin(r)*z, math.sin(r)*y+math.cos(r)*z
        x, z = math.cos(p)*x+math.sin(p)*z, -math.sin(p)*x+math.cos(p)*z
        north, east = math.cos(h)*x-math.sin(h)*y, math.sin(h)*x+math.cos(h)*y
        up = -z
        delta = [-math.sin(lon)*east-math.sin(lat)*math.cos(lon)*north+math.cos(lat)*math.cos(lon)*up,
                 math.cos(lon)*east-math.sin(lat)*math.sin(lon)*north+math.cos(lat)*math.sin(lon)*up,
                 math.cos(lat)*north+math.sin(lat)*up]
        targets[name] = [a+b for a, b in zip(xyz, delta)]
    return targets, lat, lon


def score(rows, references, city):
    errors = collections.defaultdict(list)
    unavailable = 0
    for row in rows:
        if row["status"] != "FLOAT":
            continue
        key = (int(row["week"]), Decimal(row["tow"]))
        if key not in references:
            unavailable += 1
            continue
        targets, lat, lon = reference_targets(references[key], city)
        xyz = [float(row[k]) for k in ["x_m", "y_m", "z_m"]]
        for name, target in targets.items():
            dx, dy, dz = [a-b for a, b in zip(xyz, target)]
            east = -math.sin(lon)*dx+math.cos(lon)*dy
            north = -math.sin(lat)*math.cos(lon)*dx-math.sin(lat)*math.sin(lon)*dy+math.cos(lat)*dz
            up = math.cos(lat)*math.cos(lon)*dx+math.cos(lat)*math.sin(lon)*dy+math.sin(lat)*dz
            errors[name].append((math.hypot(east, north), up, math.sqrt(dx*dx+dy*dy+dz*dz)))
    result = {}
    for name, values in errors.items():
        def quantile(data, fraction):
            data = sorted(data)
            position = (len(data)-1)*fraction
            low = int(position)
            return data[low]+(data[min(low+1, len(data)-1)]-data[low])*(position-low)
        result[name] = dict(
            epochs=len(values),
            horizontal_rms_m=math.sqrt(sum(v[0]**2 for v in values)/len(values)),
            vertical_rms_m=math.sqrt(sum(v[1]**2 for v in values)/len(values)),
            three_dimensional_rms_m=math.sqrt(sum(v[2]**2 for v in values)/len(values)),
            horizontal_p50_m=quantile([v[0] for v in values], 0.5),
            horizontal_p95_m=quantile([v[0] for v in values], 0.95),
            maximum_3d_m=max(v[2] for v in values),
        )
    return dict(valid_positions_without_reference=unavailable, hypotheses=result)


def inspect(data):
    rows = list(csv.DictReader(io.StringIO(data.decode())))
    assert rows
    status = collections.Counter(r["status"] for r in rows)
    assert set(status) <= {"FLOAT", "INVALID", "EVENT", "UNMATCHED"}
    pairs = [r for r in rows if r["kind"] == "pair"]
    assert all(r["side"] == "both" and abs(float(r["base_dt_s"])) <= 1e-7 for r in pairs)
    assert all(r["status"] != "UNMATCHED" for r in pairs)
    for row in rows:
        if row["status"] == "FLOAT":
            assert all(math.isfinite(float(row[k])) for k in [
                "x_m", "y_m", "z_m", "nis", "code_rms_m", "phase_rms_m"
            ])
            assert int(row["satellites"]) >= 4 and float(row["nis"]) >= 0
        else:
            assert all(not row[k] for k in ["x_m", "y_m", "z_m"])
    summary = dict(
        output_rows=len(rows), paired_epochs=len(pairs), status_counts=dict(status),
        rover_only=sum(r["side"] == "rover" for r in rows),
        base_only=sum(r["side"] == "base" for r in rows),
        invalid_reasons=dict(collections.Counter(r["reason"] for r in rows if r["status"] == "INVALID")),
        unmatched_continuity_notices=sum(r["reason"] == "continuity_reset_pending" for r in rows),
    )
    summary["rover_epochs"] = len(pairs)+summary["rover_only"]
    summary["base_epochs"] = len(pairs)+summary["base_only"]
    summary["float_coverage_among_pairs"] = status["FLOAT"]/len(pairs)
    summary["float_coverage_of_original_rover_epochs"] = status["FLOAT"]/summary["rover_epochs"]
    return rows, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    dataset, binary, output = args.dataset_root.resolve(), args.binary.resolve(), args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    cases = [c+"-"+r for c in ["nagoya", "tokyo"] for r in ["run1", "run2", "run3"]]
    files = [dataset/c/r/name for c in ["nagoya", "tokyo"] for r in ["run1", "run2", "run3"]
             for name in ["rover.obs", "base.obs", "base.nav", "reference.csv"]]
    sources = list((repo/"src").rglob("*.rs")) + [
        repo/"examples/ppc_rtk_float.rs", repo/"scripts/run_ppc_rtk_float.py",
        repo/"Cargo.toml", repo/"Cargo.lock", repo/"rust-toolchain.toml",
    ]
    base_geodetic = {"tokyo": [35.66633426, 139.79220181, 59.82],
                     "nagoya": [35.13470947, 136.97757427, 104.718]}
    manifest = dict(
        cases=cases, modes=["static", "kinematic"],
        source_sha256={str(p.relative_to(repo)): sha(p) for p in sources},
        input_sha256={str(p): sha(p) for p in files}, binary_sha256=sha(binary),
        base_phase_centre_geodetic_degrees_m=base_geodetic,
        base_ecef_m={city: ecef(*point) for city, point in base_geodetic.items()},
        base_position_source="PPC top-level README surveyed phase centre; no antenna correction",
        seed_source="rover RINEX approximate position, not reference.csv",
        matching="common epochs within 1e-7 s only; all unmatched epochs reported",
        unmatched_LLI_power_slip_evidence="conservative all-arc reset on next accepted pair",
        fixes_enabled=False, satellite_signals="GPS 1C only", thresholds_changed=False,
        whole_native_comparison=False, performance_accepted=False,
        reference_used_to_choose_config=False,
        reference_hypotheses="antenna; IMU with both conflicting Tokyo README lever arms",
    )
    save(output/"predeclared-inputs.json", manifest)
    report = dict(whole_native_comparison=False, performance_accepted=False, cases={})
    for case in cases:
        city, run = case.split("-")
        data = dataset/city/run
        references = {}
        with (data/"reference.csv").open() as stream:
            for raw in csv.DictReader(stream):
                row = {k.strip(): v.strip() for k, v in raw.items()}
                key = (int(row["GPS Week"]), Decimal(row["GPS TOW (s)"]))
                assert key not in references
                references[key] = row
        report["cases"][case] = {}
        ledger = None
        for mode in manifest["modes"]:
            assert all(sha(repo/name) == h for name, h in manifest["source_sha256"].items())
            assert sha(binary) == manifest["binary_sha256"]
            command = [str(binary), str(data), mode, ",".join(map(str, manifest["base_ecef_m"][city]))]
            result = subprocess.run(command, cwd=repo, capture_output=True)
            prefix = output/(case+"."+mode)
            Path(str(prefix)+".csv.gz").write_bytes(gzip.compress(result.stdout, mtime=0))
            Path(str(prefix)+".stderr.log.gz").write_bytes(gzip.compress(result.stderr, mtime=0))
            save(Path(str(prefix)+".execution.json"), dict(
                command=command, returncode=result.returncode,
                stdout_sha256=hashlib.sha256(result.stdout).hexdigest(),
                stderr_sha256=hashlib.sha256(result.stderr).hexdigest(),
            ))
            assert result.returncode == 0, result.stderr.decode()
            rows, summary = inspect(result.stdout)
            keys = [(r["kind"], r["week"], r["tow"], r["side"], r["rover_flag"], r["base_flag"]) for r in rows]
            if ledger is None:
                ledger = keys
            else:
                assert ledger == keys
            summary["reference_scores"] = score(rows, references, city)
            report["cases"][case][mode] = summary
            save(output/"report.json", report)
            print(case, mode, json.dumps({k:v for k,v in summary.items() if k!="reference_scores"}), flush=True)
    assert all(sha(Path(name)) == h for name,h in manifest["input_sha256"].items())
    assert all(sha(repo/name) == h for name,h in manifest["source_sha256"].items())


if __name__ == "__main__":
    main()
