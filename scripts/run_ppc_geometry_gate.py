"""Compare unmodified admission with an optional GDOP 50 gate on PPC FLOAT."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import runpy
import subprocess


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def geometry_summary(rows):
    valid = [r for r in rows if r["status"] == "FLOAT"]
    dop = [float(r["geometry_gdop"]) for r in valid if r["geometry_gdop"]]
    assert all(math.isfinite(v) and v > 0 for v in dop)
    assert all(bool(r["geometry_gdop"]) != bool(r["geometry_error"]) for r in valid)
    for row in valid:
        assert len(set(row["used_satellites"].split(";"))) == int(row["satellites"])
        trace = float(row["position_covariance_trace_m2"])
        diagonal = [float(row[k]) for k in ["position_qxx_m2", "position_qyy_m2", "position_qzz_m2"]]
        assert all(math.isfinite(v) and v > 0 for v in diagonal)
        assert abs(sum(diagonal) - trace) <= 1e-12 * max(1.0, trace)
    return dict(
        unavailable_geometry=len(valid)-len(dop),
        positions_above_gdop_50=sum(v > 50 for v in dop),
        maximum_gdop=max(dop, default=None),
        maximum_position_covariance_trace_m2=max(
            (float(r["position_covariance_trace_m2"]) for r in valid), default=None
        ),
    )


def worst_positions(rows, references):
    values = []
    for row in rows:
        if row["status"] != "FLOAT":
            continue
        key = (int(row["week"]), Decimal(row["tow"]))
        reference = references.get(key)
        if reference is None:
            continue
        error = math.dist(
            [float(row[k]) for k in ["x_m", "y_m", "z_m"]],
            [float(reference[k]) for k in ["ECEF X (m)", "ECEF Y (m)", "ECEF Z (m)"]],
        )
        values.append(dict(
            week=key[0], tow=str(key[1]), three_dimensional_error_m=error,
            satellites=int(row["satellites"]), used_satellites=row["used_satellites"],
            gdop=float(row["geometry_gdop"]) if row["geometry_gdop"] else None,
            posterior_position_std_norm_m=math.sqrt(float(row["position_covariance_trace_m2"])),
            code_rms_m=float(row["code_rms_m"]), phase_rms_m=float(row["phase_rms_m"]),
            nis=float(row["nis"]), ambiguity_resets=int(row["ambiguity_resets"]),
        ))
    return sorted(values, key=lambda r: r["three_dimensional_error_m"], reverse=True)[:10]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    dataset, binary, output = args.dataset_root.resolve(), args.binary.resolve(), args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    previous = repo/"docs/comparison/2026-10-21-ppc-kinematic-float"
    prior = json.loads((previous/"predeclared-inputs.json").read_text())
    helpers = runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))
    files = list((repo/"src").rglob("*.rs")) + [
        repo/"scripts/run_ppc_geometry_gate.py", repo/"scripts/run_ppc_rtk_float.py",
        repo/"examples/ppc_rtk_float.rs", repo/"Cargo.toml", repo/"Cargo.lock", repo/"rust-toolchain.toml",
    ]
    inputs = [dataset/case.split("-")[0]/case.split("-")[1]/name
              for case in prior["cases"] for name in ["rover.obs", "base.obs", "base.nav", "reference.csv"]]
    previous_outputs = {case: sha(previous/(case+".kinematic.csv.gz")) for case in prior["cases"]}
    manifest = dict(
        cases=prior["cases"], profiles=["kinematic", "kinematic-gdop50"],
        source_sha256={str(p.relative_to(repo)): sha(p) for p in files},
        input_sha256={str(p): sha(p) for p in inputs}, binary_sha256=sha(binary),
        prior_baseline_gzip_sha256=previous_outputs,
        base_ecef_m=prior["base_ecef_m"], seed_source=prior["seed_source"],
        geometry_gate_default=None, experiment_gdop_limit=50.0,
        gate_threshold_source="pinned SPPProcessor::SPPConfig::max_gdop default 50; not an upstream RTK gate",
        native_geometry_cases=16, native_whole_comparison=False, performance_accepted=False,
        reference_used_for_gate_or_solver=False, all_reference_hypotheses_retained=True,
        worst_ten_positions="post-execution diagnostics only; never a selection or gate criterion",
    )
    save(output/"predeclared-inputs.json", manifest)
    report = dict(native_whole_comparison=False, performance_accepted=False, cases={})
    worst = {}
    for case in manifest["cases"]:
        city, run = case.split("-")
        directory = dataset/city/run
        references = {}
        with (directory/"reference.csv").open() as stream:
            for raw in csv.DictReader(stream):
                row = {k.strip(): v.strip() for k, v in raw.items()}
                key = (int(row["GPS Week"]), Decimal(row["GPS TOW (s)"]))
                assert key not in references
                references[key] = row
        report["cases"][case] = {}
        worst[case] = {}
        for profile in manifest["profiles"]:
            assert all(sha(repo/name) == h for name, h in manifest["source_sha256"].items())
            assert sha(binary) == manifest["binary_sha256"]
            command = [str(binary), str(directory), profile, ",".join(map(str, manifest["base_ecef_m"][city]))]
            result = subprocess.run(command, capture_output=True, cwd=repo)
            prefix = str(output/(case+"."+profile))
            Path(prefix+".csv.gz").write_bytes(gzip.compress(result.stdout, mtime=0))
            Path(prefix+".stderr.log.gz").write_bytes(gzip.compress(result.stderr, mtime=0))
            save(Path(prefix+".execution.json"), dict(
                command=command, returncode=result.returncode,
                stdout_sha256=hashlib.sha256(result.stdout).hexdigest(),
                stderr_sha256=hashlib.sha256(result.stderr).hexdigest(),
            ))
            assert result.returncode == 0, result.stderr.decode()
            rows, summary = helpers["inspect"](result.stdout)
            if profile == "kinematic":
                assert sha(previous/(case+".kinematic.csv.gz")) == previous_outputs[case]
                old = list(csv.DictReader(io.StringIO(gzip.decompress(
                    (previous/(case+".kinematic.csv.gz")).read_bytes()
                ).decode())))
                assert len(rows) == len(old)
                assert all(all(row[k] == original[k] for k in original) for row, original in zip(rows, old))
                summary["all_legacy_columns_identical_to_e6502f4"] = True
            else:
                assert all(r["geometry_gdop"] and float(r["geometry_gdop"]) <= 50.0
                           for r in rows if r["status"] == "FLOAT")
            summary["geometry"] = geometry_summary(rows)
            summary["reference_scores"] = helpers["score"](rows, references, city)
            report["cases"][case][profile] = summary
            worst[case][profile] = worst_positions(rows, references)
            save(output/"report.json", report)
            save(output/"worst-positions.json", worst)
            print(case, profile, json.dumps({
                k: v for k, v in summary.items() if k != "reference_scores"
            }), flush=True)
    assert all(sha(Path(name)) == h for name, h in manifest["input_sha256"].items())


if __name__ == "__main__":
    main()
