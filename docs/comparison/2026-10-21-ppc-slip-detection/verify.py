"""Verify frozen code/Doppler slip runs; optionally recompute reference errors."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
import runpy


def sha(data):
    return hashlib.sha256(data).hexdigest()


def worst_positions(rows, references):
    result = []
    for r in rows:
        key = int(r["week"]), Decimal(r["tow"])
        if r["status"] != "FLOAT" or key not in references:
            continue
        reference = references[key]
        error = math.dist([float(r[k]) for k in ["x_m", "y_m", "z_m"]],
                          [float(reference[k]) for k in ["ECEF X (m)", "ECEF Y (m)", "ECEF Z (m)"]])
        result.append(dict(error_3d_m=error, output_row=r))
    return sorted(result, key=lambda r: r["error_3d_m"], reverse=True)[:10]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root", type=Path)
    parser.add_argument("--record-worst", action="store_true")
    args = parser.parse_args()
    archive = Path(__file__).resolve().parent
    repo = archive.parents[2]
    baseline = archive.parent / "2026-10-21-native-rtk-common-inputs"
    scoring = runpy.run_path(str(repo / "scripts/run_ppc_rtk_float.py"))
    slip = runpy.run_path(str(repo / "scripts/run_ppc_slip.py"))
    manifest = json.loads((archive / "predeclared-inputs.json").read_text())
    report = json.loads((archive / "report.json").read_text())
    assert not report["practical_PPC_parity_accepted"] and not manifest["other_thresholds_changed"]
    assert manifest["code_threshold_m"] == 5 and manifest["doppler_threshold_m"] == 0.2
    assert manifest["max_doppler_interval_s"] == 5
    assert set(report["cases"]) == set(manifest["cases"])
    for name, expected in manifest["source_sha256"].items():
        assert sha((repo / name).read_bytes()) == expected, name
    totals = {p: 0 for p in manifest["profiles"]}
    worst = dict(post_run_only=True, reference_hypothesis="antenna", cases={})
    for case in manifest["cases"]:
        city, run = case.split("-")
        references = None
        if args.dataset_root:
            root = args.dataset_root / city / run
            references = {}
            for name in ["rover.obs", "base.obs", "base.nav", "reference.csv"]:
                k = next(k for k in manifest["input_sha256"] if k.endswith(f"/{city}/{run}/{name}"))
                assert sha((root / name).read_bytes()) == manifest["input_sha256"][k]
            with (root / "reference.csv").open() as stream:
                for raw in csv.DictReader(stream):
                    r = {k.strip(): v.strip() for k, v in raw.items()}
                    references[(int(r["GPS Week"]), Decimal(r["GPS TOW (s)"]))] = r
        decoded = {}
        worst["cases"][case] = {}
        for profile in manifest["profiles"]:
            stem = case + "." + profile
            execution = json.loads((archive / (stem + ".execution.json")).read_text())
            data = gzip.decompress((archive / (stem + ".stdout.csv.gz")).read_bytes())
            stderr = gzip.decompress((archive / (stem + ".stderr.log.gz")).read_bytes())
            assert execution["returncode"] == 0
            assert sha(data) == execution["csv_sha256"] and sha(stderr) == execution["stderr_sha256"]
            rows, actual = scoring["inspect"](data)
            decoded[profile] = rows
            if profile == "none":
                with gzip.open(baseline / (case + ".rust.stdout.csv.gz"), "rt") as stream:
                    assert rows == list(csv.DictReader(stream))
                actual["all_columns_identical_to_2972bf5_baseline"] = True
            else:
                actual["slip_evidence"] = slip["slip_summary"](rows)
            expected = report["cases"][case][profile]
            assert actual == {k: v for k, v in expected.items() if k != "reference_scores"}
            if references is not None:
                assert scoring["score"](rows, references, city) == expected["reference_scores"]
                worst["cases"][case][profile] = worst_positions(rows, references)
            totals[profile] += actual["status_counts"].get("FLOAT", 0)
        keys = lambda rows: [(r["kind"], int(r["week"]), float(r["tow"]).hex(), r["side"]) for r in rows]
        assert all(keys(rows) == keys(decoded["none"]) for rows in decoded.values())
        # These are observed equalities, not assumptions in the solver.
        for a, b in [("none", "cmc"), ("doppler", "both")]:
            columns = decoded["none"][0].keys()
            assert all({k: r[k] for k in columns} == {k: s[k] for k in columns}
                       for r, s in zip(decoded[a], decoded[b]))
    native = json.loads((archive / "native-slip-generation.json").read_text())
    assert native["compile_returncode"] == native["generate_returncode"] == 0 and native["cases"] == 28
    assert sha((repo / "scripts/upstream_rtk_slip.cpp").read_bytes()) == native["source_sha256"]
    assert sha((repo / "tests/fixtures/upstream_rtk_slip.csv").read_bytes()) == native["fixture_sha256"]
    q = json.loads((archive / "qualification.json").read_text())
    assert sha((archive / "predeclared-inputs.json").read_bytes()) == q["predeclared_inputs_sha256"]
    assert {c["name"] for c in q["checks"]} == {"rust", "examples", "fmt", "clippy", "release"}
    for check in q["checks"]:
        data = gzip.decompress((archive / check["log"]).read_bytes())
        assert check["returncode"] == 0 and sha(data) == check["uncompressed_sha256"]
        if check["name"] == "rust":
            assert sum(map(int, re.findall(rb"test result: ok\. (\d+) passed;", data))) == check["passed_tests"] == 291
    if args.record_worst:
        assert args.dataset_root
        path = archive / "worst-positions.json"
        assert not path.exists()
        path.write_text(json.dumps(worst, indent=2) + "\n")
    else:
        if args.dataset_root:
            assert worst == json.loads((archive / "worst-positions.json").read_text())
        for name, expected in json.loads((archive / "evidence-sha256.json").read_text()).items():
            assert sha((archive / name).read_bytes()) == expected, name
    print("Verified 24 Rust runs, unchanged legacy columns, 28 native cases, 291 tests; FLOAT:", totals)


if __name__ == "__main__":
    main()
