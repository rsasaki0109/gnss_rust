"""Verify retained geometry experiments; optionally recompute reference scores."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import runpy


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root", type=Path)
    args = parser.parse_args()
    archive = Path(__file__).resolve().parent
    repo = archive.parents[2]
    helpers = runpy.run_path(str(repo/"scripts/run_ppc_geometry_gate.py"))
    scoring = runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))
    manifest = json.loads((archive/"predeclared-inputs.json").read_text())
    report = json.loads((archive/"report.json").read_text())
    worst = json.loads((archive/"worst-positions.json").read_text())
    assert not report["performance_accepted"] and not report["native_whole_comparison"]
    for name, expected in manifest["source_sha256"].items():
        assert helpers["sha"](repo/name) == expected
    totals = {}
    previous = repo/"docs/comparison/2026-10-21-ppc-kinematic-float"
    for case in manifest["cases"]:
        references = None
        if args.dataset_root:
            city, run = case.split("-")
            directory = args.dataset_root/city/run
            for name in ["rover.obs", "base.obs", "base.nav", "reference.csv"]:
                key = next(k for k in manifest["input_sha256"] if k.endswith(f"/{city}/{run}/{name}"))
                assert helpers["sha"](directory/name) == manifest["input_sha256"][key]
            references = {}
            with (directory/"reference.csv").open() as stream:
                for raw in csv.DictReader(stream):
                    row = {k.strip(): v.strip() for k, v in raw.items()}
                    references[(int(row["GPS Week"]), Decimal(row["GPS TOW (s)"]))] = row
        for profile in manifest["profiles"]:
            prefix = str(archive/(case+"."+profile))
            execution = json.loads(Path(prefix+".execution.json").read_text())
            data = gzip.decompress(Path(prefix+".csv.gz").read_bytes())
            stderr = gzip.decompress(Path(prefix+".stderr.log.gz").read_bytes())
            assert execution["returncode"] == 0
            assert hashlib.sha256(data).hexdigest() == execution["stdout_sha256"]
            assert hashlib.sha256(stderr).hexdigest() == execution["stderr_sha256"]
            rows, actual = scoring["inspect"](data)
            if profile == "kinematic":
                old = previous/(case+".kinematic.csv.gz")
                assert helpers["sha"](old) == manifest["prior_baseline_gzip_sha256"][case]
                original = list(csv.DictReader(io.StringIO(gzip.decompress(old.read_bytes()).decode())))
                assert len(rows) == len(original)
                assert all(all(r[k] == q[k] for k in q) for r, q in zip(rows, original))
                actual["all_legacy_columns_identical_to_e6502f4"] = True
            else:
                assert all(float(r["geometry_gdop"]) <= 50 for r in rows if r["status"] == "FLOAT")
            actual["geometry"] = helpers["geometry_summary"](rows)
            expected = report["cases"][case][profile]
            assert actual == {k: v for k, v in expected.items() if k != "reference_scores"}
            if references is not None:
                assert scoring["score"](rows, references, case.split("-")[0]) == expected["reference_scores"]
                assert helpers["worst_positions"](rows, references) == worst[case][profile]
            totals[profile] = totals.get(profile, 0) + actual["status_counts"].get("FLOAT", 0)
    qualification = json.loads((archive/"qualification.json").read_text())
    assert hashlib.sha256((archive/"predeclared-inputs.json").read_bytes()).hexdigest() == qualification["predeclared_inputs_sha256"]
    assert {c["name"] for c in qualification["checks"]} == {"rust", "examples", "fmt", "clippy", "release"}
    for check in qualification["checks"]:
        assert check["returncode"] == 0
        data = gzip.decompress((archive/check["log"]).read_bytes())
        assert hashlib.sha256(data).hexdigest() == check["uncompressed_sha256"]
        if check["name"] == "rust":
            assert sum(map(int, re.findall(rb"test result: ok\. (\d+) passed;", data))) == check["passed_tests"] == 276
    for name, expected in qualification["additional_source_fixture_sha256"].items():
        assert helpers["sha"](repo/name) == expected
    native = json.loads((archive/"native-geometry-generation.json").read_text())
    assert native["returncode"] == 0 and native["cases"] == 16
    assert helpers["sha"](repo/"tests/fixtures/upstream_geometry.csv") == native["fixture_sha256"]
    probes = json.loads((archive/"raw-phase-probes.json").read_text())
    assert helpers["sha"](repo/"examples/ppc_phase_diagnostics.rs") == probes["source_sha256"]
    for execution in probes["executions"]:
        assert execution["returncode"] == 0
        data = gzip.decompress((archive/(execution["case"]+".raw-phase-probe.csv.gz")).read_bytes())
        assert hashlib.sha256(data).hexdigest() == execution["stdout_sha256"]
    print("Verified 12 runs, unchanged legacy baseline columns, 276 retained test passes; FLOAT:", totals)
    print("Reference scores/worst positions recomputed." if references is not None else
          "Reference scores retained; supply --dataset-root to recompute.")


if __name__ == "__main__":
    main()
