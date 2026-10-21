"""Verify retained FLOAT evidence; optionally recompute scores from reference CSVs."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
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
    helpers = runpy.run_path(str(repo / "scripts/run_ppc_rtk_float.py"))
    manifest = json.loads((archive / "predeclared-inputs.json").read_text())
    report = json.loads((archive / "report.json").read_text())
    assert not manifest["performance_accepted"] and not report["performance_accepted"]
    assert not manifest["whole_native_comparison"] and not manifest["fixes_enabled"]
    assert set(report["cases"]) == set(manifest["cases"])
    for name, expected in manifest["source_sha256"].items():
        assert helpers["sha"](repo / name) == expected, name
    total = {}
    for case in manifest["cases"]:
        ledger = None
        references = None
        if args.dataset_root:
            city, run = case.split("-")
            directory = args.dataset_root / city / run
            for name in ["rover.obs", "base.obs", "base.nav", "reference.csv"]:
                original = next(k for k in manifest["input_sha256"] if k.endswith(f"/{city}/{run}/{name}"))
                assert helpers["sha"](directory / name) == manifest["input_sha256"][original]
            with (directory / "reference.csv").open() as stream:
                references = {}
                for raw in csv.DictReader(stream):
                    row = {k.strip(): v.strip() for k, v in raw.items()}
                    key = (int(row["GPS Week"]), Decimal(row["GPS TOW (s)"]))
                    assert key not in references
                    references[key] = row
        for mode in manifest["modes"]:
            prefix = archive / (case + "." + mode)
            execution = json.loads(Path(str(prefix) + ".execution.json").read_text())
            data = gzip.decompress(Path(str(prefix) + ".csv.gz").read_bytes())
            stderr = gzip.decompress(Path(str(prefix) + ".stderr.log.gz").read_bytes())
            assert execution["returncode"] == 0
            assert hashlib.sha256(data).hexdigest() == execution["stdout_sha256"]
            assert hashlib.sha256(stderr).hexdigest() == execution["stderr_sha256"]
            rows, summary = helpers["inspect"](data)
            expected = report["cases"][case][mode]
            assert summary == {k: v for k, v in expected.items() if k != "reference_scores"}
            keys = [(r["kind"], r["week"], r["tow"], r["side"], r["rover_flag"], r["base_flag"]) for r in rows]
            if ledger is None:
                ledger = keys
            else:
                assert keys == ledger
            if references is not None:
                assert helpers["score"](rows, references, case.split("-")[0]) == expected["reference_scores"]
            total[mode] = total.get(mode, 0) + summary["status_counts"].get("FLOAT", 0)
    qualification = json.loads((archive / "qualification.json").read_text())
    assert hashlib.sha256((archive / "predeclared-inputs.json").read_bytes()).hexdigest() == qualification["predeclared_inputs_sha256"]
    assert {r["name"] for r in qualification["checks"]} == {"rust", "examples", "fmt", "clippy", "release"}
    for check in qualification["checks"]:
        assert check["returncode"] == 0
        data = gzip.decompress((archive / check["log"]).read_bytes())
        assert hashlib.sha256(data).hexdigest() == check["uncompressed_sha256"]
        if check["name"] == "rust":
            assert sum(map(int, re.findall(rb"test result: ok\. (\d+) passed;", data))) == check["passed_tests"] == 272
    for name, expected in qualification["test_generator_fixture_sha256"].items():
        assert helpers["sha"](repo / name) == expected
    generation = json.loads((archive / "native-moving-generation.json").read_text())
    assert generation["returncode"] == 0 and not generation["native_whole_rtk"]
    for name, expected in generation["regenerated_sha256"].items():
        assert helpers["sha"](repo / "tests/fixtures" / name) == expected
    print("Verified 12 runs and 272 retained test passes; FLOAT outputs:", total)
    print("Reference scores recomputed." if references is not None else
          "Reference scores retained only; supply --dataset-root to recompute.")


if __name__ == "__main__":
    main()
