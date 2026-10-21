"""Verify Rust SPP seed experiments; optionally recompute every reference score."""
import argparse
import collections
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import runpy


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root", type=Path)
    args = parser.parse_args()
    archive = Path(__file__).resolve().parent
    repo = archive.parents[2]
    prior_archive = archive.parent / "2026-10-21-native-rtk-common-inputs"
    native_archive = archive.parent / "2026-10-21-native-default-seed"
    scoring = runpy.run_path(str(repo / "scripts/run_ppc_rtk_float.py"))
    manifest = json.loads((archive / "predeclared-inputs.json").read_text())
    report = json.loads((archive / "report.json").read_text())
    native = json.loads((native_archive / "report.json").read_text())
    assert not report["practical_PPC_parity_accepted"] and not manifest["changed_thresholds"]
    assert set(report["cases"]) == set(manifest["cases"])
    assert sha((native_archive / "predeclared-inputs.json").read_bytes()) == manifest["native_manifest_sha256"]
    assert sha((native_archive / "report.json").read_bytes()) == manifest["native_report_sha256"]
    for name, expected in manifest["source_sha256"].items():
        assert sha((repo / name).read_bytes()) == expected, name
    totals = collections.Counter()
    for case in manifest["cases"]:
        city, run = case.split("-")
        references = None
        if args.dataset_root:
            directory = args.dataset_root / city / run
            for name in ["rover.obs", "base.obs", "base.nav", "reference.csv"]:
                key = next(k for k in manifest["input_sha256"] if k.endswith(f"/{city}/{run}/{name}"))
                assert sha((directory / name).read_bytes()) == manifest["input_sha256"][key]
            references = {}
            with (directory / "reference.csv").open() as stream:
                for raw in csv.DictReader(stream):
                    r = {k.strip(): v.strip() for k, v in raw.items()}
                    references[(int(r["GPS Week"]), Decimal(r["GPS TOW (s)"]))] = r
        ledger = None
        for profile in manifest["profiles"]:
            stem = case + "." + profile
            execution = json.loads((archive / (stem + ".execution.json")).read_text())
            data = gzip.decompress((archive / (stem + ".stdout.csv.gz")).read_bytes())
            stderr = gzip.decompress((archive / (stem + ".stderr.log.gz")).read_bytes())
            assert execution["returncode"] == 0
            assert sha(data) == execution["csv_sha256"] and sha(stderr) == execution["stderr_sha256"]
            rows, actual = scoring["inspect"](data)
            keys = [(r["kind"], int(r["week"]), float(r["tow"]).hex(), r["side"]) for r in rows]
            assert ledger is None or keys == ledger
            ledger = keys
            if profile == "previous":
                with gzip.open(prior_archive / (case + ".rust.stdout.csv.gz"), "rt") as stream:
                    assert rows == list(csv.DictReader(stream))
                actual["all_columns_identical_to_2972bf5_baseline"] = True
            else:
                accepted = [r for r in rows if r["status"] == "FLOAT"]
                assert all(r["seed_source"] in {"InitialPosition", "PreviousRtk", "SinglePoint"} for r in accepted)
                assert all((r["seed_source"] == "SinglePoint") == (not r["spp_failure"]) for r in accepted)
                actual["accepted_seed_sources"] = dict(collections.Counter(r["seed_source"] for r in accepted))
                actual["accepted_spp_fallback_errors"] = dict(collections.Counter(r["spp_failure"] for r in accepted if r["spp_failure"]))
            expected = report["cases"][case][profile]
            assert actual == {k: v for k, v in expected.items() if k != "reference_scores"}
            if references is not None:
                assert scoring["score"](rows, references, city) == expected["reference_scores"]
            totals[profile] += actual["status_counts"].get("FLOAT", 0)
        assert report["cases"][case]["native_default_SPP_seed"] == native["cases"][case]["native_default_SPP_seed"]
    q = json.loads((archive / "qualification.json").read_text())
    assert sha((archive / "predeclared-inputs.json").read_bytes()) == q["predeclared_inputs_sha256"]
    assert {c["name"] for c in q["checks"]} == {"rust", "examples", "fmt", "clippy", "release"}
    for check in q["checks"]:
        data = gzip.decompress((archive / check["log"]).read_bytes())
        assert check["returncode"] == 0 and sha(data) == check["uncompressed_sha256"]
        if check["name"] == "rust":
            assert sum(map(int, re.findall(rb"test result: ok\. (\d+) passed;", data))) == check["passed_tests"] == 280
    for name, expected in json.loads((archive / "evidence-sha256.json").read_text()).items():
        assert sha((archive / name).read_bytes()) == expected, name
    print("Verified 12 Rust executions, unchanged baseline CSVs, 280 tests; FLOAT:", dict(totals))
    print("All reference hypotheses recomputed." if args.dataset_root else "Supply --dataset-root to recompute reference scores.")


if __name__ == "__main__":
    main()
