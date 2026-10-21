"""Read-only archive verification; optional reference accuracy recomputation."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import io
from itertools import zip_longest
import json
import math
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
    runner = runpy.run_path(str(repo / "scripts/run_native_ppc_comparison.py"))
    helpers = runpy.run_path(str(repo / "scripts/run_ppc_rtk_float.py"))
    manifest = json.loads((archive / "predeclared-inputs.json").read_text())
    report = json.loads((archive / "report.json").read_text())
    assert report["native_whole_RTKProcessor"] and not report["practical_PPC_parity_accepted"]
    assert set(report["cases"]) == set(manifest["cases"])
    for name, expected in manifest["source_sha256"].items():
        assert sha((repo / name).read_bytes()) == expected, name
    previous = repo / "docs/comparison/2026-10-21-ppc-geometry-gate"
    inputs = pairs = native_rtk = rust_float = 0
    for case in manifest["cases"]:
        city, run = case.split("-")
        references = None
        if args.dataset_root:
            root = args.dataset_root / city / run
            for name in ["rover.obs", "base.obs", "base.nav", "reference.csv"]:
                key = next(k for k in manifest["input_sha256"] if k.endswith(f"/{city}/{run}/{name}"))
                assert sha((root / name).read_bytes()) == manifest["input_sha256"][key]
            references = {}
            with (root / "reference.csv").open() as stream:
                for raw in csv.DictReader(stream):
                    r = {k.strip(): v.strip() for k, v in raw.items()}
                    references[(int(r["GPS Week"]), Decimal(r["GPS TOW (s)"]))] = r
        decoded = {}
        for solver in ["native", "rust"]:
            stem = case + "." + solver
            execution = json.loads((archive / (stem + ".execution.json")).read_text())
            data = gzip.decompress((archive / (stem + ".stdout.csv.gz")).read_bytes())
            stderr = gzip.decompress((archive / (stem + ".stderr.log.gz")).read_bytes())
            assert execution["returncode"] == 0
            assert sha(data) == execution["csv_sha256"] and sha(stderr) == execution["stderr_sha256"]
            if solver == "native":
                rows, scored, actual = runner["native_summary"](data)
                assert b"input_drained=1" in stderr
                assert actual["valid_outputs_with_mismatched_time"] == 0
            else:
                rows, actual = helpers["inspect"](data)
                scored = rows
                with gzip.open(previous / (case + ".kinematic.csv.gz"), "rt") as stream:
                    assert rows == list(csv.DictReader(stream))
            decoded[solver] = rows
            expected = report["cases"][case][solver]
            assert actual == {k: v for k, v in expected.items() if k != "reference_scores"}
            if references is not None:
                assert helpers["score"](scored, references, city) == expected["reference_scores"]
        keys = lambda rows: [(r["kind"], int(r["week"]), float(r["tow"]).hex(), r["side"]) for r in rows]
        assert keys(decoded["native"]) == keys(decoded["rust"])
        comparison = json.loads((archive / (case + ".input-comparison.json")).read_text())
        assert comparison == report["cases"][case]["input_comparison"]
        assert comparison["all_decoded_values_and_derived_resets_equal"]
        assert all(re.fullmatch("[0-9a-f]{64}", comparison[k]) for k in
                   ["normalized_sha256", "native_ledger_sha256", "rust_ledger_sha256"])
        # Recoveries must match both raw ledger hashes from the primary run.
        ledgers = []
        for solver in ["native", "rust"]:
            data = gzip.decompress((archive / (case + "." + solver + ".input-ledger.csv.gz")).read_bytes())
            assert sha(data) == comparison[solver + "_ledger_sha256"]
            ledgers.append(csv.DictReader(io.StringIO(data.decode())))
        normalized = hashlib.sha256()
        count = 0
        for count, (left, right) in enumerate(zip_longest(*ledgers), 1):
            assert left is not None and right is not None
            key = runner["ledger_key"](left)
            assert key == runner["ledger_key"](right)
            normalized.update((json.dumps(key, separators=(",", ":")) + "\n").encode())
        assert count == comparison["rows"]
        assert normalized.hexdigest() == comparison["normalized_sha256"]
        replay = json.loads((archive / (case + ".input-retention-replay.json")).read_text())
        assert replay["native_returncode"] == replay["rust_ledger_returncode"] == 0
        assert replay["native_output_identical_to_primary"] and replay["input_ledgers_identical_to_primary"]
        original = json.loads((archive / (case + ".native.execution.json")).read_text())
        assert replay["native_csv_sha256"] == original["csv_sha256"]
        replay_stderr = gzip.decompress((archive / (case + ".input-retention-replay.stderr.log.gz")).read_bytes())
        assert sha(replay_stderr) == replay["native_stderr_sha256"]
        inputs += comparison["rows"]
        pairs += expected["paired_epochs"]
        native_rtk += report["cases"][case]["native"]["valid_RTK_positions"]
        rust_float += expected["status_counts"].get("FLOAT", 0)
    q = json.loads((archive / "qualification.json").read_text())
    assert sha((archive / "predeclared-inputs.json").read_bytes()) == q["predeclared_inputs_sha256"]
    assert {c["name"] for c in q["checks"]} == {"rust", "examples", "fmt", "clippy", "release"}
    for c in q["checks"]:
        data = gzip.decompress((archive / c["log"]).read_bytes())
        assert c["returncode"] == 0 and sha(data) == c["uncompressed_sha256"]
        if c["name"] == "rust":
            assert sum(map(int, re.findall(rb"test result: ok\. (\d+) passed;", data))) == c["passed_tests"] == 276
    evidence = json.loads((archive / "evidence-sha256.json").read_text())
    for name, expected_hash in evidence.items():
        assert sha((archive / name).read_bytes()) == expected_hash, name
    smoke = json.loads((archive / "synthetic-smoke.json").read_text())
    data = gzip.decompress((archive / "synthetic-smoke.csv.gz").read_bytes())
    assert sha(data) == smoke["native_csv_sha256"]
    assert smoke["returncode"] == 0 and smoke["epochs"] == 64
    assert smoke["statuses"] == {"FIX": 45, "FLOAT": 19}
    for name, expected_hash in smoke["fixture_sha256"].items():
        assert sha((repo / name).read_bytes()) == expected_hash
    truth = {}
    with (repo / "tests/fixtures/synthetic_rtk_moving_truth.csv").open() as stream:
        for r in csv.DictReader(stream):
            truth[(int(r["week"]), float(r["tow"]).hex())] = r
    errors = []
    for r in csv.DictReader(io.StringIO(data.decode())):
        assert r["kind"] == "pair" and r["status"] in {"FLOAT", "FIX"}
        t = truth[(int(r["week"]), float(r["tow"]).hex())]
        errors.append(math.sqrt(sum((float(r[k]) - float(t[k]))**2 for k in ["x_m", "y_m", "z_m"])))
    assert len(errors) == smoke["epochs"] and max(errors) == smoke["maximum_3d_error_m"]
    failed = archive.parent / "2026-10-21-native-rtk-common-inputs-format-attempt"
    first_manifest = json.loads((failed / "predeclared-inputs.json").read_text())
    assert sha((failed / "runner-before-format-fix.py").read_bytes()) == first_manifest["source_sha256"]["scripts/run_native_ppc_comparison.py"]
    failure = json.loads((failed / "format-guard-failure.json").read_text())
    assert not failure["complete"] and not failure["solver_profiles_changed"]
    assert failure["first_mismatch"]["native_float_hex"] == failure["first_mismatch"]["rust_float_hex"]
    if args.dataset_root:
        analysis = runpy.run_path(str(archive / "analyze.py"))
        assert analysis["analyze"](archive, args.dataset_root) == json.loads((archive / "diagnostics.json").read_text())
    print(f"Verified 12 solver executions, {pairs} pairs, {inputs} retained input-comparison rows; native RTK={native_rtk}, Rust FLOAT={rust_float}; 276 tests.")
    print("Both retained raw input ledgers independently compared, including all values and presence flags.")
    print("All reference hypotheses recomputed." if args.dataset_root else "Supply --dataset-root to recompute reference statistics.")


if __name__ == "__main__":
    main()
