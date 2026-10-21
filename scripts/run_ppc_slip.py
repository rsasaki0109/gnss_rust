"""Freeze PPC GPS 1C runs with native-default code/Doppler slip thresholds."""
import argparse
import collections
import csv
from decimal import Decimal
import gzip
import json
from pathlib import Path
import runpy
import subprocess


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def slip_summary(rows):
    totals = collections.Counter()
    for r in rows:
        if r["status"] != "FLOAT":
            continue
        for k in ["code_slips", "doppler_slips", "code_tests", "doppler_tests"]:
            value = int(r[k])
            assert 0 <= value <= int(r["satellites"])
            totals[k] += value
        assert int(r["code_slips"]) <= int(r["code_tests"])
        assert int(r["doppler_slips"]) <= int(r["doppler_tests"])
        union = set(r["slip_satellites"].split(";")) if r["slip_satellites"] else set()
        assert max(int(r["code_slips"]), int(r["doppler_slips"])) <= len(union) <= int(r["code_slips"]) + int(r["doppler_slips"])
        assert len(union) <= int(r["ambiguity_resets"])
        totals["epochs_with_detected_slips"] += bool(union)
        totals["satellite_arc_resets_from_detection"] += len(union)
        totals["epochs_with_invalid_doppler"] += bool(r["doppler_errors"])
    return dict(totals)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    baseline = repo / "docs/comparison/2026-10-21-native-rtk-common-inputs"
    reference = repo / "docs/comparison/2026-10-21-native-default-seed"
    control = json.loads((baseline / "predeclared-inputs.json").read_text())
    scoring = runpy.run_path(str(repo / "scripts/run_ppc_rtk_float.py"))
    capture = runpy.run_path(str(repo / "scripts/run_native_ppc_comparison.py"))["retain"]
    sha = scoring["sha"]
    output, dataset, binary = args.output_dir.resolve(), args.dataset_root.resolve(), args.binary.resolve()
    output.mkdir(parents=True, exist_ok=False)
    sources = list((repo / "src").rglob("*.rs")) + [repo / n for n in [
        "scripts/run_ppc_slip.py", "scripts/run_ppc_rtk_float.py", "scripts/run_native_ppc_comparison.py",
        "examples/ppc_rtk_float.rs", "tests/rtk_float.rs", "tests/rtk_slip.rs",
        "scripts/upstream_rtk_slip.cpp", "tests/fixtures/upstream_rtk_slip.csv",
        "Cargo.toml", "Cargo.lock", "rust-toolchain.toml",
    ]]
    inputs = {}
    for case in control["cases"]:
        city, run = case.split("-")
        for name in ["rover.obs", "base.obs", "base.nav", "reference.csv"]:
            p = dataset / city / run / name
            k = next(k for k in control["input_sha256"] if k.endswith(f"/{city}/{run}/{name}"))
            assert sha(p) == control["input_sha256"][k]
            inputs[str(p)] = sha(p)
    manifest = dict(
        source_sha256={str(p.relative_to(repo)): sha(p) for p in sources}, input_sha256=inputs,
        binary_sha256=sha(binary), cases=control["cases"], base_ecef_m=control["base_ecef_m"],
        profiles={"none": "kinematic", "cmc": "kinematic-cmc", "doppler": "kinematic-doppler", "both": "kinematic-slips"},
        code_threshold_m=5.0, doppler_threshold_m=0.2, max_doppler_interval_s=5.0,
        threshold_source="pinned native RTKConfig defaults, no dataset truth tuning",
        slip_history="only accepted epoch satellites; no comparisons across known breaks, failure, loss or long phase gap",
        SPP_seeding=False, kinematic_FIX=False, other_thresholds_changed=False,
        native_reference_report_sha256=sha(reference / "report.json"),
        native_whole_solver_rerun=False, native_history_policy_equivalence_claimed=False,
        reference_used_to_choose_configuration=False, timing_memory_acceptance=False,
        practical_PPC_parity_accepted=False,
    )
    save(output / "predeclared-inputs.json", manifest)
    report = dict(practical_PPC_parity_accepted=False, cases={})
    for case in control["cases"]:
        city, run = case.split("-")
        root = dataset / city / run
        references = {}
        with (root / "reference.csv").open() as stream:
            for raw in csv.DictReader(stream):
                r = {k.strip(): v.strip() for k, v in raw.items()}
                references[(int(r["GPS Week"]), Decimal(r["GPS TOW (s)"]))] = r
        report["cases"][case] = {}
        ledger = None
        for profile, mode in manifest["profiles"].items():
            assert all(sha(repo / p) == h for p, h in manifest["source_sha256"].items())
            assert sha(binary) == manifest["binary_sha256"]
            command = [str(binary), str(root), mode, ",".join(map(str, control["base_ecef_m"][city]))]
            result = subprocess.run(command, capture_output=True, cwd=repo)
            capture(output, case + "." + profile, command, result.stdout, result.stderr, result.returncode)
            assert result.returncode == 0, result.stderr.decode()
            rows, metrics = scoring["inspect"](result.stdout)
            keys = [(r["kind"], int(r["week"]), float(r["tow"]).hex(), r["side"]) for r in rows]
            assert ledger is None or ledger == keys
            ledger = keys
            if profile == "none":
                with gzip.open(baseline / (case + ".rust.stdout.csv.gz"), "rt") as stream:
                    assert rows == list(csv.DictReader(stream))
                metrics["all_columns_identical_to_2972bf5_baseline"] = True
            else:
                metrics["slip_evidence"] = slip_summary(rows)
            metrics["reference_scores"] = scoring["score"](rows, references, city)
            report["cases"][case][profile] = metrics
            save(output / "report.json", report)
            print(case, profile, metrics["status_counts"].get("FLOAT", 0), metrics.get("slip_evidence", {}), flush=True)
    assert all(sha(Path(p)) == h for p, h in inputs.items())
    assert all(sha(repo / p) == h for p, h in manifest["source_sha256"].items())


if __name__ == "__main__":
    main()
