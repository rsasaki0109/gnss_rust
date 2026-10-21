"""Freeze and compare Rust kinematic RTK with previous-RTK and SPP seeds."""
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    previous = repo / "docs/comparison/2026-10-21-native-rtk-common-inputs"
    native_archive = repo / "docs/comparison/2026-10-21-native-default-seed"
    control = json.loads((previous / "predeclared-inputs.json").read_text())
    native_report = json.loads((native_archive / "report.json").read_text())
    scoring = runpy.run_path(str(repo / "scripts/run_ppc_rtk_float.py"))
    capture = runpy.run_path(str(repo / "scripts/run_native_ppc_comparison.py"))["retain"]
    sha = scoring["sha"]
    output, dataset, binary = args.output_dir.resolve(), args.dataset_root.resolve(), args.binary.resolve()
    output.mkdir(parents=True, exist_ok=False)
    sources = list((repo / "src").rglob("*.rs")) + [repo / n for n in [
        "scripts/run_ppc_spp_seed.py", "scripts/run_ppc_rtk_float.py",
        "scripts/run_native_ppc_comparison.py", "examples/ppc_rtk_float.rs", "tests/rtk_float.rs",
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
        source_sha256={str(p.relative_to(repo)): sha(p) for p in sources},
        input_sha256=inputs, binary_sha256=sha(binary), cases=control["cases"],
        native_manifest_sha256=sha(native_archive / "predeclared-inputs.json"),
        native_report_sha256=sha(native_archive / "report.json"),
        base_ecef_m=control["base_ecef_m"],
        profiles={"previous": "kinematic", "spp": "kinematic-spp"},
        spp_configuration="SppConfig::default(), restricted to GPS 1C; fallback previous RTK/header on error",
        covariance="kinematic RTK variance 900 m2, no SPP covariance fusion; ambiguity block retained",
        changed_thresholds=False, reference_used_to_choose_configuration=False,
        native_and_rust_SPP_algorithm_equivalence_claimed=False,
        kinematic_FIX_enabled=False, timing_memory_acceptance=False, practical_PPC_parity_accepted=False,
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
        profiles = {}
        for profile, mode in manifest["profiles"].items():
            assert all(sha(repo / p) == h for p, h in manifest["source_sha256"].items())
            assert sha(binary) == manifest["binary_sha256"]
            command = [str(binary), str(root), mode, ",".join(map(str, control["base_ecef_m"][city]))]
            result = subprocess.run(command, capture_output=True, cwd=repo)
            capture(output, case + "." + profile, command, result.stdout, result.stderr, result.returncode)
            assert result.returncode == 0, result.stderr.decode()
            rows, metrics = scoring["inspect"](result.stdout)
            metrics["reference_scores"] = scoring["score"](rows, references, city)
            if profile == "previous":
                with gzip.open(previous / (case + ".rust.stdout.csv.gz"), "rt") as stream:
                    assert rows == list(csv.DictReader(stream))
                metrics["all_columns_identical_to_2972bf5_baseline"] = True
            else:
                accepted = [r for r in rows if r["status"] == "FLOAT"]
                metrics["accepted_seed_sources"] = dict(collections.Counter(r["seed_source"] for r in accepted))
                metrics["accepted_spp_fallback_errors"] = dict(collections.Counter(r["spp_failure"] for r in accepted if r["spp_failure"]))
            profiles[profile] = metrics
        profiles["native_default_SPP_seed"] = native_report["cases"][case]["native_default_SPP_seed"]
        report["cases"][case] = profiles
        save(output / "report.json", report)
        print(case, {p: m["status_counts"].get("FLOAT", 0) for p, m in profiles.items()}, flush=True)
    assert all(sha(Path(p)) == h for p, h in inputs.items())
    assert all(sha(repo / p) == h for p, h in manifest["source_sha256"].items())


if __name__ == "__main__":
    main()
