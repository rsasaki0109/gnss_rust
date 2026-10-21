"""Verify the isolated native seed comparison and optionally rescore references."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
from pathlib import Path
import runpy


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def diagnostics(archive, dataset):
    repo = archive.parents[2]
    previous = archive.parent / "2026-10-21-native-rtk-common-inputs"
    score = runpy.run_path(str(repo / "scripts/run_ppc_rtk_float.py"))["score"]
    report = json.loads((archive / "report.json").read_text())
    result = dict(post_run_only=True, acceptance_thresholds_defined=False, cases={})
    for case in report["cases"]:
        city, run = case.split("-")
        references = {}
        with (dataset / city / run / "reference.csv").open() as stream:
            for raw in csv.DictReader(stream):
                r = {k.strip(): v.strip() for k, v in raw.items()}
                references[(int(r["GPS Week"]), Decimal(r["GPS TOW (s)"]))] = r
        key = lambda r: (int(r["week"]), float(r["tow"]).hex())
        rows = {}
        for name, path in [
            ("native_default", archive / (case + ".native-default.stdout.csv.gz")),
            ("native_previous", previous / (case + ".native.stdout.csv.gz")),
            ("rust", previous / (case + ".rust.stdout.csv.gz")),
        ]:
            with gzip.open(path, "rt") as stream:
                rows[name] = {key(r): r for r in csv.DictReader(stream) if r["status"] in {"FLOAT", "FIX"}}
        common = set.intersection(*(set(r) for r in rows.values()))
        result["cases"][case] = dict(
            common_all_three_RTK_epochs=len(common),
            common_all_three_reference_scores={
                name: score([{**r[t], "status": "FLOAT"} for t in sorted(common)], references, city)
                for name, r in rows.items()
            },
            default_status_reference_scores={
                status: score([{**r, "status": "FLOAT"} for r in rows["native_default"].values()
                               if r["status"] == status], references, city)
                for status in ["FLOAT", "FIX"]
            },
        )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root", type=Path)
    parser.add_argument("--record-diagnostics", action="store_true")
    args = parser.parse_args()
    archive = Path(__file__).resolve().parent
    repo = archive.parents[2]
    previous = archive.parent / "2026-10-21-native-rtk-common-inputs"
    runner = runpy.run_path(str(repo / "scripts/run_native_ppc_comparison.py"))
    scoring = runpy.run_path(str(repo / "scripts/run_ppc_rtk_float.py"))
    manifest = json.loads((archive / "predeclared-inputs.json").read_text())
    report = json.loads((archive / "report.json").read_text())
    prior = json.loads((previous / "report.json").read_text())
    assert not report["practical_PPC_parity_accepted"]
    assert set(report["cases"]) == set(manifest["cases"])
    assert sha(previous / "predeclared-inputs.json") == manifest["prior_manifest_sha256"]
    for name, expected in manifest["source_sha256"].items():
        assert sha(repo / name) == expected, name
    original = (repo / "scripts/native_ppc_rtk.cpp").read_text()
    variant = original.replace("config.prefer_rover_position_seed = true;", "config.prefer_rover_position_seed = false;")
    variant = variant.replace("previous_own_RTK_seed=true AR_default_continuous=true", "SPP_seed_default=true previous_RTK_fallback=true AR_default_continuous=true")
    assert (archive / "native_ppc_rtk_default_seed.cpp").read_text() == variant
    assert sha(archive / "native_ppc_rtk_default_seed.cpp") == manifest["variant_sha256"]
    build = json.loads((archive / "adapter-build.json").read_text())
    assert build["returncode"] == 0 and build["source_sha256"] == manifest["variant_sha256"]
    for name in ["stdout", "stderr"]:
        assert hashlib.sha256(gzip.decompress((archive / ("adapter-build." + name + ".log.gz")).read_bytes())).hexdigest() == build[name + "_sha256"]
    total = 0
    for case, c in report["cases"].items():
        execution = json.loads((archive / (case + ".native-default.execution.json")).read_text())
        data = gzip.decompress((archive / (case + ".native-default.stdout.csv.gz")).read_bytes())
        stderr = gzip.decompress((archive / (case + ".native-default.stderr.log.gz")).read_bytes())
        assert execution["returncode"] == 0
        assert hashlib.sha256(data).hexdigest() == execution["csv_sha256"]
        assert hashlib.sha256(stderr).hexdigest() == execution["stderr_sha256"]
        rows, scored, metrics = runner["native_summary"](data)
        expected = c["native_default_SPP_seed"]
        assert metrics == {k: v for k, v in expected.items() if k != "reference_scores"}
        assert metrics["valid_outputs_with_mismatched_time"] == 0 and b"input_drained=1" in stderr
        with gzip.open(previous / (case + ".native.stdout.csv.gz"), "rt") as stream:
            old = list(csv.DictReader(stream))
        keys = lambda r: [(x["kind"], int(x["week"]), float(x["tow"]).hex(), x["side"]) for x in r]
        assert keys(rows) == keys(old)
        assert c["unchanged_input_ledger_sha256"] == prior["cases"][case]["input_comparison"]["native_ledger_sha256"]
        assert c["native_previous_RTK_seed"] == prior["cases"][case]["native"]
        assert c["rust_previous_RTK_seed"] == prior["cases"][case]["rust"]
        if args.dataset_root:
            city, run = case.split("-")
            root = args.dataset_root / city / run
            references = {}
            for name in ["rover.obs", "base.obs", "base.nav", "reference.csv"]:
                k = next(k for k in manifest["input_sha256"] if k.endswith(f"/{city}/{run}/{name}"))
                assert sha(root / name) == manifest["input_sha256"][k]
            with (root / "reference.csv").open() as stream:
                for raw in csv.DictReader(stream):
                    r = {k.strip(): v.strip() for k, v in raw.items()}
                    references[(int(r["GPS Week"]), Decimal(r["GPS TOW (s)"]))] = r
            assert scoring["score"](scored, references, city) == expected["reference_scores"]
        total += metrics["valid_RTK_positions"]
    if args.record_diagnostics:
        assert args.dataset_root
        path = archive / "diagnostics.json"
        assert not path.exists()
        path.write_text(json.dumps(diagnostics(archive, args.dataset_root), indent=2) + "\n")
    elif args.dataset_root:
        assert diagnostics(archive, args.dataset_root) == json.loads((archive / "diagnostics.json").read_text())
    if not args.record_diagnostics:
        for name, expected in json.loads((archive / "evidence-sha256.json").read_text()).items():
            assert sha(archive / name) == expected, name
    print(f"Verified six isolated native default-SPP-seed executions; valid RTK={total}/11733.")


if __name__ == "__main__":
    main()
