"""Isolate pinned native RTK's default SPP seed from the previous-RTK seed."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
from pathlib import Path
import runpy
import subprocess
import tempfile


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root", type=Path, required=True)
    parser.add_argument("--native-build", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    previous = repo / "docs/comparison/2026-10-21-native-rtk-common-inputs"
    source = repo / "scripts/native_ppc_rtk.cpp"
    output, dataset, build = args.output_dir.resolve(), args.dataset_root.resolve(), args.native_build.resolve()
    output.mkdir(parents=True, exist_ok=False)
    runner = runpy.run_path(str(repo / "scripts/run_native_ppc_comparison.py"))
    scoring = runpy.run_path(str(repo / "scripts/run_ppc_rtk_float.py"))
    original = json.loads((previous / "predeclared-inputs.json").read_text())
    prior = json.loads((previous / "report.json").read_text())
    provenance = json.loads((previous / "native-build.json").read_text())
    sha = runner["sha"]
    for name, expected in original["source_sha256"].items():
        assert sha(repo / name) == expected, name
    for name, expected in original["native_archive_sha256"].items():
        assert sha(Path(name)) == expected, name
    text = source.read_text()
    changes = {
        "config.prefer_rover_position_seed = true;": "config.prefer_rover_position_seed = false;",
        "previous_own_RTK_seed=true AR_default_continuous=true":
            "SPP_seed_default=true previous_RTK_fallback=true AR_default_continuous=true",
    }
    for old, new in changes.items():
        assert text.count(old) == 1
        text = text.replace(old, new)
    variant = output / "native_ppc_rtk_default_seed.cpp"
    variant.write_text(text)
    binary = build / "native_ppc_rtk_default_seed"
    command = provenance["commands"][-1][:]
    command[command.index("scripts/native_ppc_rtk.cpp")] = str(variant)
    command[-1] = str(binary)
    compile_result = subprocess.run(command, capture_output=True, cwd=repo)
    (output / "adapter-build.stdout.log.gz").write_bytes(gzip.compress(compile_result.stdout, mtime=0))
    (output / "adapter-build.stderr.log.gz").write_bytes(gzip.compress(compile_result.stderr, mtime=0))
    save(output / "adapter-build.json", dict(
        command=command, returncode=compile_result.returncode,
        stdout_sha256=hashlib.sha256(compile_result.stdout).hexdigest(),
        stderr_sha256=hashlib.sha256(compile_result.stderr).hexdigest(),
        source_sha256=sha(variant), upstream_library_reused_unmodified=True,
    ))
    assert compile_result.returncode == 0, compile_result.stderr.decode()
    inputs = {}
    for case in original["cases"]:
        city, run = case.split("-")
        for name in ["rover.obs", "base.obs", "base.nav", "reference.csv"]:
            path = dataset / city / run / name
            key = next(k for k in original["input_sha256"] if k.endswith(f"/{city}/{run}/{name}"))
            assert sha(path) == original["input_sha256"][key]
            inputs[str(path)] = sha(path)
    sources = {**original["source_sha256"], str(Path(__file__).resolve().relative_to(repo)): sha(Path(__file__))}
    manifest = dict(
        source_sha256=sources, variant_sha256=sha(variant), binary_sha256=sha(binary),
        native_archive_sha256=original["native_archive_sha256"], input_sha256=inputs,
        prior_manifest_sha256=sha(previous / "predeclared-inputs.json"), cases=original["cases"],
        upstream_commit=original["upstream_commit"], base_ecef_m=original["base_ecef_m"],
        sole_solver_configuration_change="prefer_rover_position_seed true -> false (native default)",
        unchanged_wrapper_seed="initial header then previous own valid RTK, used only as native fallback when SPP fails",
        unchanged_native_QC_AR_continuity_and_GPS_1C_selection=True,
        reference_used_to_choose_configuration=False, practical_PPC_parity_accepted=False,
        timing_memory_acceptance=False,
    )
    save(output / "predeclared-inputs.json", manifest)
    report = dict(practical_PPC_parity_accepted=False, cases={})
    with tempfile.TemporaryDirectory(prefix="native-default-seed-") as temp:
        temp = Path(temp)
        for case in original["cases"]:
            city, run = case.split("-")
            for path in temp.iterdir():
                path.unlink()
            native_csv, ledger = temp / "native.csv", temp / "ledger.csv"
            command = [str(binary), str(dataset / city / run), *map(str, original["base_ecef_m"][city]),
                       str(native_csv), str(ledger)]
            result = subprocess.run(command, capture_output=True, cwd=repo)
            data = native_csv.read_bytes() if native_csv.exists() else b""
            runner["retain"](output, case + ".native-default", command, data, result.stderr, result.returncode,
                             csv_capture="explicit native output file", process_stdout=result.stdout)
            assert result.returncode == 0, result.stderr.decode()
            # Before processing, the wrapper records original measurements and
            # pending-reset notices; changing the seed cannot change these.
            expected_ledger = prior["cases"][case]["input_comparison"]["native_ledger_sha256"]
            actual_ledger = sha(ledger)
            if actual_ledger != expected_ledger:
                (output / (case + ".failed-input-ledger.csv.gz")).write_bytes(gzip.compress(ledger.read_bytes(), mtime=0))
                raise AssertionError((case, "input ledger changed", actual_ledger, expected_ledger))
            rows, rtk, metrics = runner["native_summary"](data)
            assert metrics["valid_outputs_with_mismatched_time"] == 0
            assert b"input_drained=1" in result.stderr
            with gzip.open(previous / (case + ".native.stdout.csv.gz"), "rt") as stream:
                old_rows = list(csv.DictReader(stream))
            keys = lambda r: [(x["kind"], int(x["week"]), float(x["tow"]).hex(), x["side"]) for x in r]
            assert keys(rows) == keys(old_rows)
            references = {}
            with (dataset / city / run / "reference.csv").open() as stream:
                for raw in csv.DictReader(stream):
                    r = {k.strip(): v.strip() for k, v in raw.items()}
                    references[(int(r["GPS Week"]), Decimal(r["GPS TOW (s)"]))] = r
            metrics["reference_scores"] = scoring["score"](rtk, references, city)
            report["cases"][case] = dict(
                native_default_SPP_seed=metrics,
                native_previous_RTK_seed=prior["cases"][case]["native"],
                rust_previous_RTK_seed=prior["cases"][case]["rust"],
                unchanged_input_ledger_sha256=actual_ledger,
            )
            save(output / "report.json", report)
            print(case, metrics["status_counts"], flush=True)
    assert all(sha(Path(p)) == expected for p, expected in inputs.items())
    assert all(sha(repo / p) == expected for p, expected in sources.items())
    assert sha(binary) == manifest["binary_sha256"]


if __name__ == "__main__":
    main()
