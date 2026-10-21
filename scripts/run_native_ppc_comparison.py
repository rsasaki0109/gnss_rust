"""Compare whole pinned native RTK and Rust FLOAT on decoded-equal GPS 1C inputs."""
import argparse
import collections
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
import tempfile


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def save(path, data):
    path.write_text(json.dumps(data, indent=2) + "\n")


def ledger_key(row):
    # The native Observation API has no missing SNR/LLI presence bits.
    # Normalize only those two missing fields to its documented zero default.
    values = [row["side"], int(row["week"]), float(row["tow"]).hex(), row["satellite"]]
    for present, value in [("has_code", "code_m"), ("has_phase", "phase_cycles"),
                           ("has_doppler", "doppler_hz")]:
        flag = int(row[present])
        assert flag in {0, 1}
        values.extend([flag, float(row[value]).hex() if flag else None])
    values.extend([float(row["snr_db_hz"]).hex(), int(row["lli"]), int(row["derived_reset"])])
    return values


def compare_ledgers(native, rust):
    digest = hashlib.sha256()
    count = 0
    from itertools import zip_longest
    with native.open() as a, rust.open() as b:
        for count, (left, right) in enumerate(zip_longest(csv.DictReader(a), csv.DictReader(b)), 1):
            assert left is not None and right is not None, (count, "input ledger length mismatch")
            assert ledger_key(left) == ledger_key(right), (count, left, right)
            digest.update((json.dumps(ledger_key(left), separators=(",", ":"))+"\n").encode())
    assert count > 0
    return dict(rows=count, normalized_sha256=digest.hexdigest(),
                native_ledger_sha256=sha(native), rust_ledger_sha256=sha(rust),
                code_phase_doppler_presence_preserved=True,
                missing_SNR_LLI_normalized_to_native_zero=True,
                all_decoded_values_and_derived_resets_equal=True)


def retain(output, stem, command, stdout, stderr, returncode, csv_capture="process stdout", process_stdout=None):
    (output/(stem+".stdout.csv.gz")).write_bytes(gzip.compress(stdout, mtime=0))
    (output/(stem+".stderr.log.gz")).write_bytes(gzip.compress(stderr, mtime=0))
    save(output/(stem+".execution.json"), dict(
        command=command, returncode=returncode,
        csv_sha256=hashlib.sha256(stdout).hexdigest(), csv_capture=csv_capture,
        stderr_sha256=hashlib.sha256(stderr).hexdigest(),
    ))
    if process_stdout is not None:
        (output/(stem+".process-stdout.log.gz")).write_bytes(gzip.compress(process_stdout, mtime=0))


def native_summary(data):
    rows = list(csv.DictReader(io.StringIO(data.decode())))
    assert rows and all(None not in row for row in rows)
    pairs = [r for r in rows if r["kind"] == "pair"]
    counts = collections.Counter(r["status"] for r in pairs)
    scored = []
    mismatched_time = 0
    for row in pairs:
        valid = row["is_valid"] == "1"
        assert valid == (row["status"] != "INVALID")
        if valid:
            assert all(math.isfinite(float(row[k])) for k in ["x_m", "y_m", "z_m"])
            difference = (int(row["solution_week"])-int(row["week"]))*604800 + (
                float(row["solution_tow"])-float(row["tow"])
            )
            if abs(difference) > 1e-7:
                mismatched_time += 1
                continue
            if row["status"] in {"FLOAT", "FIX"}:
                scored.append({**row, "status": "FLOAT"})
        else:
            assert all(not row[k] for k in ["x_m", "y_m", "z_m"])
    result = dict(
        output_rows=len(rows), paired_epochs=len(pairs), status_counts=dict(counts),
        rover_only=sum(r["side"] == "rover" for r in rows),
        base_only=sum(r["side"] == "base" for r in rows),
        valid_RTK_positions=len(scored), valid_outputs_with_mismatched_time=mismatched_time,
        fallback_positions=sum(r["status"] in {"SPP", "DGPS", "PROPAGATED"} for r in pairs),
    )
    return rows, scored, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root", type=Path, required=True)
    parser.add_argument("--native-binary", type=Path, required=True)
    parser.add_argument("--native-build", type=Path, required=True)
    parser.add_argument("--rust-binary", type=Path, required=True)
    parser.add_argument("--rust-ledger-binary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    dataset, output = args.dataset_root.resolve(), args.output_dir.resolve()
    native, rust, rust_ledger = args.native_binary.resolve(), args.rust_binary.resolve(), args.rust_ledger_binary.resolve()
    output.mkdir(parents=True, exist_ok=False)
    previous = repo/"docs/comparison/2026-10-21-ppc-geometry-gate"
    control = json.loads((previous/"predeclared-inputs.json").read_text())
    helpers = runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))
    sources = list((repo/"src").rglob("*.rs")) + [repo/name for name in [
        "scripts/native_ppc_rtk.cpp", "scripts/run_native_ppc_comparison.py",
        "scripts/run_ppc_rtk_float.py", "examples/ppc_rtk_float.rs", "examples/ppc_input_ledger.rs",
        "Cargo.toml", "Cargo.lock", "rust-toolchain.toml",
    ]]
    paths = [dataset/case.split("-")[0]/case.split("-")[1]/name for case in control["cases"]
             for name in ["rover.obs", "base.obs", "base.nav", "reference.csv"]]
    build = args.native_build.resolve()
    archives = [build/name for name in [
        "libgnss_rtk_runtime.a", "libgnss_lib.a",
        "third_party/ginan-iers2010/libginan_iers2010.a", "third_party/sofa/libsofa.a",
    ]]
    manifest = dict(
        cases=control["cases"], upstream_commit="72f3b7c2c4c088dc575286a6f0edf4e407bc499a",
        source_sha256={str(p.relative_to(repo)): sha(p) for p in sources},
        input_sha256={str(p): sha(p) for p in paths},
        binary_sha256={str(p): sha(p) for p in [native, rust, rust_ledger]},
        native_archive_sha256={str(p): sha(p) for p in archives},
        native_whole_RTKProcessor=True, GPS_1C_common_epochs_only=True,
        native_AR="default CONTINUOUS/EXTENDED enabled; Rust kinematic FIX not implemented",
        base_ecef_m=control["base_ecef_m"], elevation_degrees=10.0, maximum_baseline_m=10000.0,
        initial_seed="same rover RINEX header", subsequent_seed="each solver's previous own valid RTK output",
        continuity_mapping="raw matched LLI retained in ledgers; skipped-epoch notices OR LLI bit 1 into native copies",
        native_solver_defaults_otherwise_preserved=True, truth_used_by_solver=False,
        native_SPP_fallback_excluded_from_RTK_accuracy_but_counted=True,
        timing_memory_acceptance=False, practical_PPC_parity_accepted=False,
    )
    save(output/"predeclared-inputs.json", manifest)
    report = dict(native_whole_RTKProcessor=True, practical_PPC_parity_accepted=False, cases={})
    with tempfile.TemporaryDirectory(prefix="ppc-native-inputs-") as temp:
        temp = Path(temp)
        for case in manifest["cases"]:
            city, run = case.split("-")
            directory = dataset/city/run
            references = {}
            with (directory/"reference.csv").open() as stream:
                for raw in csv.DictReader(stream):
                    row = {k.strip(): v.strip() for k, v in raw.items()}
                    references[(int(row["GPS Week"]), Decimal(row["GPS TOW (s)"]))] = row
            for name, expected in manifest["source_sha256"].items():
                assert sha(repo/name) == expected
            native_output, native_inputs, rust_inputs = [temp/(case+suffix) for suffix in [".native.csv", ".native-inputs.csv", ".rust-inputs.csv"]]
            native_command = [str(native), str(directory), *map(str, manifest["base_ecef_m"][city]),
                              str(native_output), str(native_inputs)]
            result = subprocess.run(native_command, capture_output=True, cwd=repo)
            native_data = native_output.read_bytes() if native_output.exists() else b""
            retain(output, case+".native", native_command, native_data, result.stderr, result.returncode,
                   csv_capture="explicit native output file", process_stdout=result.stdout)
            assert result.returncode == 0, result.stderr.decode()
            with rust_inputs.open("wb") as destination:
                ledger_result = subprocess.run([str(rust_ledger), str(directory)], stdout=destination, stderr=subprocess.PIPE, cwd=repo)
            assert ledger_result.returncode == 0, ledger_result.stderr.decode()
            try:
                input_comparison = compare_ledgers(native_inputs, rust_inputs)
            except AssertionError as error:
                for name, path in [("native", native_inputs), ("rust", rust_inputs)]:
                    (output/(case+"."+name+".failed-input-ledger.csv.gz")).write_bytes(
                        gzip.compress(path.read_bytes(), mtime=0)
                    )
                save(output/(case+".input-comparison-failed.json"), dict(
                    all_decoded_values_and_derived_resets_equal=False, error=str(error),
                ))
                raise
            save(output/(case+".input-comparison.json"), input_comparison)
            rust_command = [str(rust), str(directory), "kinematic", ",".join(map(str, manifest["base_ecef_m"][city]))]
            rust_result = subprocess.run(rust_command, capture_output=True, cwd=repo)
            retain(output, case+".rust", rust_command, rust_result.stdout, rust_result.stderr, rust_result.returncode)
            assert rust_result.returncode == 0, rust_result.stderr.decode()
            native_rows, native_rtk, native_metrics = native_summary(native_data)
            rust_rows, rust_metrics = helpers["inspect"](rust_result.stdout)
            # Both emit lossless representations of the same binary f64.
            # C++ scientific precision and Rust fixed precision have different
            # decimal expansions at fractional unmatched epochs.
            keys = lambda rows: [(r["kind"], int(r["week"]), float(r["tow"]).hex(), r["side"]) for r in rows]
            assert keys(native_rows) == keys(rust_rows)
            assert native_metrics["paired_epochs"] == rust_metrics["paired_epochs"]
            assert "input_drained=1" in result.stderr.decode()
            assert native_metrics["valid_outputs_with_mismatched_time"] == 0
            with gzip.open(previous/(case+".kinematic.csv.gz"), "rt") as stream:
                old = list(csv.DictReader(stream))
            assert rust_rows == old
            native_metrics["reference_scores"] = helpers["score"](native_rtk, references, city)
            rust_metrics["reference_scores"] = helpers["score"](rust_rows, references, city)
            report["cases"][case] = dict(native=native_metrics, rust=rust_metrics, input_comparison=input_comparison)
            save(output/"report.json", report)
            print(case, json.dumps({
                "native_statuses": native_metrics["status_counts"], "native_RTK_positions": native_metrics["valid_RTK_positions"],
                "rust_FLOAT_positions": rust_metrics["status_counts"].get("FLOAT", 0), "input_rows_compared": input_comparison["rows"],
            }), flush=True)
    assert all(sha(Path(name)) == h for name, h in manifest["input_sha256"].items())
    assert all(sha(repo/name) == h for name, h in manifest["source_sha256"].items())


if __name__ == "__main__":
    main()
