"""Replay unchanged binaries to retain and independently verify input ledgers."""
import gzip
import hashlib
import json
from pathlib import Path
import runpy
import subprocess
import tempfile


def main():
    archive = Path(__file__).resolve().parent
    repo = archive.parents[2]
    runner = runpy.run_path(str(repo / "scripts/run_native_ppc_comparison.py"))
    manifest = json.loads((archive / "predeclared-inputs.json").read_text())
    source_hashes = manifest["source_sha256"]
    for name, expected in source_hashes.items():
        assert runner["sha"](repo / name) == expected
    for name, expected in manifest["binary_sha256"].items():
        assert runner["sha"](Path(name)) == expected
    with tempfile.TemporaryDirectory(prefix="ppc-retain-ledgers-") as temporary:
        temporary = Path(temporary)
        for case in manifest["cases"]:
            record = json.loads((archive / (case + ".native.execution.json")).read_text())
            native_command = record["command"][:]
            native_command[-2:] = [str(temporary / "native.csv"), str(temporary / "native-inputs.csv")]
            # The native adapter deliberately refuses to overwrite prior output.
            for path in temporary.iterdir():
                path.unlink()
            result = subprocess.run(native_command, capture_output=True, cwd=repo)
            assert result.returncode == 0, result.stderr.decode()
            assert runner["sha"](temporary / "native.csv") == record["csv_sha256"]
            rust_binary = next(k for k in manifest["binary_sha256"] if k.endswith("/ppc_input_ledger"))
            rust_command = [rust_binary, native_command[1]]
            rust_inputs = temporary / "rust-inputs.csv"
            with rust_inputs.open("wb") as destination:
                rust = subprocess.run(rust_command, stdout=destination, stderr=subprocess.PIPE, cwd=repo)
            assert rust.returncode == 0, rust.stderr.decode()
            comparison = runner["compare_ledgers"](temporary / "native-inputs.csv", rust_inputs)
            assert comparison == json.loads((archive / (case + ".input-comparison.json")).read_text())
            for name in ["native", "rust"]:
                (archive / (case + "." + name + ".input-ledger.csv.gz")).write_bytes(
                    gzip.compress((temporary / (name + "-inputs.csv")).read_bytes(), mtime=0)
                )
            replay_path = archive / (case + ".input-retention-replay.json")
            assert not replay_path.exists()
            replay_path.write_text(json.dumps(dict(
                native_command=native_command, rust_ledger_command=rust_command,
                native_returncode=result.returncode, rust_ledger_returncode=rust.returncode,
                native_csv_sha256=record["csv_sha256"], native_output_identical_to_primary=True,
                input_ledgers_identical_to_primary=True,
                native_stderr_sha256=hashlib.sha256(result.stderr).hexdigest(),
            ), indent=2) + "\n")
            (archive / (case + ".input-retention-replay.stderr.log.gz")).write_bytes(
                gzip.compress(result.stderr, mtime=0)
            )
            print(case, "retained", comparison["rows"], "identical input rows", flush=True)
    for name, expected in source_hashes.items():
        assert runner["sha"](repo / name) == expected


if __name__ == "__main__":
    main()
