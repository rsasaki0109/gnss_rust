"""Freeze inputs and retain Rust PPC observation/product readiness audits."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import zipfile


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def inspect_rows(mode, data):
    rows = list(csv.DictReader(io.StringIO(data.decode())))
    assert len(rows) == (6 if mode == "observations" else 9)
    if mode == "observations":
        assert {r["case"] for r in rows} == {
            c + "-" + r for c in ["nagoya", "tokyo"] for r in ["run1", "run2", "run3"]
        }
        for r in rows:
            assert int(r["rover_epochs"]) == int(r["coincident_pairs"]) + int(r["rover_only"])
            assert int(r["base_epochs"]) == int(r["coincident_pairs"]) + int(r["base_only"])
            assert 0 < int(r["common_gps_l1_ge4_pairs"]) <= int(r["coincident_pairs"])
            assert int(r["gps_nav_records"]) > 0
        return {k: sum(int(r[k]) for r in rows) for k in [
            "rover_epochs", "base_epochs", "coincident_pairs", "rover_only", "base_only",
            "common_gps_l1_ge4_pairs", "excluded_nav_records",
        ]}
    assert len({r["date"] for r in rows}) == 9
    for r in rows:
        assert int(r["orbit_samples"]) > 0 and int(r["clock_samples"]) > 0
        assert int(r["dsb_records"]) > 0 and int(r["osb_records"]) == 0
        assert int(r["midpoint_states_with_clock"]) + int(r["midpoint_unavailable"]) == int(r["orbit_satellites"])
    return {k: sum(int(r[k]) for r in rows) for k in [
        "orbit_samples", "clock_samples", "dsb_records", "osb_records",
        "receiver_clock_records", "excluded_receiver_bias_rows",
        "midpoint_states_with_clock", "midpoint_unavailable",
    ]}


def verify(repo, output):
    manifest = json.loads((output / "predeclared-inputs.json").read_text())
    for name, expected in manifest["source_sha256"].items():
        assert sha(repo / name) == expected, name
    report = dict(performance_accepted=False, RTK_solutions_produced=0, PPP_solutions_produced=0)
    for mode in ["observations", "products"]:
        execution = json.loads((output / (mode + ".execution.json")).read_text())
        stdout = gzip.decompress((output / (mode + ".stdout.csv.gz")).read_bytes())
        stderr = gzip.decompress((output / (mode + ".stderr.log.gz")).read_bytes())
        assert execution["returncode"] == 0
        assert hashlib.sha256(stdout).hexdigest() == execution["stdout_sha256"]
        assert hashlib.sha256(stderr).hexdigest() == execution["stderr_sha256"]
        report[mode] = inspect_rows(mode, stdout)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dataset-root", type=Path)
    p.add_argument("--dataset-zip", type=Path)
    p.add_argument("--products-root", type=Path)
    p.add_argument("--binary", type=Path)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--verify-only", action="store_true")
    args = p.parse_args()
    repo = Path(__file__).resolve().parents[1]
    output = args.output_dir.resolve()
    if args.verify_only:
        assert verify(repo, output) == json.loads((output / "report.json").read_text())
        print("Verified retained 6-run / 9-day readiness evidence; no solver or network execution.")
        return
    if not all([args.dataset_root, args.dataset_zip, args.products_root, args.binary]):
        p.error("new runs require dataset-root, dataset-zip, products-root and binary")
    output.mkdir(parents=True, exist_ok=False)
    dataset, products, binary = args.dataset_root.resolve(), args.products_root.resolve(), args.binary.resolve()
    inputs = [dataset / c / r / f for c in ["nagoya", "tokyo"]
              for r in ["run1", "run2", "run3"] for f in ["rover.obs", "base.obs", "base.nav"]]
    inputs += [dataset / "LICENSE", dataset / "README.md"]
    acquisition = json.loads((products / "acquisition.json").read_text())
    assert acquisition["acquired_count"] == 27
    for product in acquisition["products"]:
        assert product["acquired"]
        for key, hash_key in [("filename", "compressed_sha256"), ("expanded_filename", "expanded_sha256")]:
            path = products / product[key]
            assert sha(path) == product[hash_key]
            inputs.append(path)
    input_hashes = {str(path): sha(path) for path in inputs}
    # Verify the retained raw observations/navigation against the cached ZIP.
    with zipfile.ZipFile(args.dataset_zip) as archive:
        for path in inputs[:20]:
            member = "PPC-Dataset/" + str(path.relative_to(dataset))
            with archive.open(member) as stream:
                assert hashlib.file_digest(stream, "sha256").hexdigest() == input_hashes[str(path)], member
    sources = list((repo / "src").rglob("*.rs"))
    sources += [repo / name for name in [
        "Cargo.toml", "Cargo.lock", "rust-toolchain.toml",
        "examples/ppc_input_audit.rs", "scripts/acquire_ppc_precise.py",
        "scripts/run_ppc_input_readiness.py", "tests/rinex_navigation.rs",
        "tests/precise.rs", "tests/bias.rs",
    ]]
    manifest = dict(
        source_sha256={str(path.relative_to(repo)): sha(path) for path in sources},
        input_sha256=input_hashes, binary_sha256=sha(binary),
        dataset_repository="https://github.com/taroz/PPC-Dataset",
        cached_dataset_zip_sha256=sha(args.dataset_zip),
        raw_inputs_match_cached_zip=True, constellation_selection="GPS navigation only",
        excluded_navigation="raw records retained by Rust selection API, counts in output",
        receiver_bias_rows="retained uninterpreted, never applied",
        reference_frame_and_lever_arm_evaluation=False,
        interpolation_probe="one 450-second-after-midnight epoch per product day",
        RTK_solution=False, PPP_solution=False, performance_accepted=False,
    )
    save(output / "predeclared-inputs.json", manifest)
    (output / "product-acquisition.json").write_bytes((products / "acquisition.json").read_bytes())
    for mode, root in [("observations", dataset), ("products", products)]:
        assert all(sha(repo / name) == h for name, h in manifest["source_sha256"].items())
        command = [str(binary), mode, str(root)]
        result = subprocess.run(command, cwd=repo, capture_output=True)
        (output / (mode + ".stdout.csv.gz")).write_bytes(gzip.compress(result.stdout, mtime=0))
        (output / (mode + ".stderr.log.gz")).write_bytes(gzip.compress(result.stderr, mtime=0))
        save(output / (mode + ".execution.json"), dict(
            command=command, returncode=result.returncode,
            stdout_sha256=hashlib.sha256(result.stdout).hexdigest(),
            stderr_sha256=hashlib.sha256(result.stderr).hexdigest(),
        ))
        assert result.returncode == 0, result.stderr.decode()
    assert all(sha(Path(name)) == expected for name, expected in input_hashes.items())
    report = verify(repo, output)
    save(output / "report.json", report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
