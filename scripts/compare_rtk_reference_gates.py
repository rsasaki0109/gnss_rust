"""Compare selected/all-reference local gates on four fixed synthetic inputs."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import subprocess


def sha(data):
    return hashlib.sha256(data).hexdigest()


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def perturb(text, cycles, epochs):
    """Change only G01 L1C numeric fields; preserve tracking flags and times."""
    assert "G    2 C1C L1C" in text
    lines = text.splitlines(keepends=True)
    epoch = -1
    changed = 0
    for index, line in enumerate(lines):
        if line.startswith(">"):
            epoch += 1
        if epoch in epochs and line.startswith("G01"):
            value = f"{Decimal(line[19:33]) + cycles:14.4f}"
            assert len(value) == 14
            lines[index] = line[:19] + value + line[33:]
            assert lines[index][:19] == line[:19] and lines[index][33:] == line[33:]
            changed += 1
    assert epoch == 31 and changed == len(epochs)
    return "".join(lines).encode()


def summarize_pair(selected_data, all_data):
    selected = list(csv.DictReader(io.StringIO(selected_data.decode())))
    complete = list(csv.DictReader(io.StringIO(all_data.decode())))
    assert len(selected) == len(complete) == 32
    stamps = [(int(r["week"]), Decimal(r["tow"])) for r in selected]
    assert stamps == sorted(set(stamps))
    unchanged = [
        "week", "tow", "satellites", "reference", "float_update_nis", "iterations",
        "ambiguity_resets", "candidate_ratio", "candidate_pairs",
    ]
    limits = {
        "all_reference_phase_rms_max_m": 0.01,
        "all_reference_phase_abs_max_m": 0.03,
        "all_reference_code_rms_max_m": 1.0,
        "all_reference_code_phase_abs_max_m": 0.5,
    }
    metrics = {
        "epochs": 32, "selected_fixed": 0, "all_reference_fixed": 0,
        "all_reference_audited_epochs": 0, "all_reference_phase_rejections": 0,
        "maximum_common_fixed_position_difference_m": 0.0,
        "changed_epochs": [],
    }
    for original, all_row in zip(selected, complete):
        assert all(original[k] == all_row[k] for k in unchanged)
        assert original["status"] in {"FLOAT", "FIX"}
        assert all_row["status"] in {"FLOAT", "FIX"}
        for row in [original, all_row]:
            assert all(math.isfinite(float(row[k])) for k in ["x_m", "y_m", "z_m"])
            assert (row["status"] == "FIX") == (row["fix_decision"] == "accepted")
        metrics["selected_fixed"] += original["status"] == "FIX"
        metrics["all_reference_fixed"] += all_row["status"] == "FIX"
        count = int(all_row["audited_references"])
        assert count in {0, 9}
        if count:
            metrics["all_reference_audited_epochs"] += 1
            for field in limits:
                value = float(all_row[field])
                assert math.isfinite(value) and value >= 0
            if all_row["status"] == "FIX":
                assert all(float(all_row[k]) <= v + 5e-7 for k, v in limits.items())
        else:
            assert all(not all_row[field] for field in limits)
            assert all_row["status"] == "FLOAT"
        if all_row["fix_decision"] == "phase_residuals":
            metrics["all_reference_phase_rejections"] += 1
        if original["status"] == all_row["status"] == "FIX":
            difference = math.sqrt(sum(
                (float(original[k]) - float(all_row[k])) ** 2
                for k in ["x_m", "y_m", "z_m"]
            ))
            assert difference <= 1e-6
            metrics["maximum_common_fixed_position_difference_m"] = max(
                metrics["maximum_common_fixed_position_difference_m"], difference
            )
        if original["fix_decision"] != all_row["fix_decision"]:
            metrics["changed_epochs"].append({
                "week": int(original["week"]), "tow": original["tow"],
                "selected_decision": original["fix_decision"],
                "all_reference_decision": all_row["fix_decision"],
                "all_reference_phase_rms_max_m": all_row["all_reference_phase_rms_max_m"],
            })
    return metrics


def verify_archive(repo, output):
    manifest = json.loads((output / "predeclared-inputs.json").read_text())
    for name, expected in manifest["repository_sha256"].items():
        assert sha((repo / name).read_bytes()) == expected, name
    for name, expected in manifest["generated_input_sha256"].items():
        assert sha((output / name).read_bytes()) == expected, name
    report = {
        "synthetic_only": True, "performance_accepted": False,
        "new_PPC_execution": False, "new_native_whole_execution": False, "cases": {},
    }
    for case in manifest["cases"]:
        streams = {}
        for mode in ["selected", "all"]:
            prefix = output / (case + "." + mode)
            execution = json.loads(prefix.with_suffix(prefix.suffix + ".execution.json").read_text())
            assert execution["returncode"] == 0
            stdout = gzip.decompress(prefix.with_suffix(prefix.suffix + ".csv.gz").read_bytes())
            stderr = gzip.decompress(prefix.with_suffix(prefix.suffix + ".stderr.log.gz").read_bytes())
            assert sha(stdout) == execution["stdout_sha256"]
            assert sha(stderr) == execution["stderr_sha256"]
            streams[mode] = stdout
        report["cases"][case] = summarize_pair(streams["selected"], streams["all"])
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    output = args.output_dir.resolve()
    if args.verify_only:
        actual = verify_archive(repo, output)
        assert actual == json.loads((output / "report.json").read_text())
        print("Verified 8 retained CLI runs / 256 epoch rows; no solver executed.")
        return
    if not args.binary:
        parser.error("--binary is required for new executions")
    binary = args.binary.resolve()
    output.mkdir(parents=True, exist_ok=False)
    original = (repo / "tests/fixtures/synthetic_rtk_rover.obs").read_text()
    generated = {
        "phase-bias.obs": perturb(original, Decimal("0.15"), set(range(10, 32))),
        "transient.obs": perturb(original, Decimal("0.25"), {6}),
    }
    for name, data in generated.items():
        (output / name).write_bytes(data)
    truth = (repo / "tests/fixtures/synthetic_rtk_truth.csv").read_text()
    base_position = ",".join(truth.splitlines()[1].split(",")[:3])
    # Freeze every Rust source and fixed input before the first solver run.
    paths = list((repo / "src").rglob("*.rs"))
    paths += [repo / name for name in [
        "Cargo.toml", "Cargo.lock", "rust-toolchain.toml",
        "scripts/compare_rtk_reference_gates.py",
        "tests/rtk_fix.rs", "tests/cli.rs",
        "tests/fixtures/synthetic_rtk.nav", "tests/fixtures/synthetic_rtk_truth.csv",
        "tests/fixtures/synthetic_rtk_rover.obs", "tests/fixtures/synthetic_rtk_base.obs",
        "tests/fixtures/synthetic_rtk_rover_noisy.obs", "tests/fixtures/synthetic_rtk_base_noisy.obs",
    ]]
    manifest = {
        "repository_sha256": {str(p.relative_to(repo)): sha(p.read_bytes()) for p in paths},
        "generated_input_sha256": {name: sha(data) for name, data in generated.items()},
        "binary_sha256": sha(binary.read_bytes()),
        "cases": ["clean", "noisy", "phase-bias", "transient"],
        "phase_bias": {"satellite": "G01", "cycles": "0.15", "zero_based_epochs": list(range(10, 32))},
        "transient": {"satellite": "G01", "cycles": "0.25", "zero_based_epochs": [6]},
        "thresholds_changed": False, "best_reference_selected": False,
        "truth_used_to_select_candidate": False, "synthetic_only": True,
        "new_PPC_execution": False, "new_native_whole_execution": False,
        "performance_accepted": False,
    }
    save_json(output / "predeclared-inputs.json", manifest)
    for case in manifest["cases"]:
        rover = (output / (case + ".obs")) if case in {"phase-bias", "transient"} else (
            repo / ("tests/fixtures/synthetic_rtk_rover" + ("_noisy" if case == "noisy" else "") + ".obs")
        )
        base = repo / ("tests/fixtures/synthetic_rtk_base" + ("_noisy" if case == "noisy" else "") + ".obs")
        for mode in ["selected", "all"]:
            assert all(sha((repo / name).read_bytes()) == h for name, h in manifest["repository_sha256"].items())
            assert all(sha((output / name).read_bytes()) == h for name, h in manifest["generated_input_sha256"].items())
            assert sha(binary.read_bytes()) == manifest["binary_sha256"]
            command = [str(binary), "rtk", "--rover", str(rover), "--base", str(base),
                       "--nav", str(repo / "tests/fixtures/synthetic_rtk.nav"),
                       "--base-position", base_position]
            if mode == "all":
                command.append("--all-reference-gates")
            result = subprocess.run(command, cwd=repo, capture_output=True)
            prefix = output / (case + "." + mode)
            prefix.with_suffix(prefix.suffix + ".csv.gz").write_bytes(gzip.compress(result.stdout, mtime=0))
            prefix.with_suffix(prefix.suffix + ".stderr.log.gz").write_bytes(gzip.compress(result.stderr, mtime=0))
            save_json(prefix.with_suffix(prefix.suffix + ".execution.json"), {
                "command": command, "returncode": result.returncode,
                "stdout_sha256": sha(result.stdout), "stderr_sha256": sha(result.stderr),
            })
            assert result.returncode == 0, result.stderr.decode()
    report = verify_archive(repo, output)
    save_json(output / "report.json", report)
    for case, result in report["cases"].items():
        print(case, json.dumps({k: v for k, v in result.items() if k != "changed_epochs"}))


if __name__ == "__main__":
    main()
