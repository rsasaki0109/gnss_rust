"""Run three fixed synthetic cases and verify equivalent integer candidates."""
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
import subprocess


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def verify_rows(data, truth):
    rows = list(csv.DictReader(io.StringIO(data.decode())))
    groups = collections.OrderedDict()
    for row in rows:
        groups.setdefault((int(row["week"]), Decimal(row["tow"])), []).append(row)
    assert len(groups) == 32 and len(rows) == 288
    assert list(groups) == sorted(groups)
    sensitive = correct = original_pass = 0
    maximum_spread = maximum_covariance = maximum_score = 0.0
    decisions = []
    for stamp, entries in groups.items():
        assert {r["reference"] for r in entries} == set(truth)
        assert len(entries) == len(truth)
        source = entries[0]["source_reference"]
        assert all(r["source_reference"] == source for r in entries)
        original = next(r for r in entries if r["reference"] == source)
        values = {source: 0, **dict(zip(original["targets"].split(";"),
                                       map(int, original["integers"].split(";"))))}
        changed = len({r["local_rejection"] for r in entries}) > 1
        assert all(r["reference_sensitive_local_rejection"] == str(changed).lower()
                   for r in entries)
        sensitive += changed
        original_pass += original["local_rejection"] == "pass"
        correct += original["candidate_correct_for_synthetic_truth"] == "true"
        for row in entries:
            reference = row["reference"]
            targets = row["targets"].split(";")
            integers = list(map(int, row["integers"].split(";")))
            assert len(targets) == len(integers) == len(truth)-1
            assert set(targets) == set(truth)-{reference}
            assert integers == [values[s]-values[reference] for s in targets]
            truth_correct = integers == [truth[reference]-truth[s] for s in targets]
            assert row["candidate_correct_for_synthetic_truth"] == str(truth_correct).lower()
            spread = math.sqrt(sum((float(row[k])-float(original[k]))**2
                                  for k in ["baseline_x_m", "baseline_y_m", "baseline_z_m"]))
            assert spread <= 1e-6
            maximum_spread = max(maximum_spread, spread)
            q = [[float(row[k]) for k in keys]
                 for keys in [["qxx", "qxy", "qxz"], ["qxy", "qyy", "qyz"],
                              ["qxz", "qyz", "qzz"]]]
            l = [[0.0]*3 for _ in range(3)]
            for i in range(3):
                for j in range(i+1):
                    value = q[i][j]-sum(l[i][k]*l[j][k] for k in range(j))
                    assert math.isfinite(value)
                    if i == j:
                        assert value > 0
                        l[i][j] = math.sqrt(value)
                    else:
                        l[i][j] = value/l[j][j]
            for field in ["qxx", "qxy", "qxz", "qyy", "qyz", "qzz"]:
                difference = abs(float(row[field])-float(original[field]))
                assert difference <= 1e-10*max(1.0, abs(float(original[field])))
                maximum_covariance = max(maximum_covariance, difference)
            for field in ["ambiguity_cost", "nominal_postfit_score_per_row"]:
                actual = float(row[field])
                expected = float(original[field])
                assert math.isfinite(actual) and actual >= 0
                assert abs(actual-expected) <= 1e-6*max(1.0, abs(expected))
                maximum_score = max(maximum_score, abs(actual-expected))
            for field in ["phase_rms_m", "phase_max_m", "code_rms_m",
                          "code_phase_rms_m", "code_phase_max_m"]:
                assert math.isfinite(float(row[field])) and float(row[field]) >= 0
        decisions.append(dict(week=stamp[0], tow=str(stamp[1]), original_reference=source,
                              original_local_rejection=original["local_rejection"],
                              candidate_correct_for_synthetic_truth=original["candidate_correct_for_synthetic_truth"],
                              reference_sensitive=changed,
                              reference_rejections={r["reference"]: r["local_rejection"] for r in entries}))
    return dict(epochs=len(groups), reference_rows=len(rows),
                reference_sensitive_epochs=sensitive,
                synthetically_correct_candidate_epochs=correct,
                original_basis_local_pass_epochs=original_pass,
                maximum_conditioned_position_spread_m=maximum_spread,
                maximum_conditioned_covariance_difference_m2=maximum_covariance,
                maximum_nominal_score_difference=maximum_score), decisions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    binary = args.binary.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    files = ["src/rtk.rs", "src/rtk_fix.rs", "src/filter.rs", "src/lambda.rs",
             "examples/rtk_candidate_reference_audit.rs",
             "scripts/run_rtk_candidate_reference_audit.py",
             "tests/fixtures/synthetic_rtk.nav", "tests/fixtures/synthetic_rtk_truth.csv"]
    files += [f"tests/fixtures/synthetic_rtk_{receiver}{suffix}.obs"
              for receiver in ["rover", "base"] for suffix in ["", "_noisy"]]
    frozen = {name: sha(repo/name) for name in files}
    manifest = dict(source_input_sha256=frozen, binary_sha256=sha(binary),
                    cases=["clean", "noisy", "phase-bias"],
                    fixed_phase_bias_cycles=0.15, first_zero_based_fault_epoch=10,
                    candidate_source="one original-basis LAMBDA candidate, re-expressed exactly",
                    best_reference_selected=False, truth_used_to_select_candidate=False,
                    truth_used_to_report_synthetic_correctness_only=True,
                    synthetic_only=True, new_PPC_execution=False,
                    new_native_whole_execution=False, performance_accepted=False)
    (output/"predeclared-inputs.json").write_text(json.dumps(manifest, indent=2)+"\n")
    truth = {line.split(",")[0]: int(line.split(",")[1])
             for line in (repo/"tests/fixtures/synthetic_rtk_truth.csv").read_text().splitlines()[2:]}
    report = dict(cases={}, synthetic_only=True, performance_accepted=False)
    for mode in manifest["cases"]:
        assert all(sha(repo/name) == expected for name, expected in frozen.items())
        assert sha(binary) == manifest["binary_sha256"]
        command = [str(binary), mode]
        result = subprocess.run(command, cwd=repo, capture_output=True)
        (output/(mode+".csv.gz")).write_bytes(gzip.compress(result.stdout, mtime=0))
        (output/(mode+".stderr.log.gz")).write_bytes(gzip.compress(result.stderr, mtime=0))
        execution = dict(command=command, returncode=result.returncode,
                         stdout_sha256=hashlib.sha256(result.stdout).hexdigest(),
                         stderr_sha256=hashlib.sha256(result.stderr).hexdigest())
        (output/(mode+".execution.json")).write_text(json.dumps(execution, indent=2)+"\n")
        assert result.returncode == 0, result.stderr.decode()
        metrics, decisions = verify_rows(result.stdout, truth)
        report["cases"][mode] = metrics
        (output/(mode+".decisions.json")).write_text(json.dumps(decisions, indent=2)+"\n")
        (output/"report.json").write_text(json.dumps(report, indent=2)+"\n")
        print(mode, json.dumps(metrics), flush=True)
    assert all(sha(repo/name) == expected for name, expected in frozen.items())


if __name__ == "__main__":
    main()
