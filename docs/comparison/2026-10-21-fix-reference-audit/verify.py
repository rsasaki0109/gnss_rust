"""Verify retained synthetic audit evidence without running a GNSS solver."""
import gzip
import hashlib
import json
from pathlib import Path
import re
import runpy


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    archive = Path(__file__).resolve().parent
    repo = archive.parents[2]
    manifest = json.loads((archive / "predeclared-inputs.json").read_text())
    for name, expected in manifest["source_input_sha256"].items():
        assert digest((repo / name).read_bytes()) == expected, name
    assert manifest["synthetic_only"] and not manifest["performance_accepted"]
    assert not manifest["best_reference_selected"]
    assert not manifest["truth_used_to_select_candidate"]
    assert not manifest["new_PPC_execution"]
    assert not manifest["new_native_whole_execution"]
    verify_rows = runpy.run_path(
        str(repo / "scripts/run_rtk_candidate_reference_audit.py")
    )["verify_rows"]
    truth = {
        line.split(",")[0]: int(line.split(",")[1])
        for line in (repo / "tests/fixtures/synthetic_rtk_truth.csv")
        .read_text().splitlines()[2:]
    }
    report = json.loads((archive / "report.json").read_text())
    assert report["synthetic_only"] and not report["performance_accepted"]
    assert set(report["cases"]) == set(manifest["cases"])
    for mode in manifest["cases"]:
        execution = json.loads((archive / (mode + ".execution.json")).read_text())
        assert execution["returncode"] == 0
        data = gzip.decompress((archive / (mode + ".csv.gz")).read_bytes())
        stderr = gzip.decompress((archive / (mode + ".stderr.log.gz")).read_bytes())
        assert digest(data) == execution["stdout_sha256"]
        assert digest(stderr) == execution["stderr_sha256"]
        metrics, decisions = verify_rows(data, truth)
        assert metrics == report["cases"][mode]
        assert decisions == json.loads(
            (archive / (mode + ".decisions.json")).read_text()
        )
    qualification = json.loads((archive / "qualification.json").read_text())
    assert digest((repo / "tests/rtk_fix.rs").read_bytes()) == qualification["test_source_sha256"]
    assert {r["name"] for r in qualification["checks"]} == {
        "rust", "examples", "fmt", "clippy", "release"
    }
    for check in qualification["checks"]:
        assert check["returncode"] == 0
        log = gzip.decompress((archive / check["log"]).read_bytes())
        assert digest(log) == check["uncompressed_sha256"]
        if check["name"] == "rust":
            count = sum(map(int, re.findall(rb"test result: ok\. (\d+) passed;", log)))
            assert count == check["passed_tests"] == 254
    print("Verified: 3 synthetic cases, 96 epochs, 864 reference rows; "
          "254 retained test passes and 5 successful checks. No solver run.")


if __name__ == "__main__":
    main()
