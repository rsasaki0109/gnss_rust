"""Check retained CLI evidence and qualification logs without solver execution."""
import gzip
import hashlib
import json
from pathlib import Path
import re
import runpy


def main():
    archive = Path(__file__).resolve().parent
    repo = archive.parents[2]
    helpers = runpy.run_path(str(repo / "scripts/compare_rtk_reference_gates.py"))
    report = helpers["verify_archive"](repo, archive)
    assert report == json.loads((archive / "report.json").read_text())
    qualification = json.loads((archive / "qualification.json").read_text())
    assert hashlib.sha256((archive / "predeclared-inputs.json").read_bytes()).hexdigest() == (
        qualification["predeclared_inputs_sha256"]
    )
    assert len(qualification["checks"]) == 5
    assert {r["name"] for r in qualification["checks"]} == {
        "rust", "examples", "fmt", "clippy", "release"
    }
    for check in qualification["checks"]:
        assert check["returncode"] == 0
        data = gzip.decompress((archive / check["log"]).read_bytes())
        assert hashlib.sha256(data).hexdigest() == check["uncompressed_sha256"]
        if check["name"] == "rust":
            count = sum(map(int, re.findall(rb"test result: ok\. (\d+) passed;", data)))
            assert count == check["passed_tests"] == 260
    print("Verified 8 CLI runs / 256 epoch rows, 260 retained test passes "
          "and 5 successful checks; no solver executed.")


if __name__ == "__main__":
    main()
