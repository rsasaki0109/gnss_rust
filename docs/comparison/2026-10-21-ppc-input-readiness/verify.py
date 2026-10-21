"""Verify retained input-readiness evidence without a solver or network call."""
import gzip
import hashlib
import json
from pathlib import Path
import re
import runpy


def main():
    archive = Path(__file__).resolve().parent
    repo = archive.parents[2]
    helpers = runpy.run_path(str(repo / "scripts/run_ppc_input_readiness.py"))
    assert helpers["verify"](repo, archive) == json.loads((archive / "report.json").read_text())
    qualification = json.loads((archive / "qualification.json").read_text())
    assert hashlib.sha256((archive / "predeclared-inputs.json").read_bytes()).hexdigest() == (
        qualification["predeclared_inputs_sha256"]
    )
    assert {r["name"] for r in qualification["checks"]} == {
        "rust", "examples", "fmt", "clippy", "release"
    }
    for check in qualification["checks"]:
        assert check["returncode"] == 0
        log = gzip.decompress((archive / check["log"]).read_bytes())
        assert hashlib.sha256(log).hexdigest() == check["uncompressed_sha256"]
        if check["name"] == "rust":
            assert sum(map(int, re.findall(rb"test result: ok\. (\d+) passed;", log))) == check["passed_tests"] == 266
    acquisition = json.loads((archive / "product-acquisition.json").read_text())
    assert acquisition["acquired_count"] == len(acquisition["products"]) == 27
    assert all(p["acquired"] for p in acquisition["products"])
    assert acquisition["source_sha256"] == helpers["sha"](repo / "scripts/acquire_ppc_precise.py")
    print("Verified 6-run / 9-day readiness evidence, 27 acquired products, "
          "266 retained test passes and 5 checks. No positioning acceptance claim.")


if __name__ == "__main__":
    main()
