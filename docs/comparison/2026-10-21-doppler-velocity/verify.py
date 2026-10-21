"""Verify single-system Doppler evidence without rerunning solvers."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import runpy


def sha(data):return hashlib.sha256(data).hexdigest()
def read(path):return gzip.decompress(path.read_bytes())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root",type=Path)
    args=parser.parse_args();archive=Path(__file__).resolve().parent;repo=archive.parents[2]
    runner=runpy.run_path(str(repo/"scripts/run_ppc_base_velocity.py"))
    manifest=json.loads((archive/"predeclared-inputs.json").read_text());report=json.loads((archive/"report.json").read_text())
    previous=repo/"docs/comparison/2026-10-21-gps-galileo-float"
    assert not report["practical_PPC_parity_accepted"] and not manifest["velocity_fused_into_RTK"] and not manifest["FIX_policy_evaluated"]
    assert sha((previous/"predeclared-inputs.json").read_bytes())==manifest["prior_manifest_sha256"]
    for path,h in manifest["source_sha256"].items():assert sha((repo/path).read_bytes())==h,path
    total=0
    for case in manifest["cases"]:
        execution=json.loads((archive/(case+".execution.json")).read_text())
        data=read(archive/(case+".stdout.csv.gz"));stderr=read(archive/(case+".stderr.log.gz"))
        assert execution["returncode"]==0 and sha(data)==execution["csv_sha256"] and sha(stderr)==execution["stderr_sha256"]
        metrics=runner["summary"](data);metrics["all_base_epoch_keys_identical_to_prior_ledger"]=True
        assert metrics==report["cases"][case];total+=metrics["epochs"]
        assert b"velocity_not_fused=true" in stderr and b"no_FIX_policy=true" in stderr
        expected=[(r["week"],float(r["tow"]).hex()) for r in csv.DictReader(io.StringIO(read(previous/(case+".gps.stdout.csv.gz")).decode())) if r["kind"]=="pair" or r["side"]=="base"]
        assert [(r["week"],float(r["tow"]).hex()) for r in csv.DictReader(io.StringIO(data.decode()))]==expected
        if args.dataset_root:
            city,iteration=case.split("-")
            for name in ["base.obs","base.nav"]:
                path=args.dataset_root/city/iteration/name;expected=next(h for k,h in manifest["input_sha256"].items() if k.endswith(f"/{city}/{iteration}/{name}"))
                assert sha(path.read_bytes())==expected
    qualification=json.loads((archive/"qualification.json").read_text())
    assert qualification["predeclared_inputs_sha256"]==sha((archive/"predeclared-inputs.json").read_bytes())
    assert {c["name"] for c in qualification["checks"]}=={"rust","examples","fmt","clippy","release"}
    for c in qualification["checks"]:
        data=read(archive/c["log"]);assert c["returncode"]==0 and sha(data)==c["uncompressed_sha256"]
        if c["name"]=="rust":assert sum(map(int,re.findall(rb"test result: ok\. (\d+) passed;",data)))==c["passed_tests"]==337
    component=json.loads((archive/"component-generation.json").read_text())
    assert component["compile_returncode"]==component["run_returncode"]==0 and component["regenerated_bytes_identical"]
    assert sha((repo/"scripts/upstream_velocity.cpp").read_bytes())==component["generator_sha256"]
    assert sha((repo/"tests/fixtures/upstream_velocity.csv").read_bytes())==component["fixture_sha256"]
    assert sha(read(archive/"component-build.log.gz"))==component["compile_log_sha256"]
    for path,h in json.loads((archive/"evidence-sha256.json").read_text()).items():assert sha((archive/path).read_bytes())==h,path
    print("Verified six stationary-base PPC runs:",total,"epochs; 12 native cases; 337 tests.")


if __name__=="__main__":main()
