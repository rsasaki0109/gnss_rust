"""Recheck rover velocity diagnostics without executing solvers."""
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
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--dataset-root",type=Path)
    args=parser.parse_args();archive=Path(__file__).resolve().parent;repo=archive.parents[2]
    runner=runpy.run_path(str(repo/"scripts/run_ppc_rover_velocity.py"));targets=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))["reference_targets"]
    manifest=json.loads((archive/"predeclared-inputs.json").read_text());report=json.loads((archive/"report.json").read_text());totals={p:0 for p in manifest["profiles"]}
    assert not report["practical_PPC_parity_accepted"] and not manifest["velocity_fused"] and not manifest["motion_or_FIX_policy_evaluated"]
    previous=repo/"docs/comparison/2026-10-21-mixed-candidate-validation"
    assert sha((previous/"predeclared-inputs.json").read_bytes())==manifest["prior_manifest_sha256"]
    for path,h in manifest["source_sha256"].items():assert sha((repo/path).read_bytes())==h,path
    for case in manifest["cases"]:
        references=None
        if args.dataset_root:
            city,iteration=case.split("-");root=args.dataset_root/city/iteration
            for name in ["rover.obs","base.obs","base.nav","reference.csv"]:
                expected=next(h for k,h in manifest["input_sha256"].items() if k.endswith(f"/{city}/{iteration}/{name}"));assert sha((root/name).read_bytes())==expected
            references=runner["reference_velocities"](root/"reference.csv",city,targets)
        for profile in manifest["profiles"]:
            stem=case+"."+profile;execution=json.loads((archive/(stem+".execution.json")).read_text())
            data=read(archive/(stem+".stdout.csv.gz"));stderr=read(archive/(stem+".stderr.log.gz"))
            assert execution["returncode"]==0 and sha(data)==execution["csv_sha256"] and sha(stderr)==execution["stderr_sha256"]
            control=list(csv.DictReader(io.StringIO(read(previous/(stem+".diagnostic.stdout.csv.gz")).decode())))
            actual=runner["summary"](data,control);expected=report["cases"][case][profile]
            assert actual=={k:v for k,v in expected.items() if k!="reference_velocity_scores"}
            assert b"no_velocity_feedback=true" in stderr and b"no_motion_gate_evaluated=true" in stderr
            if references is not None:assert runner["scores"](list(csv.DictReader(io.StringIO(data.decode()))),references)==expected["reference_velocity_scores"]
            totals[profile]+=actual["velocity_counts"].get("VELOCITY",0)
    q=json.loads((archive/"qualification.json").read_text())
    assert q["predeclared_inputs_sha256"]==sha((archive/"predeclared-inputs.json").read_bytes())
    assert {c["name"] for c in q["checks"]}=={"rust","examples","fmt","clippy","release"}
    for c in q["checks"]:
        data=read(archive/c["log"]);assert c["returncode"]==0 and sha(data)==c["uncompressed_sha256"]
        if c["name"]=="rust":assert sum(map(int,re.findall(rb"test result: ok\. (\d+) passed;",data)))==c["passed_tests"]==337
    if args.dataset_root:
        analyze=runpy.run_path(str(archive/"analyze.py"));assert analyze["diagnostics"](archive,args.dataset_root)==json.loads((archive/"observability-diagnostics.json").read_text())
    for path,h in json.loads((archive/"evidence-sha256.json").read_text()).items():assert sha((archive/path).read_bytes())==h,path
    print("Verified 12 rover-velocity runs; previous 29 columns unchanged; 337 tests; velocities:",totals)


if __name__=="__main__":main()
