"""Recheck frozen multi-clock velocity evidence without running GNSS solvers."""
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
def rows(data):return list(csv.DictReader(io.StringIO(data.decode())))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--dataset-root",type=Path)
    args=parser.parse_args();archive=Path(__file__).resolve().parent;repo=archive.parents[2]
    runner=runpy.run_path(str(repo/"scripts/run_ppc_multiclock_velocity.py"));rover=runpy.run_path(str(repo/"scripts/run_ppc_rover_velocity.py"));targets=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))["reference_targets"]
    manifest=json.loads((archive/"predeclared-inputs.json").read_text());report=json.loads((archive/"report.json").read_text());totals={p:{c:0 for c in manifest["conditions"]} for p in manifest["profiles"]}
    assert not report["practical_PPC_parity_accepted"] and not manifest["velocity_fused"] and not manifest["thresholds_changed"] and not manifest["default_changed"]
    previous=repo/"docs/comparison/2026-10-21-rover-velocity";candidates=repo/"docs/comparison/2026-10-21-mixed-candidate-validation"
    assert sha((previous/"predeclared-inputs.json").read_bytes())==manifest["prior_manifest_sha256"]
    for path,h in manifest["source_sha256"].items():assert sha((repo/path).read_bytes())==h,path
    for case in manifest["cases"]:
        references=None
        if args.dataset_root:
            city,iteration=case.split("-");root=args.dataset_root/city/iteration
            for name in ["rover.obs","base.obs","base.nav","reference.csv"]:
                expected=next(h for k,h in manifest["input_sha256"].items() if k.endswith(f"/{city}/{iteration}/{name}"));assert sha((root/name).read_bytes())==expected
            references=rover["reference_velocities"](root/"reference.csv",city,targets)
        for profile in manifest["profiles"]:
            control=rows(read(candidates/(case+"."+profile+".diagnostic.stdout.csv.gz")));all_rows={}
            for condition in manifest["conditions"]:
                stem=case+"."+profile+"."+condition;execution=json.loads((archive/(stem+".execution.json")).read_text())
                data=read(archive/(stem+".stdout.csv.gz"));stderr=read(archive/(stem+".stderr.log.gz"))
                assert execution["returncode"]==0 and sha(data)==execution["csv_sha256"] and sha(stderr)==execution["stderr_sha256"]
                decoded=rows(data);actual=rover["summary"](data,control)
                if condition=="GPS":
                    assert decoded==rows(read(previous/(case+"."+profile+".stdout.csv.gz")))
                    actual["all_51_columns_identical_to_prior_GPS_control"]=True
                else:
                    actual["mixed"]=runner["mixed_summary"](decoded)
                    assert b"GPS_Galileo_separate_clocks" in stderr
                expected=report["cases"][case][profile][condition]
                assert actual=={k:v for k,v in expected.items() if k!="reference_velocity_scores"}
                if references is not None:assert rover["scores"](decoded,references)==expected["reference_velocity_scores"]
                all_rows[condition]=decoded;totals[profile][condition]+=actual["velocity_counts"].get("VELOCITY",0)
            if references is not None:assert runner["compare"](all_rows["GPS"],all_rows["mixed"],references,rover["scores"])==report["cases"][case][profile]["common_epochs"]
    q=json.loads((archive/"qualification.json").read_text());assert q["predeclared_inputs_sha256"]==sha((archive/"predeclared-inputs.json").read_bytes())
    assert {c["name"] for c in q["checks"]}=={"rust","examples","fmt","clippy","release"}
    for c in q["checks"]:
        data=read(archive/c["log"]);assert c["returncode"]==0 and sha(data)==c["uncompressed_sha256"]
        if c["name"]=="rust":assert sum(map(int,re.findall(rb"test result: ok\. (\d+) passed;",data)))==c["passed_tests"]==343
    component=json.loads((archive/"component-generation.json").read_text());assert component["compile_returncode"]==component["run_returncode"]==0 and component["regenerated_bytes_identical"]
    assert sha((repo/"scripts/upstream_multiclock_velocity.cpp").read_bytes())==component["generator_sha256"]
    assert sha((repo/"scripts/upstream_velocity.cpp").read_bytes())==component["helper_sha256"]
    assert sha((repo/"tests/fixtures/upstream_multiclock_velocity.csv").read_bytes())==component["fixture_sha256"]
    assert sha(read(archive/"component-build.log.gz"))==component["compile_log_sha256"]
    if args.dataset_root:
        analyze=runpy.run_path(str(archive/"analyze.py"));assert analyze["diagnostics"](archive,args.dataset_root)==json.loads((archive/"observability-diagnostics.json").read_text())
    for path,h in json.loads((archive/"evidence-sha256.json").read_text()).items():assert sha((archive/path).read_bytes())==h,path
    print("Verified 24 PPC runs, unchanged GPS controls, 16 independent Eigen cases and 343 tests; velocities:",totals)


if __name__=="__main__":main()
