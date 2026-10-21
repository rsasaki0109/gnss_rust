"""Verify locked-subset policy evidence without running positioning solvers."""
import argparse
import csv
from decimal import Decimal
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
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--dataset-root",type=Path);args=parser.parse_args()
    archive=Path(__file__).resolve().parent;repo=archive.parents[2];previous=repo/"docs/comparison/2026-10-21-kinematic-fix"
    m=json.loads((archive/"predeclared-inputs.json").read_text());report=json.loads((archive/"report.json").read_text());runner=runpy.run_path(str(repo/"scripts/run_ppc_locked_subset.py"))
    score=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))["score"]
    assert not m["practical_PPC_parity_accepted"] and not report["practical_PPC_parity_accepted"] and not m["native_rerun"] and not m["held_out_validation"]
    assert not m["numeric_thresholds_changed"] and not m["default_static_or_full_set_policy_changed"] and not m["ratio_driven_subset_search"] and not m["integer_hold"] and not m["FLOAT_feedback"]
    assert m["min_lock_epochs"]==5 and m["min_ambiguity_pairs"]==4 and m["kinematic_default_motion_std_mps"]==1 and m["predeclared_motion_std_sensitivity"]==[1,2]
    assert len(m["cases"])==6 and len(m["profiles"])==16
    assert sha((previous/"predeclared-inputs.json").read_bytes())==m["prior_manifest_sha256"]
    for path,h in m["source_sha256"].items():assert sha((repo/path).read_bytes())==h,path
    for case in m["cases"]:
        city,iteration=case.split("-");references=None;controls={};lock_ledgers={}
        if args.dataset_root:
            root=args.dataset_root/city/iteration
            for name in ["rover.obs","base.obs","base.nav","reference.csv"]:
                h=next(v for p,v in m["input_sha256"].items() if p.endswith(f"/{city}/{iteration}/{name}"));assert sha((root/name).read_bytes())==h
            with (root/"reference.csv").open() as f:refs=[{k.strip():v.strip() for k,v in r.items()} for r in csv.DictReader(f)]
            references={(int(r["GPS Week"]),Decimal(r["GPS TOW (s)"])):r for r in refs}
        for profile,c in m["profiles"].items():
            old=c["control_profile"]
            if old not in controls:controls[old]=rows(read(previous/(case+"."+old+".stdout.csv.gz")))
            stem=case+"."+profile;e=json.loads((archive/(stem+".execution.json")).read_text());data=read(archive/(stem+".stdout.csv.gz"));stderr=read(archive/(stem+".stderr.log.gz"))
            assert e["returncode"]==0 and sha(data)==e["csv_sha256"] and sha(stderr)==e["stderr_sha256"]
            assert ("locked_subset="+str(c["locked_subset"]).lower()).encode() in stderr and b"fixed_or_velocity_feedback=false" in stderr and b"integer_holding=false" in stderr
            decoded,valid,fixed,actual=runner["summarize"](data,controls[old],c["locked_subset"],c["motion_std_mps"],c["mode"].startswith("all"))
            ledger=[r["satellite_lock_epochs"] for r in decoded]
            assert ledger==lock_ledgers.setdefault(c["base_profile"],ledger)
            expected=report["cases"][case][profile];assert actual=={k:v for k,v in expected.items() if k not in {"all_RTK_reference_scores","FIX_reference_scores"}}
            if references is not None:
                for label,part in [("all_RTK",valid),("FIX",fixed)]:assert score([{**r,"status":"FLOAT"} for r in part],references,city)==expected[label+"_reference_scores"]
            history={}
            for r in valid:
                t=Decimal(r["week"])*604800+Decimal(r["tow"])
                if r["velocity_status"]=="VELOCITY":history[t]=r
                for check in r["motion_checks"].split(";") if r["motion_checks"] else []:
                    f=check.split("@");dt=float(f[1]);past=history[t-Decimal(str(dt))]
                    for j,k in enumerate(["vx_mps","vy_mps","vz_mps"]):assert (.5*float(past[k])+.5*float(r[k]))*dt==float(f[2+j])
    analysis=json.loads((archive/"analysis.json").read_text());actual=runpy.run_path(str(archive/"analyze.py"))["analyze"](archive,args.dataset_root)
    if args.dataset_root is None:
        for case in analysis["cases"].values():
            for p in case["profiles"].values():p.pop("common_RTK_reference_scores",None)
    assert analysis==actual
    native=json.loads((archive/"native-generator.json").read_text());fixture=repo/"tests/fixtures/upstream_rtk_locked_subset.csv"
    assert native["returncode"]==0 and native["cases"]==24 and len(fixture.read_text().splitlines())==25
    assert native["source_sha256"]==sha((repo/"scripts/upstream_rtk_locked_subset.cpp").read_bytes()) and native["build_log_sha256"]==sha((archive/"native-generator-build.log.gz").read_bytes())
    assert all(r["returncode"]==0 and r["fixture_sha256"]==sha(fixture.read_bytes()) for r in native["repeat_runs"])
    baseline=json.loads((repo/"docs/comparison/2026-10-21-gps-galileo-float/native-comparison/predeclared-inputs.json").read_text())
    assert all(baseline["native_archive_sha256"][path]==h for path,h in native["native_library_sha256"].items())
    q=json.loads((archive/"qualification.json").read_text());assert q["predeclared_inputs_sha256"]==sha((archive/"predeclared-inputs.json").read_bytes())
    assert {c["name"] for c in q["checks"]}=={"rust","examples","fmt","clippy","release"}
    for c in q["checks"]:
        data=read(archive/c["log"]);assert c["returncode"]==0 and sha(data)==c["uncompressed_sha256"]
        if c["name"]=="rust":assert sum(map(int,re.findall(rb"test result: ok\. (\d+) passed;",data)))==c["passed_tests"]==357
    sealed=json.loads((archive/"evidence-sha256.json").read_text());assert set(sealed)=={str(p.relative_to(archive)) for p in archive.rglob("*") if p.is_file() and p.name!="evidence-sha256.json"}
    for path,h in sealed.items():assert sha((archive/path).read_bytes())==h,path
    print("Verified 96 PPC runs, 48 unchanged full-set controls, 24 native component cases and 357 Rust tests")


if __name__=="__main__":main()
