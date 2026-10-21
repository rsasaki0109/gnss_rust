"""Verify retained kinematic policy evidence without rerunning positioning solvers."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
import runpy


def sha(data):return hashlib.sha256(data).hexdigest()
def read(path):return gzip.decompress(path.read_bytes())
def rows(data):return list(csv.DictReader(io.StringIO(data.decode())))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--dataset-root",type=Path);args=parser.parse_args()
    archive=Path(__file__).resolve().parent;repo=archive.parents[2]
    manifest=json.loads((archive/"predeclared-inputs.json").read_text());report=json.loads((archive/"report.json").read_text())
    runner=runpy.run_path(str(repo/"scripts/run_ppc_kinematic_fix.py"));score=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))["score"]
    assert not manifest["practical_PPC_parity_accepted"] and not report["practical_PPC_parity_accepted"] and not manifest["native_rerun"]
    assert manifest["full_set_only"] and not manifest["partial_AR"] and not manifest["integer_hold"] and not manifest["default_static_policy_changed"]
    assert manifest["kinematic_default_motion_std_mps"]==1 and manifest["predeclared_motion_std_sensitivity_after_shape_only_smoke"]==[1,2,4]
    assert not manifest["numeric_local_integer_thresholds_changed"] and not manifest["held_out_validation"] and not manifest["reference_used_for_configuration"]
    for path,h in manifest["source_sha256"].items():assert sha((repo/path).read_bytes())==h,path
    previous=repo/"docs/comparison/2026-10-21-carrier-lock-velocity";oldfloat=repo/"docs/comparison/2026-10-21-gps-galileo-float"
    assert sha((previous/"predeclared-inputs.json").read_bytes())==manifest["prior_manifest_sha256"]
    assert len(manifest["cases"])==6 and len(manifest["profiles"])==12
    for case in manifest["cases"]:
        city,iteration=case.split("-");references=None
        if args.dataset_root:
            root=args.dataset_root/city/iteration
            for name in ["rover.obs","base.obs","base.nav","reference.csv"]:
                expected=next(h for p,h in manifest["input_sha256"].items() if p.endswith(f"/{city}/{iteration}/{name}"));assert sha((root/name).read_bytes())==expected
            with (root/"reference.csv").open() as f:
                refs=[{k.strip():v.strip() for k,v in r.items()} for r in csv.DictReader(f)]
            references={(int(r["GPS Week"]),Decimal(r["GPS TOW (s)"])):r for r in refs}
        controls={p:(rows(read(oldfloat/(case+"."+p+".stdout.csv.gz"))),rows(read(previous/(case+"."+p+".mixed_locked.stdout.csv.gz")))) for p in ["galileo","galileo_slips"]}
        for profile,config in manifest["profiles"].items():
            stem=case+"."+profile;e=json.loads((archive/(stem+".execution.json")).read_text());data=read(archive/(stem+".stdout.csv.gz"));stderr=read(archive/(stem+".stderr.log.gz"))
            assert e["returncode"]==0 and sha(data)==e["csv_sha256"] and sha(stderr)==e["stderr_sha256"]
            assert b"fixed_or_velocity_feedback=false" in stderr and b"integer_holding=false" in stderr and b"performance_accepted=false" in stderr
            decoded,valid,fixed,actual=runner["summary"](data,*controls[config["base_profile"]],config["motion_std_mps"],config["mode"].startswith("all"))
            expected=report["cases"][case][profile];assert actual=={k:v for k,v in expected.items() if k not in {"all_RTK_reference_scores","FIX_reference_scores"}}
            if references is not None:
                assert runner["scores"](valid,references,city,score)==expected["all_RTK_reference_scores"]
                assert runner["scores"](fixed,references,city,score)==expected["FIX_reference_scores"]
            # Independently check every published prediction. Candidate positions
            # are intentionally not reconstructed from FLOAT positions.
            history={}
            for r in decoded:
                if r["status"] not in {"FLOAT","FIX"}:continue
                t=Decimal(r["week"])*604800+Decimal(r["tow"])
                if r["velocity_status"]=="VELOCITY":history[t]=r
                for check in r["motion_checks"].split(";") if r["motion_checks"] else []:
                    f=check.split("@");dt=float(f[1]);past=history[t-Decimal(str(dt))]
                    for j,k in enumerate(["vx_mps","vy_mps","vz_mps"]):
                        prediction=(.5*float(past[k])+.5*float(r[k]))*dt
                        assert prediction==float(f[2+j])
                        if f[0]=="LastFixed" and r["status"]=="FIX":
                            assert past["status"]=="FIX"
                            residual=float(r[["x_m","y_m","z_m"][j]])-float(past[["x_m","y_m","z_m"][j]])-prediction
                            # Published ECEF has nine decimal places; each end
                            # can round by 0.5 nm, plus ECEF arithmetic roundoff.
                            assert math.isclose(residual,float(f[5+j]),rel_tol=0,abs_tol=3e-9)
    analysis=json.loads((archive/"analysis.json").read_text());recomputed=runpy.run_path(str(archive/"analyze.py"))["analyze"](archive,args.dataset_root)
    if args.dataset_root is None:
        for c in analysis["cases"].values():
            for p in c["profiles"].values():p.pop("common_RTK_reference_scores",None)
    assert analysis==recomputed
    q=json.loads((archive/"qualification.json").read_text());assert q["predeclared_inputs_sha256"]==sha((archive/"predeclared-inputs.json").read_bytes())
    assert {c["name"] for c in q["checks"]}=={"rust","examples","fmt","clippy","release"}
    for c in q["checks"]:
        data=read(archive/c["log"]);assert c["returncode"]==0 and sha(data)==c["uncompressed_sha256"]
        if c["name"]=="rust":assert sum(map(int,re.findall(rb"test result: ok\. (\d+) passed;",data)))==c["passed_tests"]==352
    failed=archive/"initial-incomplete";fm=json.loads((failed/"predeclared-inputs.json").read_text());fr=json.loads((failed/"report.json").read_text())
    assert sha((failed/"runner-before-comparator-fix.py").read_bytes())==fm["source_sha256"]["scripts/run_ppc_kinematic_fix.py"]
    assert sum(map(len,fr["cases"].values()))==12 and len(list(failed.glob("*.execution.json")))==13
    for path in failed.glob("*.execution.json"):
        e=json.loads(path.read_text());stem=path.name.removesuffix(".execution.json")
        assert e["returncode"]==0 and sha(read(failed/(stem+".stdout.csv.gz")))==e["csv_sha256"] and sha(read(failed/(stem+".stderr.log.gz")))==e["stderr_sha256"]
    sealed=json.loads((archive/"evidence-sha256.json").read_text())
    assert set(sealed)=={str(p.relative_to(archive)) for p in archive.rglob("*") if p.is_file() and p.name!="evidence-sha256.json"}
    for path,h in sealed.items():assert sha((archive/path).read_bytes())==h,path
    print("Verified 72 kinematic policy runs, unchanged FLOAT/Doppler controls, frozen native comparison and 352 tests")


if __name__=="__main__":main()
