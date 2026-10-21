"""Verify carrier-lock admission evidence without rerunning positioning solvers."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import runpy
import statistics
import subprocess


def sha(data):return hashlib.sha256(data).hexdigest()
def read(path):return gzip.decompress(path.read_bytes())
def rows(data):return list(csv.DictReader(io.StringIO(data.decode())))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--dataset-root",type=Path);args=parser.parse_args()
    archive=Path(__file__).resolve().parent;repo=archive.parents[2]
    runner=runpy.run_path(str(repo/"scripts/run_ppc_carrier_lock_velocity.py"));rover=runpy.run_path(str(repo/"scripts/run_ppc_rover_velocity.py"));multi=runpy.run_path(str(repo/"scripts/run_ppc_multiclock_velocity.py"));targets=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))["reference_targets"]
    manifest=json.loads((archive/"predeclared-inputs.json").read_text());report=json.loads((archive/"report.json").read_text());totals={p:{c:0 for c in manifest["conditions"]} for p in manifest["profiles"]}
    assert not report["practical_PPC_parity_accepted"] and not manifest["velocity_fused"] and not manifest["numeric_thresholds_changed"] and not manifest["default_changed"]
    assert manifest["admission_chosen_after_post_run_raw_diagnostic"] and not manifest["held_out_validation"]
    previous=repo/"docs/comparison/2026-10-21-multiclock-velocity";candidates=repo/"docs/comparison/2026-10-21-mixed-candidate-validation"
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
                stem=case+"."+profile+"."+condition;execution=json.loads((archive/(stem+".execution.json")).read_text());data=read(archive/(stem+".stdout.csv.gz"));stderr=read(archive/(stem+".stderr.log.gz"))
                assert execution["returncode"]==0 and sha(data)==execution["csv_sha256"] and sha(stderr)==execution["stderr_sha256"]
                decoded=rows(data);actual=rover["summary"](data,control)
                if condition in {"GPS","mixed"}:
                    assert decoded==rows(read(previous/(case+"."+profile+"."+condition+".stdout.csv.gz")));actual["all_columns_identical_to_prior_control"]=True
                if condition.startswith("mixed"):actual["mixed"]=multi["mixed_summary"](decoded)
                if condition.endswith("_locked"):assert b"velocity_require_carrier_lock=true" in stderr
                expected=report["cases"][case][profile][condition];assert actual=={k:v for k,v in expected.items() if k!="reference_velocity_scores"}
                if references is not None:assert rover["scores"](decoded,references)==expected["reference_velocity_scores"]
                all_rows[condition]=decoded;totals[profile][condition]+=actual["velocity_counts"].get("VELOCITY",0)
            if references is not None:
                for a,b in [("GPS","GPS_locked"),("mixed","mixed_locked"),("GPS_locked","mixed_locked")]:
                    assert runner["comparison"](all_rows[a],all_rows[b],references,rover["scores"],multi,a,b)==report["cases"][case][profile]["common_epochs"][a+"_vs_"+b]
    probe=archive/"read-only-probe";pr=json.loads((probe/"report.json").read_text());old=json.loads((previous/"observability-diagnostics.json").read_text())
    for path,h in pr["source_sha256"].items():
        data=subprocess.check_output(["git","show",pr["base_commit"]+":"+path],cwd=repo) if path.startswith("src/") else (repo/path).read_bytes();assert sha(data)==h,path
    cross=json.loads((probe/"prior-frozen-input-crosscheck.json").read_text());assert cross["prior_manifest_sha256"]==sha((previous/"predeclared-inputs.json").read_bytes())
    for case,c in pr["cases"].items():
        data=read(probe/(case+".csv.gz"));assert c["returncode"]==0 and sha(data)==c["stdout_sha256"] and sha(read(probe/(case+".stderr.log.gz")))==c["stderr_sha256"]
        assert c["input_sha256"]==cross["verified_input_sha256"][case]
        raw=rows(data);used=dict(v.split(":") for v in old["cases"][case]["galileo"]["worst_reference_is_antenna"][0]["row"]["velocity_tracking_codes"].split(";"))
        for system,g in c["groups"].items():
            selected=[r for r in raw if used.get(r["satellite"])==r["tracking"] and r["satellite"][0]==system]
            values=[float(r["raw_phase_minus_Doppler_m"]) for r in selected if r["raw_phase_minus_Doppler_m"]];center=statistics.median(values) if values else None
            assert g==dict(rows=[dict(r,centered_residual_m=float(r["raw_phase_minus_Doppler_m"])-center if r["raw_phase_minus_Doppler_m"] else None) for r in selected],median_raw_m=center)
    q=json.loads((archive/"qualification.json").read_text());assert q["predeclared_inputs_sha256"]==sha((archive/"predeclared-inputs.json").read_bytes());assert {c["name"] for c in q["checks"]}=={"rust","examples","fmt","clippy","release"}
    for c in q["checks"]:
        data=read(archive/c["log"]);assert c["returncode"]==0 and sha(data)==c["uncompressed_sha256"]
        if c["name"]=="rust":assert sum(map(int,re.findall(rb"test result: ok\. (\d+) passed;",data)))==c["passed_tests"]==344
    for path,h in json.loads((archive/"evidence-sha256.json").read_text()).items():assert sha((archive/path).read_bytes())==h,path
    print("Verified 48 velocity runs, six read-only raw probes and 344 tests; velocities:",totals)


if __name__=="__main__":main()
