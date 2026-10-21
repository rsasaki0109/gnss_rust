"""Freeze GPS vs GPS/Galileo independent-clock velocity on matching FLOAT snapshots."""
import argparse
import collections
import csv
import gzip
import io
import json
from pathlib import Path
import runpy
import subprocess


def mixed_summary(rows):
    counts=collections.Counter();tracking=collections.Counter()
    for r in rows:
        if r["status"]!="FLOAT":continue
        choices=[v.split(":") for v in r["velocity_tracking_codes"].split(";") if v]
        assert len(choices)==int(r["velocity_prepared_count"])
        for sat,code in choices:
            assert sat[0] in {"G","E"} and code in {"1C","1X"}
            if sat[0]=="G":assert code=="1C"
            tracking[sat[0]+":"+code]+=1
        if r["velocity_status"]!="VELOCITY":continue
        systems=r["velocity_clock_systems"].split(";");drifts=r["velocity_clock_drifts_sps"].split(";")
        used=r["velocity_used_satellites"].split(";")
        assert set(systems)=={s[0] for s in used} and len(systems)==len(set(systems))==len(drifts)
        assert set(used)=={s for s,_ in choices}
        assert len(used)>=3+len(systems)
        counts["two_clock_epochs"]+=len(systems)==2;counts["single_clock_epochs"]+=len(systems)==1
    return dict(counts=dict(counts),tracking_rows=dict(tracking))


def compare(gps,mixed,references,score):
    key=lambda r:(r["week"],r["tow"])
    a={key(r):r for r in gps if r["velocity_status"]=="VELOCITY"};b={key(r):r for r in mixed if r["velocity_status"]=="VELOCITY"}
    common=sorted(a.keys()&b.keys())
    return dict(common_velocity_epochs=len(common),GPS_only=len(a.keys()-b.keys()),mixed_only=len(b.keys()-a.keys()),
        common_reference_scores=dict(GPS=score([a[k] for k in common],references),mixed=score([b[k] for k in common],references)))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--dataset-root",type=Path,required=True);parser.add_argument("--binary",type=Path,required=True);parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args();repo=Path(__file__).resolve().parents[1]
    previous=repo/"docs/comparison/2026-10-21-rover-velocity";candidates=repo/"docs/comparison/2026-10-21-mixed-candidate-validation"
    prior=json.loads((previous/"predeclared-inputs.json").read_text())
    helpers=runpy.run_path(str(repo/"scripts/run_native_ppc_comparison.py"));sha,save,retain=(helpers[k] for k in ["sha","save","retain"])
    rover=runpy.run_path(str(repo/"scripts/run_ppc_rover_velocity.py"));targets=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))["reference_targets"]
    dataset=args.dataset_root.resolve();binary=args.binary.resolve();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    sources=list((repo/"src").rglob("*.rs"))+[repo/p for p in ["scripts/run_ppc_multiclock_velocity.py","scripts/run_ppc_rover_velocity.py","scripts/run_ppc_rtk_float.py","scripts/run_native_ppc_comparison.py","examples/ppc_rtk_candidates.rs","tests/velocity.rs","tests/multiclock_velocity.rs","scripts/upstream_multiclock_velocity.cpp","scripts/upstream_velocity.cpp","tests/fixtures/upstream_multiclock_velocity.csv","Cargo.toml","Cargo.lock","rust-toolchain.toml"]]
    inputs={}
    for case in prior["cases"]:
        city,iteration=case.split("-")
        for name in ["rover.obs","base.obs","base.nav","reference.csv"]:
            p=dataset/city/iteration/name;h=next(h for k,h in prior["input_sha256"].items() if k.endswith(f"/{city}/{iteration}/{name}"));assert sha(p)==h;inputs[str(p)]=h
    manifest=dict(base_commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=repo,text=True).strip(),source_sha256={str(p.relative_to(repo)):sha(p) for p in sources},input_sha256=inputs,binary_sha256=sha(binary),
        cases=prior["cases"],base_ecef_m=prior["base_ecef_m"],profiles=prior["profiles"],conditions=["GPS","mixed"],
        changed_setting="velocity only: add Galileo E1, separate receiver clock-range-rate per system; GPS code/phase FLOAT unchanged",
        Galileo_selection="valid Doppler first, valid code next, then 1C before 1X; alternatives retained",thresholds_changed=False,default_changed=False,
        reference_used_for_configuration=False,velocity_fused=False,motion_or_FIX_policy_evaluated=False,native_whole_solver_rerun=False,practical_PPC_parity_accepted=False,
        prior_manifest_sha256=sha(previous/"predeclared-inputs.json"))
    save(out/"predeclared-inputs.json",manifest);report=dict(practical_PPC_parity_accepted=False,cases={})
    for case in manifest["cases"]:
        city,iteration=case.split("-");root=dataset/city/iteration;report["cases"][case]={}
        references=rover["reference_velocities"](root/"reference.csv",city,targets)
        for profile,mode in manifest["profiles"].items():
            with gzip.open(candidates/(case+"."+profile+".diagnostic.stdout.csv.gz"),"rt") as stream:control=list(csv.DictReader(stream))
            report["cases"][case][profile]={};all_rows={}
            for condition in manifest["conditions"]:
                assert all(sha(repo/p)==h for p,h in manifest["source_sha256"].items()) and sha(binary)==manifest["binary_sha256"]
                selected_mode=mode if condition=="GPS" else mode+"-mixed"
                command=[str(binary),str(root),selected_mode,",".join(map(str,manifest["base_ecef_m"][city]))]
                r=subprocess.run(command,capture_output=True,cwd=repo);retain(out,case+"."+profile+"."+condition,command,r.stdout,r.stderr,r.returncode);assert r.returncode==0,r.stderr.decode()
                rows=list(csv.DictReader(io.StringIO(r.stdout.decode())));metrics=rover["summary"](r.stdout,control)
                if condition=="GPS":
                    with gzip.open(previous/(case+"."+profile+".stdout.csv.gz"),"rt") as stream:assert rows==list(csv.DictReader(stream))
                    metrics["all_51_columns_identical_to_prior_GPS_control"]=True
                else:metrics["mixed"]=mixed_summary(rows)
                metrics["reference_velocity_scores"]=rover["scores"](rows,references)
                report["cases"][case][profile][condition]=metrics;all_rows[condition]=rows;save(out/"report.json",report)
                print(case,profile,condition,metrics["velocity_counts"],flush=True)
            report["cases"][case][profile]["common_epochs"]=compare(all_rows["GPS"],all_rows["mixed"],references,rover["scores"])
            save(out/"report.json",report)
    assert all(sha(Path(p))==h for p,h in inputs.items())


if __name__=="__main__":main()
