"""Freeze full-set kinematic FIX and a declared motion-uncertainty sensitivity grid."""
import argparse
import collections
import csv
from decimal import Decimal
import gzip
import io
import json
import math
from pathlib import Path
import runpy
import subprocess


def summary(data,float_control,velocity_control,motion_std,all_references):
    rows=list(csv.DictReader(io.StringIO(data.decode())));assert len(rows)==len(float_control)==len(velocity_control)
    decisions=collections.Counter();fixed=[];valid=[];checks=0
    for r,b,v in zip(rows,float_control,velocity_control):
        for k in ["kind","week","tow","side"]:assert r[k]==b[k]==v[k]
        if b["status"]!="FLOAT":
            assert r["status"]==b["status"] and all(not r[k] for k in ["x_m","y_m","z_m","float_x_m","float_y_m","float_z_m"]);continue
        assert r["status"] in {"FLOAT","FIX"}
        for k in ["x","y","z"]:assert r["float_"+k+"_m"]==b[k+"_m"]
        for k in ["used_satellites","reference_satellites","ambiguity_resets"]:assert r[k]==b[k]
        assert r["ambiguity_pairs"]==b["dd_phase_rows"]
        for k in ["code_rms_m","phase_rms_m"]:assert r["float_"+k]==b[k]
        for axis in ["xx","yy","zz"]:assert r["float_q"+axis+"_m2"]==b["position_q"+axis+"_m2"]
        assert r["velocity_status"]==v["velocity_status"]
        if r["velocity_status"]=="VELOCITY":
            for k in ["vx_mps","vy_mps","vz_mps"]:assert r[k]==v[k]
            for a in ["xx","xy","xz","yx","yy","yz","zx","zy","zz"]:assert r["velocity_q"+a]==v["vq"+a]
            std=math.sqrt(max(float(v[k]) for k in ["vqxx","vqyy","vqzz"]));assert float(r["velocity_max_std_mps"])==std
        decisions[r["decision"]]+=1;valid.append(r)
        if r["motion_checks"]:
            for item in r["motion_checks"].split(";"):
                f=item.split("@");assert len(f)==10;dt=float(f[1]);assert 0<dt<=2
                residual=[float(x) for x in f[5:8]];assert math.isclose(math.hypot(*residual),float(f[8]),abs_tol=1e-14)
                assert float(f[9])==.1;checks+=1
        if r["status"]=="FIX":
            assert r["decision"]=="accepted" and int(r["confirmation_count"])>=2 and float(r["ratio"])>=3
            assert int(r["ambiguity_pairs"])>=4 and float(r["velocity_max_std_mps"])<=motion_std
            assert r["motion_checks"] and float(r["motion_max_residual_m"])<=.1
            if all_references:
                counts=collections.Counter(s[0] for s in r["used_satellites"].split(";"));assert int(r["system_reference_sets"])==math.prod(counts.values())
            fixed.append(r)
        else:
            assert r["decision"]!="accepted"
            for k in ["x_m","y_m","z_m"]:assert r[k]==r["float_"+k]
    pairs=[r for r in rows if r["kind"]=="pair"]
    return rows,valid,fixed,dict(output_rows=len(rows),paired_epochs=len(pairs),status_counts=dict(collections.Counter(r["status"] for r in rows)),
        valid_RTK_positions=len(valid),FIX_positions=len(fixed),FIX_rate_among_pairs=len(fixed)/len(pairs),
        decision_counts=dict(decisions),motion_checks=checks,all_FLOAT_positions_residuals_covariance_satellites_and_resets_unchanged=True,
        all_Doppler_velocities_and_covariances_identical_to_previous_carrier_lock_diagnostic=True)


def scores(rows,references,city,score):
    return score([{**r,"status":"FLOAT"} for r in rows],references,city)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--dataset-root",type=Path,required=True);parser.add_argument("--binary",type=Path,required=True);parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args();repo=Path(__file__).resolve().parents[1]
    previous=repo/"docs/comparison/2026-10-21-carrier-lock-velocity";float_archive=repo/"docs/comparison/2026-10-21-gps-galileo-float"
    prior=json.loads((previous/"predeclared-inputs.json").read_text())
    helpers=runpy.run_path(str(repo/"scripts/run_native_ppc_comparison.py"));sha,save,retain=(helpers[k] for k in ["sha","save","retain"])
    score=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))["score"]
    dataset=args.dataset_root.resolve();binary=args.binary.resolve();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    sources=list((repo/"src").rglob("*.rs"))+[repo/p for p in ["scripts/run_ppc_kinematic_fix.py","scripts/run_ppc_rtk_float.py","scripts/run_native_ppc_comparison.py","examples/ppc_rtk_kinematic_fix.rs","tests/rtk_kinematic_fix.rs","tests/rtk_fix.rs","tests/rtk_multisystem.rs","docs/kinematic-fix-policy.md","Cargo.toml","Cargo.lock","rust-toolchain.toml"]]
    inputs={}
    for case in prior["cases"]:
        city,iteration=case.split("-")
        for name in ["rover.obs","base.obs","base.nav","reference.csv"]:
            p=dataset/city/iteration/name;h=next(h for k,h in prior["input_sha256"].items() if k.endswith(f"/{city}/{iteration}/{name}"));assert sha(p)==h;inputs[str(p)]=h
    profiles={}
    for base_profile in ["galileo","galileo_slips"]:
        for gates in ["selected","all"]:
            for cap in [1,2,4]:profiles[base_profile+"_"+gates+"_std"+str(cap)]=dict(base_profile=base_profile,mode=gates+("-slips" if base_profile.endswith("_slips") else ""),motion_std_mps=cap)
    manifest=dict(base_commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=repo,text=True).strip(),source_sha256={str(p.relative_to(repo)):sha(p) for p in sources},input_sha256=inputs,binary_sha256=sha(binary),
        cases=prior["cases"],base_ecef_m=prior["base_ecef_m"],profiles=profiles,motion_max_gap_s=2,displacement_residual_m=.1,
        velocity="GPS/Galileo independent clocks; finite carrier and clear LLI; no covariance/position fusion",
        full_set_only=True,partial_AR=False,integer_hold=False,default_static_policy_changed=False,kinematic_default_motion_std_mps=1,
        predeclared_motion_std_sensitivity_after_shape_only_smoke=[1,2,4],reference_read_in_shape_smoke=False,numeric_local_integer_thresholds_changed=False,
        reference_used_for_configuration=False,held_out_validation=False,native_rerun=False,practical_PPC_parity_accepted=False,
        prior_manifest_sha256=sha(previous/"predeclared-inputs.json"))
    save(out/"predeclared-inputs.json",manifest);report=dict(practical_PPC_parity_accepted=False,cases={})
    for case in manifest["cases"]:
        city,iteration=case.split("-");root=dataset/city/iteration;report["cases"][case]={};references={}
        with (root/"reference.csv").open() as stream:
            for raw in csv.DictReader(stream):
                r={k.strip():v.strip() for k,v in raw.items()};references[(int(r["GPS Week"]),Decimal(r["GPS TOW (s)"]))]=r
        controls={}
        for base in ["galileo","galileo_slips"]:
            with gzip.open(float_archive/(case+"."+base+".stdout.csv.gz"),"rt") as stream:f=list(csv.DictReader(stream))
            with gzip.open(previous/(case+"."+base+".mixed_locked.stdout.csv.gz"),"rt") as stream:v=list(csv.DictReader(stream))
            controls[base]=(f,v)
        for profile,config in profiles.items():
            assert all(sha(repo/p)==h for p,h in manifest["source_sha256"].items()) and sha(binary)==manifest["binary_sha256"]
            command=[str(binary),str(root),config["mode"],",".join(map(str,manifest["base_ecef_m"][city])),str(config["motion_std_mps"])]
            r=subprocess.run(command,capture_output=True,cwd=repo);retain(out,case+"."+profile,command,r.stdout,r.stderr,r.returncode);assert r.returncode==0,r.stderr.decode()
            decoded,valid,fixed,metrics=summary(r.stdout,*controls[config["base_profile"]],config["motion_std_mps"],config["mode"].startswith("all"))
            metrics["all_RTK_reference_scores"]=scores(valid,references,city,score);metrics["FIX_reference_scores"]=scores(fixed,references,city,score)
            report["cases"][case][profile]=metrics;save(out/"report.json",report)
            print(case,profile,metrics["FIX_positions"],metrics["decision_counts"],flush=True)
    assert all(sha(Path(p))==h for p,h in inputs.items())


if __name__=="__main__":main()
