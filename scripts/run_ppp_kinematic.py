"""Freeze native-component/synthetic moving PPP and unchanged static solver controls."""
import argparse
import collections
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import runpy
import subprocess


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rows(data):return list(csv.DictReader(io.StringIO(data.decode())))
def moving_metrics(data,truth,moving):
    rr=rows(data);assert len(rr)==64
    statuses=collections.Counter(r["status"] for r in rr);assert set(statuses)<={"PPP_FLOAT","INVALID","EVENT"}
    errors=[];velocities=[];exceeds=[]
    for i,r in enumerate(rr):
        t=truth[i];assert (int(r["week"]),Decimal(r["tow"]))==(int(t["week"]),Decimal(t["tow"]))
        if r["status"]!="PPP_FLOAT":assert all(not r[k] for k in ["x_m","y_m","z_m","vx_mps","vy_mps","vz_mps"]);continue
        assert not r["error"] and all(math.isfinite(float(r[k])) for k in ["x_m","y_m","z_m","clock_s","zenith_m","position_covariance_trace_m2","code_rms_m","phase_rms_m","nis"])
        error=math.sqrt(sum((float(r[k])-float(t[k]))**2 for k in ["x_m","y_m","z_m"]));errors.append((i,error))
        if error>=.002:exceeds.append(dict(epoch=i,week=r["week"],tow=r["tow"],error_m=error))
        if moving:
            assert all(math.isfinite(float(r[k])) for k in ["vx_mps","vy_mps","vz_mps","velocity_covariance_trace_m2_s2"])
            v=math.sqrt(sum((float(r[k])-float(t[k]))**2 for k in ["vx_mps","vy_mps","vz_mps"]));velocities.append((i,v))
        else:assert all(not r[k] for k in ["vx_mps","vy_mps","vz_mps","velocity_covariance_trace_m2_s2"])
    def summary(values):return dict(epochs=len(values),rms=math.sqrt(sum(v*v for _,v in values)/len(values)),maximum=max(v for _,v in values)) if values else None
    return dict(input_epochs=64,status_counts=dict(statuses),position_3D_m=summary(errors),position_after_first_four_epochs_m=summary([(i,v) for i,v in errors if i>=4]),velocity_3D_mps=summary(velocities),velocity_after_first_two_epochs_mps=summary([(i,v) for i,v in velocities if i>=2]),full_trace_2mm_target_met=len(errors)==64 and not exceeds,two_mm_exceedances=exceeds,terminal_error_reasons=dict(collections.Counter(r["error"] for r in rr if r["status"]=="INVALID")),arc_reset_total=sum(int(r["arc_resets"]) for r in rr if r["status"]=="PPP_FLOAT"))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--moving-binary",type=Path,required=True);p.add_argument("--static-binary",type=Path,required=True);p.add_argument("--old-static-root",type=Path,required=True);p.add_argument("--output-dir",type=Path,required=True);a=p.parse_args()
    repo=Path(__file__).resolve().parents[1];out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True);assert not (out/"predeclared-inputs.json").exists()
    moving=a.moving_binary.resolve();static=a.static_binary.resolve();oldroot=a.old_static_root.resolve();old=oldroot/"ppp";build=json.loads((oldroot/"build.json").read_text());assert sha(old)==build["binary_sha256"]
    sources=list((repo/"src").rglob("*.rs"))+[repo/f for f in ["examples/ppp.rs","examples/ppp_kinematic.rs","scripts/run_ppp_kinematic.py","scripts/build_ppp_static_control.py","scripts/upstream_ppp_kinematic.cpp","scripts/upstream_precise_transmit.cpp","scripts/upstream_precise.cpp","scripts/upstream_navigation.cpp","scripts/run_native_ppc_comparison.py","tests/ppp.rs","tests/ppp_kinematic.rs","docs/kinematic-ppp.md","Cargo.toml","Cargo.lock","rust-toolchain.toml"]]+list((repo/"tests/fixtures").glob("*ppp*"))
    profiles={}
    for model in ["static","kinematic"]:
        for variant in ["clean","noisy"]:profiles[model+"_"+variant]=dict(arguments=[model,variant],moving=model=="kinematic",velocity_variance=100 if model=="kinematic" else None)
    for variance in [10,1000]:
        for variant in ["clean","noisy"]:profiles[f"kinematic_{variant}_v{variance}"]=dict(arguments=["kinematic",variant,str(variance)],moving=True,velocity_variance=variance)
    for variant in ["unrounded-phase","unrounded-code","unrounded-both"]:profiles["kinematic_"+variant.replace("-","_")]=dict(arguments=["kinematic",variant],moving=True,velocity_variance=100)
    controls={}
    for mode,args in [("plain",[]),("windup",["windup"]),("antenna",["antenna"]),("satellite",["satellite"]),("legacy_tide",["legacy-tide"])]:
        for noise in [False,True]:controls[mode+("_noisy" if noise else "_clean")]=(["noisy"] if noise else [])+args
    manifest=dict(base_commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=repo,text=True).strip(),source_sha256={str(s.relative_to(repo)):sha(s) for s in sources},binary_sha256=dict(moving=sha(moving),static=sha(static),old_static=sha(old)),moving_profiles=profiles,static_profiles=controls,old_static_build_sha256=sha(oldroot/"build.json"),original_full_trace_clean_target_m=.002,default_position_model_changed=False,native_predict_component_calls=24,whole_native_PPP_executed=False,real_PPC_PPP_executed=False,truth_used_by_solver=False,existing_thresholds_changed=False,practical_parity_accepted=False)
    helper=runpy.run_path(str(repo/"scripts/run_native_ppc_comparison.py"));save,retain=(helper[k] for k in ["save","retain"])
    save(out/"predeclared-inputs.json",manifest);(out/"old-static-build.json").write_bytes((oldroot/"build.json").read_bytes())
    for check in build["checks"]:(out/("old-static-"+check["log"])).write_bytes((oldroot/check["log"]).read_bytes())
    (out/"old-static-ppp.rs.gz").write_bytes(gzip.compress((oldroot/"src/ppp.rs").read_bytes(),mtime=0))
    truth=rows((repo/"tests/fixtures/synthetic_ppp_kinematic_truth.csv").read_bytes());report=dict(moving={},static_controls={},whole_native_PPP_executed=False,real_PPC_PPP_executed=False,practical_parity_accepted=False)
    def execute(stem,binary,arguments):
        assert all(sha(repo/p)==h for p,h in manifest["source_sha256"].items())
        command=[str(binary)]+arguments;r=subprocess.run(command,capture_output=True,cwd=repo);retain(out,stem,command,r.stdout,r.stderr,r.returncode);assert r.returncode==0,r.stderr.decode();return r.stdout
    for name,c in profiles.items():
        data=execute(name,moving,c["arguments"]);metrics=moving_metrics(data,truth,c["moving"]);report["moving"][name]=metrics;save(out/"report.json",report);print(name,metrics["status_counts"],metrics["position_3D_m"],flush=True)
    for name,arguments in controls.items():
        prior=execute("static-old-"+name,old,arguments);current=execute("static-new-"+name,static,arguments);assert prior==current,(name,"static CSV mismatch");rr=rows(current);assert len(rr)==64 and all(r["status"]=="PPP_FLOAT" for r in rr)
        report["static_controls"][name]=dict(epochs=64,all_CSV_bytes_identical=True,csv_sha256=hashlib.sha256(current).hexdigest());save(out/"report.json",report);print('static-control',name,'all bytes identical',flush=True)
    assert sha(moving)==manifest["binary_sha256"]["moving"] and sha(static)==manifest["binary_sha256"]["static"] and sha(old)==manifest["binary_sha256"]["old_static"]


if __name__=="__main__":main()
