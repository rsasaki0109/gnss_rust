"""Freeze rover Doppler diagnostics at immutable PPC FLOAT snapshot positions."""
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


def statistics(values):
    if not values:return None
    ordered=sorted(values)
    def quantile(f):
        p=(len(ordered)-1)*f;i=int(p);return ordered[i]+(ordered[min(i+1,len(ordered)-1)]-ordered[i])*(p-i)
    return dict(epochs=len(values),rms=math.sqrt(sum(v*v for v in values)/len(values)),p50=quantile(.5),p95=quantile(.95),maximum=max(values))


def summary(data,control):
    rows=list(csv.DictReader(io.StringIO(data.decode())))
    assert len(rows)==len(control)
    counts=collections.Counter();failures=collections.Counter();std=[];speed=[]
    for r,b in zip(rows,control):
        assert {k:r[k] for k in b}==b
        if r["status"]!="FLOAT":assert not any(r[k] for k in ["velocity_status","vx_mps","vy_mps","vz_mps"]);continue
        counts[r["velocity_status"]]+=1
        if r["velocity_status"]!="VELOCITY":failures[r["velocity_reason"].split(" {")[0]]+=1;continue
        assert all(math.isfinite(float(r[k])) for k in ["vx_mps","vy_mps","vz_mps","vqxx","vqyy","vqzz","velocity_chi_square_per_dof"])
        assert 0<=float(r["velocity_chi_square_per_dof"])<=25
        std.append(max(float(r[k]) for k in ["vqxx","vqyy","vqzz"])**.5)
        speed.append(math.sqrt(sum(float(r[k])**2 for k in ["vx_mps","vy_mps","vz_mps"])))
        assert len(r["velocity_used_satellites"].split(";"))==int(r["velocity_prepared_count"])>=4
    return dict(all_previous_29_columns_identical=True,velocity_counts=dict(counts),velocity_failures=dict(failures),
        maximum_axis_velocity_std_mps=statistics(std),speed_mps=statistics(speed))


def reference_velocities(path,city,targets):
    with path.open() as stream:raws=[{k.strip():v.strip() for k,v in r.items()} for r in csv.DictReader(stream)]
    keys=[(int(r["GPS Week"]),Decimal(r["GPS TOW (s)"])) for r in raws]
    times=[Decimal(k[0])*604800+k[1] for k in keys]
    offsets=[];geometries=[]
    for r in raws:
        xyz=[float(r[k]) for k in ["ECEF X (m)","ECEF Y (m)","ECEF Z (m)"]]
        positions,lat,lon=targets(r,city)
        offsets.append({name:[v-x for v,x in zip(p,xyz)] for name,p in positions.items()})
        geometries.append((lat,lon))
    result={}
    for i,(r,key) in enumerate(zip(raws,keys)):
        lat,lon=geometries[i];e,n,u=[float(r[k]) for k in ["East Velocity (m/s)","North Velocity (m/s)","Up Velocity (m/s)"]]
        velocity=[-math.sin(lon)*e-math.sin(lat)*math.cos(lon)*n+math.cos(lat)*math.cos(lon)*u,
                  math.cos(lon)*e-math.sin(lat)*math.sin(lon)*n+math.cos(lat)*math.sin(lon)*u,
                  math.cos(lat)*n+math.sin(lat)*u]
        hypotheses={"reference_is_antenna":velocity}
        lo=max(0,i-1);hi=min(len(raws)-1,i+1);dt=float(times[hi]-times[lo])
        if hi>lo and 0<dt<=1.0:
            for name in offsets[i]:
                if name=="reference_is_antenna":continue
                derivative=[(offsets[hi][name][j]-offsets[lo][name][j])/dt for j in range(3)]
                hypotheses[name]=[a+b for a,b in zip(velocity,derivative)]
        result[key]=(hypotheses,lat,lon)
    return result


def scores(rows,references):
    errors=collections.defaultdict(lambda:collections.defaultdict(list));missing=0
    for r in rows:
        if r["velocity_status"]!="VELOCITY":continue
        key=(int(r["week"]),Decimal(r["tow"]))
        if key not in references:missing+=1;continue
        hypotheses,lat,lon=references[key];v=[float(r[k]) for k in ["vx_mps","vy_mps","vz_mps"]]
        for name,target in hypotheses.items():
            x,y,z=[a-b for a,b in zip(v,target)]
            e=-math.sin(lon)*x+math.cos(lon)*y
            n=-math.sin(lat)*math.cos(lon)*x-math.sin(lat)*math.sin(lon)*y+math.cos(lat)*z
            u=math.cos(lat)*math.cos(lon)*x+math.cos(lat)*math.sin(lon)*y+math.sin(lat)*z
            errors[name]["horizontal_error_mps"].append(math.hypot(e,n))
            errors[name]["vertical_error_mps"].append(abs(u));errors[name]["3d_error_mps"].append(math.sqrt(x*x+y*y+z*z))
    return dict(valid_velocities_without_reference=missing,hypotheses={name:{k:statistics(v) for k,v in components.items()} for name,components in errors.items()})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root",type=Path,required=True);parser.add_argument("--binary",type=Path,required=True);parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args();repo=Path(__file__).resolve().parents[1]
    previous=repo/"docs/comparison/2026-10-21-mixed-candidate-validation"
    prior=json.loads((previous/"predeclared-inputs.json").read_text())
    helpers=runpy.run_path(str(repo/"scripts/run_native_ppc_comparison.py"));sha,save,retain=(helpers[k] for k in ["sha","save","retain"])
    targets=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))["reference_targets"]
    dataset=args.dataset_root.resolve();binary=args.binary.resolve();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    sources=list((repo/"src").rglob("*.rs"))+[repo/p for p in ["scripts/run_ppc_rover_velocity.py","scripts/run_ppc_rtk_float.py","scripts/run_native_ppc_comparison.py","examples/ppc_rtk_candidates.rs","tests/velocity.rs","Cargo.toml","Cargo.lock","rust-toolchain.toml"]]
    inputs={}
    for case in prior["cases"]:
        city,iteration=case.split("-")
        for name in ["rover.obs","base.obs","base.nav","reference.csv"]:
            p=dataset/city/iteration/name;expected=next(h for k,h in prior["input_sha256"].items() if k.endswith(f"/{city}/{iteration}/{name}"));assert sha(p)==expected;inputs[str(p)]=expected
    manifest=dict(base_commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=repo,text=True).strip(),source_sha256={str(p.relative_to(repo)):sha(p) for p in sources},input_sha256=inputs,binary_sha256=sha(binary),
        cases=prior["cases"],base_ecef_m=prior["base_ecef_m"],profiles={"galileo":"galileo-velocity","galileo_slips":"galileo-slips-velocity"},
        velocity_signal="GPS 1C only; at current FLOAT ECEF position",velocity_sigma_mps=.5,elevation_degrees=10,max_chi_square_per_dof=25,max_unscaled_velocity_variance_m2s2=400,
        no_velocity_at_unmatched_or_failed_FLOAT_epochs=True,velocity_fused=False,motion_or_FIX_policy_evaluated=False,
        reference_used_for_solver_configuration=False,reference_velocity="provided ENU; add finite-difference lever-arm-offset rate for each IMU hypothesis (central, endpoints one-sided, <=1s span)",
        native_rerun=False,practical_PPC_parity_accepted=False,prior_manifest_sha256=sha(previous/"predeclared-inputs.json"))
    save(out/"predeclared-inputs.json",manifest);report=dict(practical_PPC_parity_accepted=False,cases={})
    for case in manifest["cases"]:
        city,iteration=case.split("-");root=dataset/city/iteration;report["cases"][case]={}
        references=reference_velocities(root/"reference.csv",city,targets)
        for profile,mode in manifest["profiles"].items():
            assert all(sha(repo/p)==h for p,h in manifest["source_sha256"].items()) and sha(binary)==manifest["binary_sha256"]
            command=[str(binary),str(root),mode,",".join(map(str,manifest["base_ecef_m"][city]))]
            r=subprocess.run(command,capture_output=True,cwd=repo);retain(out,case+"."+profile,command,r.stdout,r.stderr,r.returncode);assert r.returncode==0,r.stderr.decode()
            with gzip.open(previous/(case+"."+profile+".diagnostic.stdout.csv.gz"),"rt") as stream:control=list(csv.DictReader(stream))
            metrics=summary(r.stdout,control);rows=list(csv.DictReader(io.StringIO(r.stdout.decode())));metrics["reference_velocity_scores"]=scores(rows,references)
            report["cases"][case][profile]=metrics;save(out/"report.json",report)
            print(case,profile,metrics["velocity_counts"],metrics["velocity_failures"],flush=True)
    assert all(sha(Path(p))==h for p,h in inputs.items())


if __name__=="__main__":main()
