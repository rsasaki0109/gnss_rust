"""Freeze stationary-base GPS 1C Doppler diagnostics on all six PPC runs."""
import argparse
import collections
import csv
import gzip
import io
import json
import math
from pathlib import Path
import runpy
import subprocess


def summary(data):
    rows=list(csv.DictReader(io.StringIO(data.decode())))
    good=[r for r in rows if r["status"]=="VELOCITY"]
    failures=collections.Counter(r["reason"].split(" {")[0] for r in rows if r["status"]!="VELOCITY")
    for r in rows:
        if r["status"]!="VELOCITY":
            assert not any(r[k] for k in ["vx_mps","vy_mps","vz_mps","speed_mps"])
            continue
        assert all(math.isfinite(float(r[k])) for k in ["vx_mps","vy_mps","vz_mps","speed_mps","chi_square_per_dof","velocity_covariance_trace_m2s2"])
        assert 0<=float(r["chi_square_per_dof"])<=25
        assert float(r["covariance_scale"])>=1 and float(r["max_unscaled_velocity_variance_m2s2"])<=400
        assert len(r["used_satellites"].split(";"))==int(r["prepared_count"])>=4
    speed=sorted(float(r["speed_mps"]) for r in good)
    def quantile(f):
        p=(len(speed)-1)*f;lo=int(p);return speed[lo]+(speed[min(lo+1,len(speed)-1)]-speed[lo])*(p-lo)
    return dict(epochs=len(rows),status_counts=dict(collections.Counter(r["status"] for r in rows)),failures=dict(failures),
        prepared_measured_zero_rows=sum(int(r["zero_doppler_count"] or 0) for r in rows),
        missing_code_fallback_rows=sum(int(r["missing_code_count"] or 0) for r in rows),
        stationary_base_speed_mps=dict(epochs=len(good),rms=math.sqrt(sum(v*v for v in speed)/len(speed)),p50=quantile(.5),p95=quantile(.95),maximum=max(speed)) if speed else None,
        mean_velocity_mps={k:sum(float(r[k]) for r in good)/len(good) for k in ["vx_mps","vy_mps","vz_mps"]} if good else None)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root",type=Path,required=True)
    parser.add_argument("--binary",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args();repo=Path(__file__).resolve().parents[1]
    previous=repo/"docs/comparison/2026-10-21-gps-galileo-float"
    prior=json.loads((previous/"predeclared-inputs.json").read_text());prior_report=json.loads((previous/"report.json").read_text())
    helpers=runpy.run_path(str(repo/"scripts/run_native_ppc_comparison.py"));sha,save,retain=(helpers[k] for k in ["sha","save","retain"])
    dataset=args.dataset_root.resolve();binary=args.binary.resolve();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    sources=list((repo/"src").rglob("*.rs"))+[repo/p for p in ["scripts/run_ppc_base_velocity.py","scripts/run_native_ppc_comparison.py","scripts/upstream_velocity.cpp","tests/fixtures/upstream_velocity.csv","tests/velocity.rs","examples/ppc_base_velocity.rs","Cargo.toml","Cargo.lock","rust-toolchain.toml"]]
    inputs={}
    for case in prior["cases"]:
        city,iteration=case.split("-")
        for name in ["base.obs","base.nav"]:
            p=dataset/city/iteration/name;h=next(h for k,h in prior["input_sha256"].items() if k.endswith(f"/{city}/{iteration}/{name}"))
            assert sha(p)==h;inputs[str(p)]=h
    manifest=dict(base_commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=repo,text=True).strip(),
        source_sha256={str(p.relative_to(repo)):sha(p) for p in sources},input_sha256=inputs,binary_sha256=sha(binary),
        cases=prior["cases"],base_ecef_m=prior["base_ecef_m"],signal="GPS 1C only",sigma_mps=.5,elevation_degrees=10,
        max_ephemeris_age_s=14400,min_satellites=4,max_chi_square_per_dof=25,max_unscaled_velocity_variance_m2s2=400,
        native_zero_sentinel_replaced_by_explicit_missingness=True,base_stationary_hypothesis=True,
        rover_positions_or_truth_read=False,velocity_fused_into_RTK=False,FIX_policy_evaluated=False,
        native_whole_solver_rerun=False,practical_PPC_parity_accepted=False,prior_manifest_sha256=sha(previous/"predeclared-inputs.json"))
    save(out/"predeclared-inputs.json",manifest);report=dict(practical_PPC_parity_accepted=False,cases={})
    for case in manifest["cases"]:
        city,iteration=case.split("-");root=dataset/city/iteration
        assert all(sha(repo/p)==h for p,h in manifest["source_sha256"].items()) and sha(binary)==manifest["binary_sha256"]
        command=[str(binary),str(root),",".join(map(str,manifest["base_ecef_m"][city]))]
        r=subprocess.run(command,capture_output=True,cwd=repo);retain(out,case,command,r.stdout,r.stderr,r.returncode);assert r.returncode==0,r.stderr.decode()
        metrics=summary(r.stdout);assert metrics["epochs"]==prior_report["cases"][case]["gps"]["base_epochs"]
        with gzip.open(previous/(case+".gps.stdout.csv.gz"),"rt") as stream:expected=[(row["week"],float(row["tow"]).hex()) for row in csv.DictReader(stream) if row["kind"]=="pair" or row["side"]=="base"]
        actual=[(row["week"],float(row["tow"]).hex()) for row in csv.DictReader(io.StringIO(r.stdout.decode()))];assert actual==expected
        metrics["all_base_epoch_keys_identical_to_prior_ledger"]=True;report["cases"][case]=metrics;save(out/"report.json",report)
        print(case,metrics["status_counts"],metrics["stationary_base_speed_mps"],flush=True)
    assert all(sha(Path(p))==h for p,h in inputs.items())


if __name__=="__main__":main()
