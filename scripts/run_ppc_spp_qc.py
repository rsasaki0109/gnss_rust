"""Freeze mixed-clock SPP MAD/RAIM refits with exact native-variance controls."""
import argparse
import collections
import csv
from decimal import Decimal
import gzip
import json
import math
from pathlib import Path
import runpy
import subprocess


def seed_summary(rows,mixed):
    used=collections.Counter();clocks=collections.Counter();counts=collections.Counter()
    for r in rows:
        if r["status"]!="FLOAT":continue
        counts[r["seed_source"]]+=1
        if mixed:
            assert r["spp_scope"]=="GpsGalileo"
            if r["seed_source"]=="SinglePoint":
                g,e=bool(r["spp_gps_clock_bias_s"]),bool(r["spp_galileo_clock_bias_s"]);assert g or e
                assert int(r["spp_satellites"])>=3+g+e and not r["spp_failure"]
                for k in ["spp_gps_clock_bias_s","spp_galileo_clock_bias_s"]:
                    if r[k]:assert math.isfinite(float(r[k]))
                clocks[("G" if g else "")+("E" if e else "")]+=1;used[r["spp_satellites"]]+=1
            else:assert r["spp_satellites"]=="0" and not r["spp_gps_clock_bias_s"] and not r["spp_galileo_clock_bias_s"] and r["spp_failure"]
    return dict(accepted_seed_sources=dict(counts),successful_seed_clock_groups=dict(clocks),successful_seed_satellite_counts=dict(used),scope="accepted FLOAT epochs only; not all attempted SPP epochs")


def qc_summary(rows):
    branches=collections.Counter();fallbacks=collections.Counter();total=collections.Counter()
    for r in rows:
        if r["status"]!="FLOAT":continue
        if r["seed_source"]!="SinglePoint":
            fallbacks[r["spp_failure"]]+=1
            assert all(not r[k] for k in ["spp_pre_qc_satellites","spp_pre_qc_rms_m","spp_final_rms_m","spp_outlier_rejections","spp_raim_attempts","spp_raim_rejections","spp_independent_clock_gdop"])
            continue
        pre=int(r["spp_pre_qc_satellites"]);mad=int(r["spp_outlier_rejections"]);raim=int(r["spp_raim_rejections"]);attempts=int(r["spp_raim_attempts"])
        rms0=float(r["spp_pre_qc_rms_m"]);rms1=float(r["spp_final_rms_m"]);gdop=float(r["spp_independent_clock_gdop"])
        assert pre>=int(r["spp_satellites"])+mad+raim and mad>=0 and raim in {0,1} and attempts>=0 and not(mad and raim)
        assert all(math.isfinite(v) for v in [rms0,rms1,gdop]) and rms0>=0 and rms1>=0 and 0<gdop<=50
        if mad:assert attempts==0
        if raim:assert attempts==pre and pre>=5+bool(r["spp_gps_clock_bias_s"])+bool(r["spp_galileo_clock_bias_s"]) and rms0-rms1+1e-9>=max(1,rms0*0.25)
        branches["MAD_refit" if mad else "RAIM_refit" if raim else "baseline"]+=1
        total["outlier_rejected_satellites"]+=mad;total["RAIM_attempts"]+=attempts;total["RAIM_rejected_satellites"]+=raim
    return dict(accepted_FLOAT_branches=dict(branches),accepted_FLOAT_fallback_reasons=dict(fallbacks),accepted_FLOAT_totals=dict(total),scope="accepted FLOAT epochs only; terminal/early failures may have no seed diagnostics")


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--dataset-root",type=Path,required=True);parser.add_argument("--binary",type=Path,required=True);parser.add_argument("--output-dir",type=Path,required=True);args=parser.parse_args()
    repo=Path(__file__).resolve().parents[1];previous=repo/"docs/comparison/2026-10-21-spp-variance";prior=json.loads((previous/"predeclared-inputs.json").read_text())
    scoring=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"));inspect,score,sha=(scoring[k] for k in ["inspect","score","sha"]);mixed_summary=runpy.run_path(str(repo/"scripts/run_ppc_galileo.py"))["mixed_summary"]
    helpers=runpy.run_path(str(repo/"scripts/run_native_ppc_comparison.py"));save,retain=(helpers[k] for k in ["save","retain"])
    dataset=args.dataset_root.resolve();binary=args.binary.resolve();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    sources=list((repo/"src").rglob("*.rs"))+[repo/p for p in ["scripts/run_ppc_spp_qc.py","scripts/run_ppc_galileo.py","scripts/run_ppc_rtk_float.py","scripts/run_native_ppc_comparison.py","examples/ppc_rtk_float.rs","tests/spp_multiclock.rs","tests/spp_quality.rs","tests/spp_qc.rs","scripts/upstream_spp_quality.cpp","tests/fixtures/upstream_spp_quality.csv","tests/spp.rs","tests/velocity.rs","tests/multiclock_velocity.rs","tests/fixtures/upstream_velocity.csv","tests/fixtures/upstream_multiclock_velocity.csv","scripts/upstream_rtk_mixed.cpp","scripts/upstream_rtk_float.cpp","tests/fixtures/synthetic_rtk_mixed_rover.obs","tests/fixtures/synthetic_rtk_mixed_base.obs","tests/fixtures/synthetic_rtk_mixed.nav","tests/fixtures/synthetic_rtk_mixed_truth.csv","docs/spp-measurement-variance.md","docs/spp-quality-control.md","Cargo.toml","Cargo.lock","rust-toolchain.toml"]]
    inputs={}
    for case in prior["cases"]:
        city,iteration=case.split("-")
        for name in ["rover.obs","base.obs","base.nav","reference.csv"]:
            p=dataset/city/iteration/name;h=next(v for k,v in prior["input_sha256"].items() if k.endswith(f"/{city}/{iteration}/{name}"));assert sha(p)==h;inputs[str(p)]=h
    profiles={}
    for old in ["galileo_spp_ggspp_nativevar","galileo_spp_slips_ggspp_nativevar","galileo_spp_slips_gate_ggspp_nativevar"]:
        for enabled in [False,True]:profiles[old+("_sppqc" if enabled else "")]=dict(control_profile=old,mode=prior["profiles"][old]["mode"]+("-sppqc" if enabled else ""),mixed_seed=True,quality_control=enabled)
    m=dict(base_commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=repo,text=True).strip(),source_sha256={str(p.relative_to(repo)):sha(p) for p in sources},input_sha256=inputs,binary_sha256=sha(binary),cases=prior["cases"],base_ecef_m=prior["base_ecef_m"],profiles=profiles,
        independent_receiver_clocks=True,variance_SNR_reference_dbhz=45,variance_iono_scale=0.5,variance_trop_scale=0.1,default_variance_model_changed=False,MAD_refit_connected=True,RAIM_refits_connected=True,recursive_QC=False,refit_initialization="baseline position and all system clocks (native stored MAD seed / empty local ISB differ)",geometry_clock_model="independent systems (native uses common clock)",quality_max_gdop=50,default_quality_control_changed=False,legacy_primary_group_delay=True,phase_required_for_SPP=False,SPP_covariance_fused=False,existing_thresholds_changed=False,default_seed_scope_changed=False,native_whole_rerun=False,held_out_validation=False,truth_used_by_solver=False,practical_PPC_parity_accepted=False,prior_manifest_sha256=sha(previous/"predeclared-inputs.json"))
    save(out/"predeclared-inputs.json",m);report=dict(practical_PPC_parity_accepted=False,cases={})
    for case in m["cases"]:
        city,iteration=case.split("-");root=dataset/city/iteration;report["cases"][case]={};ledger=None
        with (root/"reference.csv").open() as f:refs=[{k.strip():v.strip() for k,v in r.items()} for r in csv.DictReader(f)]
        references={(int(r["GPS Week"]),Decimal(r["GPS TOW (s)"])):r for r in refs}
        for profile,c in profiles.items():
            assert all(sha(repo/p)==h for p,h in m["source_sha256"].items()) and sha(binary)==m["binary_sha256"]
            command=[str(binary),str(root),c["mode"],",".join(map(str,m["base_ecef_m"][city]))]
            r=subprocess.run(command,capture_output=True,cwd=repo);retain(out,case+"."+profile,command,r.stdout,r.stderr,r.returncode);assert r.returncode==0,r.stderr.decode()
            rows,metrics=inspect(r.stdout);keys=[(r["kind"],r["week"],r["tow"],r["side"]) for r in rows];assert ledger is None or ledger==keys;ledger=keys
            if not c["quality_control"]:
                with gzip.open(previous/(case+"."+profile+".stdout.csv.gz"),"rt") as f:assert rows==list(csv.DictReader(f))
                metrics["all_columns_identical_to_prior_native_variance_control"]=True
            if c["quality_control"]:metrics["SPP_QC"]=qc_summary(rows)
            metrics["mixed_system"]=mixed_summary(rows);metrics["SPP_seed"]=seed_summary(rows,c["mixed_seed"]);metrics["reference_scores"]=score(rows,references,city)
            seeds=[{**r,"x_m":r["seed_x_m"],"y_m":r["seed_y_m"],"z_m":r["seed_z_m"]} for r in rows if r["status"]=="FLOAT" and r["seed_source"]=="SinglePoint"]
            metrics["accepted_SPP_seed_reference_scores"]=score(seeds,references,city)
            report["cases"][case][profile]=metrics;save(out/"report.json",report);print(case,profile,metrics["status_counts"],metrics["SPP_seed"]["accepted_seed_sources"],flush=True)
    assert all(sha(Path(p))==h for p,h in inputs.items())


if __name__=="__main__":main()
