"""Freeze optional native SPP measurement variance with unchanged mixed-clock controls."""
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


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--dataset-root",type=Path,required=True);parser.add_argument("--binary",type=Path,required=True);parser.add_argument("--output-dir",type=Path,required=True);args=parser.parse_args()
    repo=Path(__file__).resolve().parents[1];previous=repo/"docs/comparison/2026-10-21-multiclock-spp";prior=json.loads((previous/"predeclared-inputs.json").read_text())
    scoring=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"));inspect,score,sha=(scoring[k] for k in ["inspect","score","sha"]);mixed_summary=runpy.run_path(str(repo/"scripts/run_ppc_galileo.py"))["mixed_summary"]
    helpers=runpy.run_path(str(repo/"scripts/run_native_ppc_comparison.py"));save,retain=(helpers[k] for k in ["save","retain"])
    dataset=args.dataset_root.resolve();binary=args.binary.resolve();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    sources=list((repo/"src").rglob("*.rs"))+[repo/p for p in ["scripts/run_ppc_spp_variance.py","scripts/run_ppc_galileo.py","scripts/run_ppc_rtk_float.py","scripts/run_native_ppc_comparison.py","examples/ppc_rtk_float.rs","tests/spp_multiclock.rs","tests/spp_quality.rs","scripts/upstream_spp_quality.cpp","tests/fixtures/upstream_spp_quality.csv","tests/spp.rs","tests/velocity.rs","tests/multiclock_velocity.rs","tests/fixtures/upstream_velocity.csv","tests/fixtures/upstream_multiclock_velocity.csv","scripts/upstream_rtk_mixed.cpp","scripts/upstream_rtk_float.cpp","tests/fixtures/synthetic_rtk_mixed_rover.obs","tests/fixtures/synthetic_rtk_mixed_base.obs","tests/fixtures/synthetic_rtk_mixed.nav","tests/fixtures/synthetic_rtk_mixed_truth.csv","docs/spp-measurement-variance.md","Cargo.toml","Cargo.lock","rust-toolchain.toml"]]
    inputs={}
    for case in prior["cases"]:
        city,iteration=case.split("-")
        for name in ["rover.obs","base.obs","base.nav","reference.csv"]:
            p=dataset/city/iteration/name;h=next(v for k,v in prior["input_sha256"].items() if k.endswith(f"/{city}/{iteration}/{name}"));assert sha(p)==h;inputs[str(p)]=h
    profiles={}
    for old in ["galileo_spp_ggspp","galileo_spp_slips_ggspp","galileo_spp_slips_gate_ggspp"]:
        for native in [False,True]:profiles[old+("_nativevar" if native else "")]=dict(control_profile=old,mode=prior["profiles"][old]["mode"]+("-nativevar" if native else ""),mixed_seed=True,native_variance=native)
    m=dict(base_commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=repo,text=True).strip(),source_sha256={str(p.relative_to(repo)):sha(p) for p in sources},input_sha256=inputs,binary_sha256=sha(binary),cases=prior["cases"],base_ecef_m=prior["base_ecef_m"],profiles=profiles,
        independent_receiver_clocks=True,variance_SNR_reference_dbhz=45,variance_iono_scale=0.5,variance_trop_scale=0.1,default_variance_model_changed=False,robust_mask_component_only=True,RAIM_selector_component_only=True,legacy_primary_group_delay=True,phase_required_for_SPP=False,SPP_covariance_fused=False,existing_thresholds_changed=False,default_seed_scope_changed=False,native_whole_rerun=False,held_out_validation=False,truth_used_by_solver=False,practical_PPC_parity_accepted=False,prior_manifest_sha256=sha(previous/"predeclared-inputs.json"))
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
            if not c["native_variance"]:
                with gzip.open(previous/(case+"."+profile+".stdout.csv.gz"),"rt") as f:assert rows==list(csv.DictReader(f))
                metrics["all_columns_identical_to_prior_mixed_seed_control"]=True
            metrics["mixed_system"]=mixed_summary(rows);metrics["SPP_seed"]=seed_summary(rows,c["mixed_seed"]);metrics["reference_scores"]=score(rows,references,city)
            seeds=[{**r,"x_m":r["seed_x_m"],"y_m":r["seed_y_m"],"z_m":r["seed_z_m"]} for r in rows if r["status"]=="FLOAT" and r["seed_source"]=="SinglePoint"]
            metrics["accepted_SPP_seed_reference_scores"]=score(seeds,references,city)
            report["cases"][case][profile]=metrics;save(out/"report.json",report);print(case,profile,metrics["status_counts"],metrics["SPP_seed"]["accepted_seed_sources"],flush=True)
    assert all(sha(Path(p))==h for p,h in inputs.items())


if __name__=="__main__":main()
