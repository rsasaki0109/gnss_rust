"""Freeze Doppler carrier-phase/LLI admission checks; no motion or FIX fusion."""
import argparse
import csv
import gzip
import io
import json
from pathlib import Path
import runpy
import subprocess


def comparison(a,b,references,score,multi,a_name,b_name):
    c=multi["compare"](a,b,references,score)
    return dict(a_condition=a_name,b_condition=b_name,common_velocity_epochs=c["common_velocity_epochs"],
                a_only=c["GPS_only"],b_only=c["mixed_only"],common_reference_scores={a_name:c["common_reference_scores"]["GPS"],b_name:c["common_reference_scores"]["mixed"]})


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--dataset-root",type=Path,required=True);parser.add_argument("--binary",type=Path,required=True);parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args();repo=Path(__file__).resolve().parents[1]
    previous=repo/"docs/comparison/2026-10-21-multiclock-velocity";candidates=repo/"docs/comparison/2026-10-21-mixed-candidate-validation"
    prior=json.loads((previous/"predeclared-inputs.json").read_text())
    helpers=runpy.run_path(str(repo/"scripts/run_native_ppc_comparison.py"));sha,save,retain=(helpers[k] for k in ["sha","save","retain"])
    rover=runpy.run_path(str(repo/"scripts/run_ppc_rover_velocity.py"));multi=runpy.run_path(str(repo/"scripts/run_ppc_multiclock_velocity.py"));targets=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))["reference_targets"]
    dataset=args.dataset_root.resolve();binary=args.binary.resolve();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    assert not (out/"predeclared-inputs.json").exists()
    sources=list((repo/"src").rglob("*.rs"))+[repo/p for p in ["scripts/run_ppc_carrier_lock_velocity.py","scripts/run_ppc_rover_velocity.py","scripts/run_ppc_multiclock_velocity.py","scripts/run_ppc_rtk_float.py","scripts/run_native_ppc_comparison.py","examples/ppc_rtk_candidates.rs","examples/ppc_phase_doppler_probe.rs","tests/velocity.rs","tests/multiclock_velocity.rs","Cargo.toml","Cargo.lock","rust-toolchain.toml"]]
    inputs={}
    for case in prior["cases"]:
        city,iteration=case.split("-")
        for name in ["rover.obs","base.obs","base.nav","reference.csv"]:
            p=dataset/city/iteration/name;h=next(h for k,h in prior["input_sha256"].items() if k.endswith(f"/{city}/{iteration}/{name}"));assert sha(p)==h;inputs[str(p)]=h
    manifest=dict(base_commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=repo,text=True).strip(),source_sha256={str(p.relative_to(repo)):sha(p) for p in sources},input_sha256=inputs,binary_sha256=sha(binary),
        cases=prior["cases"],base_ecef_m=prior["base_ecef_m"],profiles=prior["profiles"],conditions=["GPS","mixed","GPS_locked","mixed_locked"],
        changed_setting="require_carrier_lock false -> true: finite carrier phase and no LLI loss/half-cycle; selected tracking unchanged",
        carrier_values_used_in_velocity_equation=False,minimum_lock_or_temporal_continuity_inferred=False,
        admission_chosen_after_post_run_raw_diagnostic=True,numeric_thresholds_chosen_using_reference=False,held_out_validation=False,
        numeric_thresholds_changed=False,default_changed=False,velocity_fused=False,motion_or_FIX_policy_evaluated=False,native_rerun=False,practical_PPC_parity_accepted=False,
        prior_manifest_sha256=sha(previous/"predeclared-inputs.json"))
    save(out/"predeclared-inputs.json",manifest);report=dict(practical_PPC_parity_accepted=False,cases={})
    for case in manifest["cases"]:
        city,iteration=case.split("-");root=dataset/city/iteration;report["cases"][case]={};references=rover["reference_velocities"](root/"reference.csv",city,targets)
        for profile,mode in manifest["profiles"].items():
            with gzip.open(candidates/(case+"."+profile+".diagnostic.stdout.csv.gz"),"rt") as stream:control=list(csv.DictReader(stream))
            all_rows={};report["cases"][case][profile]={}
            for condition in manifest["conditions"]:
                assert all(sha(repo/p)==h for p,h in manifest["source_sha256"].items()) and sha(binary)==manifest["binary_sha256"]
                selected=mode+("-mixed" if condition.startswith("mixed") else "")+("-locked" if condition.endswith("_locked") else "")
                command=[str(binary),str(root),selected,",".join(map(str,manifest["base_ecef_m"][city]))]
                r=subprocess.run(command,capture_output=True,cwd=repo);retain(out,case+"."+profile+"."+condition,command,r.stdout,r.stderr,r.returncode);assert r.returncode==0,r.stderr.decode()
                rows=list(csv.DictReader(io.StringIO(r.stdout.decode())));metrics=rover["summary"](r.stdout,control)
                if condition in {"GPS","mixed"}:
                    with gzip.open(previous/(case+"."+profile+"."+condition+".stdout.csv.gz"),"rt") as stream:assert rows==list(csv.DictReader(stream))
                    metrics["all_columns_identical_to_prior_control"]=True
                if condition.startswith("mixed"):metrics["mixed"]=multi["mixed_summary"](rows)
                metrics["reference_velocity_scores"]=rover["scores"](rows,references)
                report["cases"][case][profile][condition]=metrics;all_rows[condition]=rows;save(out/"report.json",report)
                print(case,profile,condition,metrics["velocity_counts"],flush=True)
            report["cases"][case][profile]["common_epochs"]={}
            for a,b in [("GPS","GPS_locked"),("mixed","mixed_locked"),("GPS_locked","mixed_locked")]:
                report["cases"][case][profile]["common_epochs"][a+"_vs_"+b]=comparison(all_rows[a],all_rows[b],references,rover["scores"],multi,a,b)
            save(out/"report.json",report)
    assert all(sha(Path(p))==h for p,h in inputs.items())


if __name__=="__main__":main()
