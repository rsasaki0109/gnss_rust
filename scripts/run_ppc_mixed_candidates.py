"""Freeze unchanged PPC FLOAT controls and read-only mixed integer diagnostics."""
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


def candidate_summary(data,control):
    rows=list(csv.DictReader(io.StringIO(data.decode())))
    assert len(rows)==len(control)
    counts=collections.Counter();local=collections.Counter();audit=collections.Counter();sets=0
    for r,b in zip(rows,control):
        for k in ["kind","week","tow","side","status","x_m","y_m","z_m","used_satellites","reference_satellites","ambiguity_resets"]:assert r[k]==b[k],(k,r,b)
        assert r["status"]!="FIX"
        if r["status"]!="FLOAT":
            assert all(not r[k] for k in ["conditioned_x_m","conditioned_y_m","conditioned_z_m","local_status","audit_status"])
            continue
        assert int(r["ambiguity_pairs"])==int(b["dd_phase_rows"])
        counts["FLOAT_snapshots"]+=1
        if r["search_status"]!="OK":counts["search_errors"]+=1;continue
        counts["search_successes"]+=1
        assert len(r["integers"].split(";"))==int(r["ambiguity_pairs"])
        ratio=bool(r["ratio"]) and math.isfinite(float(r["ratio"])) and float(r["ratio"])>=3.0
        counts["ratio_at_least_3"]+=ratio
        local[r["local_rejection"] or r["local_status"]]+=1
        if r["local_status"]=="PASS_LOCAL_GATES":
            counts["local_gate_pass"]+=1
            if ratio:
                counts["ratio_and_local_gate_pass"]+=1
                counts["and_at_least_four_pairs"]+=int(r["ambiguity_pairs"])>=4
                assert r["audit_status"] in {"COMPLETE","ERROR"}
                if r["audit_status"]=="COMPLETE":
                    counts["complete_reference_audits"]+=1;sets+=int(r["reference_sets"])
                    audit[r["conservative_rejection"]]+=1
                    if r["conservative_rejection"]=="PASS_LOCAL_GATES":counts["ratio_and_all_reference_gates_pass"]+=1
        else:assert not r["audit_status"]
    return dict(counts=dict(counts),local_results=dict(local),conservative_results=dict(audit),reference_sets_evaluated=sets,
                all_FLOAT_statuses_positions_satellites_and_resets_identical=True,FIX_declared=0)


def selections(rows):
    numeric=[r for r in rows if r["status"]=="FLOAT" and r["local_status"] not in {"","ERROR"}]
    passing=[r for r in numeric if r["local_status"]=="PASS_LOCAL_GATES" and r["ratio"] and math.isfinite(float(r["ratio"])) and float(r["ratio"])>=3]
    return dict(all_numerically_conditioned=numeric,ratio_and_local_gates=passing,
                ratio_and_all_reference_gates=[r for r in passing if r["conservative_rejection"]=="PASS_LOCAL_GATES"])


def candidate_scores(rows,references,city,score):
    return {name:score([{**r,"x_m":r["conditioned_x_m"],"y_m":r["conditioned_y_m"],"z_m":r["conditioned_z_m"]} for r in selected],references,city) for name,selected in selections(rows).items()}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root",type=Path,required=True)
    parser.add_argument("--float-binary",type=Path,required=True)
    parser.add_argument("--candidate-binary",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args();repo=Path(__file__).resolve().parents[1]
    previous=repo/"docs/comparison/2026-10-21-gps-galileo-float"
    prior=json.loads((previous/"predeclared-inputs.json").read_text())
    helpers=runpy.run_path(str(repo/"scripts/run_native_ppc_comparison.py"));sha,save,retain=(helpers[k] for k in ["sha","save","retain"])
    scoring=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))
    dataset=args.dataset_root.resolve();float_binary=args.float_binary.resolve();binary=args.candidate_binary.resolve();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    sources=list((repo/"src").rglob("*.rs"))+[repo/p for p in ["scripts/run_ppc_mixed_candidates.py","scripts/run_ppc_rtk_float.py","scripts/run_native_ppc_comparison.py","examples/ppc_rtk_candidates.rs","examples/ppc_rtk_float.rs","tests/rtk_multisystem.rs","tests/rtk_fix.rs","Cargo.toml","Cargo.lock","rust-toolchain.toml"]]
    inputs={}
    for case in prior["cases"]:
        city,iteration=case.split("-")
        for name in ["rover.obs","base.obs","base.nav","reference.csv"]:
            p=dataset/city/iteration/name;expected=next(h for k,h in prior["input_sha256"].items() if k.endswith(f"/{city}/{iteration}/{name}"))
            assert sha(p)==expected;inputs[str(p)]=expected
    manifest=dict(base_commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=repo,text=True).strip(),
        source_sha256={str(p.relative_to(repo)):sha(p) for p in sources},input_sha256=inputs,
        binaries_sha256={str(float_binary):sha(float_binary),str(binary):sha(binary)},cases=prior["cases"],base_ecef_m=prior["base_ecef_m"],
        controls={"gps":"kinematic","galileo":"kinematic-galileo","galileo_slips":"kinematic-galileo-slips"},
        diagnostic_modes={"galileo":"galileo","galileo_slips":"galileo-slips"},
        reference_audit_limit=4096,all_reference_combinations_required=True,
        thresholds_changed=False,FIX_policy_enabled=False,lock_confirmation_motion_history_evaluated=False,
        native_rerun=False,truth_used_for_solver_configuration=False,practical_PPC_parity_accepted=False,
        prior_manifest_sha256=sha(previous/"predeclared-inputs.json"))
    save(out/"predeclared-inputs.json",manifest);report=dict(practical_PPC_parity_accepted=False,cases={})
    for case in manifest["cases"]:
        city,iteration=case.split("-");root=dataset/city/iteration;references={};report["cases"][case]={}
        with (root/"reference.csv").open() as stream:
            for raw in csv.DictReader(stream):
                r={k.strip():v.strip() for k,v in raw.items()};references[(int(r["GPS Week"]),Decimal(r["GPS TOW (s)"]))]=r
        for profile,mode in manifest["controls"].items():
            assert all(sha(repo/p)==h for p,h in manifest["source_sha256"].items())
            assert all(sha(Path(p))==h for p,h in manifest["binaries_sha256"].items())
            command=[str(float_binary),str(root),mode,",".join(map(str,manifest["base_ecef_m"][city]))]
            result=subprocess.run(command,capture_output=True,cwd=repo);retain(out,case+"."+profile+".control",command,result.stdout,result.stderr,result.returncode);assert result.returncode==0
            control,metrics=scoring["inspect"](result.stdout)
            with gzip.open(previous/(case+"."+profile+".stdout.csv.gz"),"rt") as stream:assert control==list(csv.DictReader(stream))
            metrics["all_columns_identical_to_prior_control"]=True;report["cases"][case][profile]=dict(control=metrics)
            if profile in manifest["diagnostic_modes"]:
                command=[str(binary),str(root),manifest["diagnostic_modes"][profile],",".join(map(str,manifest["base_ecef_m"][city]))]
                result=subprocess.run(command,capture_output=True,cwd=repo);retain(out,case+"."+profile+".diagnostic",command,result.stdout,result.stderr,result.returncode);assert result.returncode==0,result.stderr.decode()
                summary=candidate_summary(result.stdout,control);decoded=list(csv.DictReader(io.StringIO(result.stdout.decode())))
                summary["candidate_reference_scores"]=candidate_scores(decoded,references,city,scoring["score"])
                report["cases"][case][profile]["diagnostic"]=summary
                print(case,profile,summary["counts"],flush=True)
            save(out/"report.json",report)
    assert all(sha(Path(p))==h for p,h in inputs.items())


if __name__=="__main__":main()
