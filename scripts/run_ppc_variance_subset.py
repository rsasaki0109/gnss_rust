"""Compare native relative-variance preselection with unchanged lock-qualified controls."""
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


def rows(data):return list(csv.DictReader(io.StringIO(data.decode())))


def summarize(data,control,variance,cap,all_references):
    decoded=rows(data);assert len(decoded)==len(control);valid=[];fixed=[];decisions=collections.Counter();checks=0;partial=0;filtered=0;removed=0
    raw=["float_x_m","float_y_m","float_z_m","used_satellites","reference_satellites","ambiguity_pairs","ambiguity_resets","velocity_status","vx_mps","vy_mps","vz_mps","velocity_max_std_mps","velocity_reason","float_code_rms_m","float_phase_rms_m","float_qxx_m2","float_qyy_m2","float_qzz_m2","satellite_lock_epochs"]+["velocity_q"+a for a in ["xx","xy","xz","yx","yy","yz","zx","zy","zz"]]
    for r,b in zip(decoded,control):
        for k in ["kind","week","tow","side"]:assert r[k]==b[k]
        if not variance:assert {k:r[k] for k in b}==b
        if b["status"] not in {"FLOAT","FIX"}:
            assert r["status"]==b["status"] and not r["x_m"] and not r["ar_satellites"];continue
        assert all(r[k]==b[k] for k in raw);assert r["status"] in {"FLOAT","FIX"}
        sats=r["ar_satellites"].split(";") if r["ar_satellites"] else [];refs=r["ar_reference_satellites"].split(";") if r["ar_reference_satellites"] else []
        locks={s:int(n) for s,n in [item.split(":") for item in r["satellite_lock_epochs"].split(";")]};used=r["used_satellites"].split(";");assert set(locks)==set(used)
        groups=collections.Counter(s[0] for s in sats);assert len(refs)==len(groups) and all(n>=2 for n in groups.values()) and set(refs)<=set(sats)
        assert len(set(s[0] for s in refs))==len(refs) and int(r["ar_ambiguity_pairs"])==len(sats)-len(refs)
        qualified=[s for s in used if locks[s]>=5];counts=collections.Counter(s[0] for s in qualified);qualified=[s for s in qualified if counts[s[0]]>=2]
        expected=qualified
        if variance and qualified:
            source_refs=r["variance_source_references"].split(";");targets=r["variance_source_targets"].split(";");values=[float(v) for v in r["variance_dd_cycles2"].split(";")]
            assert targets==[s for s in qualified if s not in source_refs] and set(source_refs)|set(targets)==set(qualified) and len(values)==len(targets)
            assert len({s[0] for s in source_refs})==len(source_refs) and len(source_refs)==len({s[0] for s in qualified})
            assert all(math.isfinite(v) and v>=0 for v in values)
            median=sorted(values)[len(values)//2];threshold=max(median*10,1e-4);assert float(r["variance_median_cycles2"])==median and float(r["variance_threshold_cycles2"])==threshold
            indices=[i for i,v in enumerate(values) if v<=threshold]
            if len(indices)<4:indices=[]
            assert r["variance_selected_indices"]==";".join(map(str,indices))
            applied=bool(indices) and len(indices)<len(values);assert r["variance_applied"]==str(applied).lower()
            if applied:
                selected_targets=[targets[i] for i in indices];systems={s[0] for s in selected_targets};kept=set(selected_targets)|{s for s in source_refs if s[0] in systems}
                expected=[s for s in qualified if s in kept];filtered+=1;removed+=len(values)-len(indices)
        else:assert not any(r[k] for k in r if k.startswith("variance_"))
        assert sats==expected;partial+=len(sats)<len(used)
        decisions[r["decision"]]+=1;valid.append(r)
        if r["motion_checks"]:
            for check in r["motion_checks"].split(";"):
                f=check.split("@");assert len(f)==10 and 0<float(f[1])<=2 and float(f[9])==.1
                assert math.isclose(math.hypot(*(float(v) for v in f[5:8])),float(f[8]),abs_tol=1e-14);checks+=1
        if r["status"]=="FIX":
            assert r["decision"]=="accepted" and int(r["confirmation_count"])>=2 and int(r["ar_ambiguity_pairs"])>=4
            assert float(r["ratio"])>=3 and float(r["velocity_max_std_mps"])<=cap and all(locks[s]>=5 for s in sats)
            assert r["motion_checks"] and float(r["motion_max_residual_m"])<=.1 and float(r["fixed_correction_m"])<=1 and float(r["fixed_phase_rms_m"])<=.01 and float(r["fixed_code_rms_m"])<=1
            assert all(float(r["solution_q"+axis+"_m2"])<=.05**2 for axis in ["xx","yy","zz"])
            assert len(r["integers"].split(";"))==len(sats)-len(refs)
            if all_references:assert int(r["system_reference_sets"])==math.prod(groups.values())
            fixed.append(r)
        else:assert r["decision"]!="accepted" and all(r[k+"_m"]==r["float_"+k+"_m"] for k in ["x","y","z"])
    paired=sum(r["kind"]=="pair" for r in decoded)
    return decoded,valid,fixed,dict(output_rows=len(decoded),paired_epochs=paired,status_counts=dict(collections.Counter(r["status"] for r in decoded)),valid_RTK_positions=len(valid),FIX_positions=len(fixed),FIX_rate_among_pairs=len(fixed)/paired,
        decision_counts=dict(decisions),motion_checks=checks,epochs_with_fewer_AR_than_FLOAT_satellites=partial,all_FLOAT_and_Doppler_columns_identical_to_previous_control=True,all_original_53_columns_identical_to_prior_control=not variance,variance_filtered_epochs=filtered,removed_DD_pairs=removed)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--dataset-root",type=Path,required=True);parser.add_argument("--binary",type=Path,required=True);parser.add_argument("--output-dir",type=Path,required=True);args=parser.parse_args()
    repo=Path(__file__).resolve().parents[1];previous=repo/"docs/comparison/2026-10-21-locked-subset-ar";prior=json.loads((previous/"predeclared-inputs.json").read_text())
    h=runpy.run_path(str(repo/"scripts/run_native_ppc_comparison.py"));sha,save,retain=(h[k] for k in ["sha","save","retain"]);score=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))["score"]
    dataset=args.dataset_root.resolve();binary=args.binary.resolve();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    paths=list((repo/"src").rglob("*.rs"))+[repo/p for p in ["scripts/run_ppc_variance_subset.py","scripts/run_ppc_rtk_float.py","scripts/run_native_ppc_comparison.py","scripts/upstream_rtk_ar_selection.cpp","examples/ppc_rtk_kinematic_fix.rs","tests/rtk_kinematic_fix.rs","tests/fixtures/upstream_rtk_ar_selection.csv","tests/rtk_ar_selection.rs","docs/variance-subset-ar.md","Cargo.toml","Cargo.lock","rust-toolchain.toml"]]
    inputs={}
    for case in prior["cases"]:
        city,iteration=case.split("-")
        for name in ["rover.obs","base.obs","base.nav","reference.csv"]:
            p=dataset/city/iteration/name;expected=next(v for k,v in prior["input_sha256"].items() if k.endswith(f"/{city}/{iteration}/{name}"));assert sha(p)==expected;inputs[str(p)]=expected
    profiles={}
    for base in ["galileo","galileo_slips"]:
        for gates in ["selected","all"]:
            for cap in [1,2]:
                old=base+"_"+gates+"_std"+str(cap)+"_locked"
                for variance in [False,True]:profiles[old+("_variance" if variance else "")]=dict(base_profile=base,control_profile=old,mode=gates+("-slips" if base.endswith("_slips") else "")+("-variance" if variance else "-locked"),motion_std_mps=cap,relative_variance=variance)
    manifest=dict(base_commit=subprocess.check_output(["git","rev-parse","HEAD"],text=True,cwd=repo).strip(),source_sha256={str(p.relative_to(repo)):sha(p) for p in paths},input_sha256=inputs,binary_sha256=sha(binary),cases=prior["cases"],base_ecef_m=prior["base_ecef_m"],profiles=profiles,
        min_lock_epochs=5,min_ambiguity_pairs=4,motion_max_gap_s=2,displacement_residual_m=.1,kinematic_default_motion_std_mps=1,predeclared_motion_std_sensitivity=[1,2],
        subset_selection="one lock-qualified subset with optional native relative-variance prefilter",variance_multiplier=10,variance_minimum_cycles2=1e-4,minimum_filter_pairs=4,progressive_paths_component_only=True,ratio_driven_subset_search=False,integer_hold=False,FLOAT_feedback=False,
        numeric_local_and_motion_thresholds_changed=False,default_static_or_full_set_policy_changed=False,truth_used_by_solver=False,held_out_validation=False,native_rerun=False,practical_PPC_parity_accepted=False,prior_manifest_sha256=sha(previous/"predeclared-inputs.json"))
    save(out/"predeclared-inputs.json",manifest);report=dict(practical_PPC_parity_accepted=False,cases={})
    for case in manifest["cases"]:
        city,iteration=case.split("-");root=dataset/city/iteration;controls={};report["cases"][case]={}
        with (root/"reference.csv").open() as f:reference_rows=[{k.strip():v.strip() for k,v in r.items()} for r in csv.DictReader(f)]
        references={(int(r["GPS Week"]),Decimal(r["GPS TOW (s)"])):r for r in reference_rows}
        for name,c in profiles.items():
            assert all(sha(repo/p)==v for p,v in manifest["source_sha256"].items()) and sha(binary)==manifest["binary_sha256"]
            old=c["control_profile"]
            if old not in controls:controls[old]=rows(gzip.decompress((previous/(case+"."+old+".stdout.csv.gz")).read_bytes()))
            command=[str(binary),str(root),c["mode"],",".join(map(str,manifest["base_ecef_m"][city])),str(c["motion_std_mps"])]
            r=subprocess.run(command,capture_output=True,cwd=repo);retain(out,case+"."+name,command,r.stdout,r.stderr,r.returncode);assert r.returncode==0,r.stderr.decode()
            decoded,valid,fixed,metrics=summarize(r.stdout,controls[old],c["relative_variance"],c["motion_std_mps"],c["mode"].startswith("all"))
            def scores(data):return score([{**r,"status":"FLOAT"} for r in data],references,city)
            metrics["all_RTK_reference_scores"]=scores(valid);metrics["FIX_reference_scores"]=scores(fixed);report["cases"][case][name]=metrics;save(out/"report.json",report)
            print(case,name,metrics["FIX_positions"],metrics["decision_counts"],flush=True)
    assert all(sha(Path(p))==v for p,v in inputs.items())


if __name__=="__main__":main()
