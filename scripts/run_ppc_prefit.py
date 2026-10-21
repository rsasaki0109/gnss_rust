"""Evaluate frozen 3 m whole-epoch prefit gates with RTK/SPP seed controls."""
import argparse
import collections
import csv
from decimal import Decimal
import gzip
import json
from pathlib import Path
import runpy
import subprocess


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def prefit_summary(rows):
    counts=collections.Counter()
    sources=collections.Counter()
    failures=collections.Counter()
    for r in rows:
        if not r["prefit_threshold_m"]:
            assert r["status"]!="FLOAT" and r["reason"]!="PrefitResidualRejected"
            continue
        assert r["kind"]=="pair" and float(r["prefit_threshold_m"])==3.0
        phase,code=map(float,[r["prefit_phase_max_m"],r["prefit_code_max_m"]])
        rejected=sum(int(r[k]) for k in ["prefit_rejected_phase_rows","prefit_rejected_code_rows"])
        pairs=r["prefit_rejected_pairs"].split(";") if r["prefit_rejected_pairs"] else []
        assert len(pairs)==rejected
        assert (rejected>0)==(phase>3.0 or code>3.0)
        assert (rejected>0)==(r["reason"]=="PrefitResidualRejected")
        counts["screened_epochs"]+=1
        counts["prefit_rejected_epochs"]+=rejected>0
        counts["passed_screen_then_failed_update"]+=rejected==0 and r["status"]!="FLOAT"
        counts["screened_accepted_FLOAT"]+=r["status"]=="FLOAT"
        counts["rejected_phase_rows"]+=int(r["prefit_rejected_phase_rows"])
        counts["rejected_code_rows"]+=int(r["prefit_rejected_code_rows"])
        for k in ["code_slips","doppler_slips","code_tests","doppler_tests"]:
            counts[k+"_among_screened_attempts"]+=int(r[k])
        sources[r["seed_source"]]+=1
        if r["spp_failure"]:failures[r["spp_failure"]]+=1
    return dict(counts=dict(counts),seed_sources=dict(sources),spp_failures=dict(failures))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root",type=Path,required=True)
    parser.add_argument("--binary",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    repo=Path(__file__).resolve().parents[1]
    baseline=repo/"docs/comparison/2026-10-21-native-rtk-common-inputs"
    spp_baseline=repo/"docs/comparison/2026-10-21-rust-spp-seed"
    control=json.loads((baseline/"predeclared-inputs.json").read_text())
    scoring=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))
    capture=runpy.run_path(str(repo/"scripts/run_native_ppc_comparison.py"))["retain"]
    sha=scoring["sha"]
    output,dataset,binary=args.output_dir.resolve(),args.dataset_root.resolve(),args.binary.resolve()
    output.mkdir(parents=True,exist_ok=False)
    sources=list((repo/"src").rglob("*.rs"))+[repo/n for n in [
        "scripts/run_ppc_prefit.py","scripts/run_ppc_rtk_float.py","scripts/run_native_ppc_comparison.py",
        "examples/ppc_rtk_float.rs","tests/rtk_float.rs","tests/rtk_quality.rs",
        "scripts/upstream_rtk_prefit.cpp","tests/fixtures/upstream_rtk_prefit.csv",
        "Cargo.toml","Cargo.lock","rust-toolchain.toml",
    ]]
    inputs={}
    for case in control["cases"]:
        city,run=case.split("-")
        for name in ["rover.obs","base.obs","base.nav","reference.csv"]:
            p=dataset/city/run/name
            k=next(k for k in control["input_sha256"] if k.endswith(f"/{city}/{run}/{name}"))
            assert sha(p)==control["input_sha256"][k]
            inputs[str(p)]=sha(p)
    manifest=dict(
        source_sha256={str(p.relative_to(repo)):sha(p) for p in sources},input_sha256=inputs,
        binary_sha256=sha(binary),cases=control["cases"],base_ecef_m=control["base_ecef_m"],
        profiles={"none":"kinematic","spp":"kinematic-spp","gate":"kinematic-prefit3",
                  "spp_gate":"kinematic-spp-prefit3","slips_gate":"kinematic-slips-prefit3",
                  "spp_slips_gate":"kinematic-spp-slips-prefit3"},
        max_prefit_residual_m=3.0,threshold_source="pinned native RTKConfig::outlier_threshold",
        action="reject whole epoch; no native H/v zeroing and no covariance mutation",
        classification="once at initial seed, after ambiguity resets, before iterated update",
        native_whole_solver_rerun=False, native_whole_policy_equivalence_claimed=False,
        other_thresholds_changed=False,reference_used_to_choose_configuration=False,
        kinematic_FIX_enabled=False,timing_memory_acceptance=False,practical_PPC_parity_accepted=False,
    )
    save(output/"predeclared-inputs.json",manifest)
    report=dict(practical_PPC_parity_accepted=False,cases={})
    for case in control["cases"]:
        city,run=case.split("-")
        root=dataset/city/run
        references={}
        with (root/"reference.csv").open() as stream:
            for raw in csv.DictReader(stream):
                r={k.strip():v.strip() for k,v in raw.items()}
                references[(int(r["GPS Week"]),Decimal(r["GPS TOW (s)"]))]=r
        report["cases"][case]={}
        ledger=None
        for profile,mode in manifest["profiles"].items():
            assert all(sha(repo/p)==h for p,h in manifest["source_sha256"].items())
            assert sha(binary)==manifest["binary_sha256"]
            command=[str(binary),str(root),mode,",".join(map(str,control["base_ecef_m"][city]))]
            result=subprocess.run(command,capture_output=True,cwd=repo)
            capture(output,case+"."+profile,command,result.stdout,result.stderr,result.returncode)
            assert result.returncode==0,result.stderr.decode()
            rows,metrics=scoring["inspect"](result.stdout)
            keys=[(r["kind"],int(r["week"]),float(r["tow"]).hex(),r["side"]) for r in rows]
            assert ledger is None or ledger==keys
            ledger=keys
            if profile in {"none","spp"}:
                path=(baseline/(case+".rust.stdout.csv.gz") if profile=="none" else spp_baseline/(case+".spp.stdout.csv.gz"))
                with gzip.open(path,"rt") as stream:assert rows==list(csv.DictReader(stream))
                metrics["all_columns_identical_to_prior_control"]=True
            else:metrics["prefit"]=prefit_summary(rows)
            metrics["reference_scores"]=scoring["score"](rows,references,city)
            report["cases"][case][profile]=metrics
            save(output/"report.json",report)
            print(case,profile,metrics["status_counts"].get("FLOAT",0),metrics.get("prefit",{}).get("counts",{}),flush=True)
    assert all(sha(Path(p))==h for p,h in inputs.items())
    assert all(sha(repo/p)==h for p,h in manifest["source_sha256"].items())


if __name__=="__main__":main()
