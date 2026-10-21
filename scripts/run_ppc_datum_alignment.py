"""Freeze all seven PPC controls and their per-system SD datum alignment variants."""
import argparse
import collections
import csv
from decimal import Decimal
import gzip
import json
from pathlib import Path
import runpy
import subprocess


def alignment_summary(rows):
    counts=collections.Counter();systems=collections.Counter();maximum=0.0
    for r in rows:
        if r["status"]!="FLOAT":continue
        retained=r["datum_retained_satellites"].split(";") if r["datum_retained_satellites"] else []
        reset=r["datum_initialized_satellites"].split(";") if r["datum_initialized_satellites"] else []
        used=r["used_satellites"].split(";")
        assert set(retained).isdisjoint(reset) and set(retained+reset)==set(used)
        assert len(reset)==int(r["ambiguity_resets"])
        offsets={s.split(":")[0]:float(s.split(":")[1]) for s in r["datum_offsets_cycles"].split(";")}
        assert set(offsets)=={s[0] for s in used}
        for system,offset in offsets.items():
            old=[s for s in retained if s[0]==system];new=[s for s in reset if s[0]==system]
            if not old:assert offset==0.0
            counts["system_rows"]+=1
            counts["retained_and_new_system_rows"]+=bool(old and new)
            counts["retained_satellite_rows"]+=len(old)
            counts["initialized_satellite_rows"]+=len(new)
            systems[system]+=1;maximum=max(maximum,abs(offset))
    return dict(counts=dict(counts),systems=dict(systems),maximum_abs_applied_offset_cycles=maximum)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root",type=Path,required=True)
    parser.add_argument("--binary",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args();repo=Path(__file__).resolve().parents[1]
    previous=repo/"docs/comparison/2026-10-21-gps-galileo-float"
    prior=json.loads((previous/"predeclared-inputs.json").read_text())
    helpers=runpy.run_path(str(repo/"scripts/run_native_ppc_comparison.py"))
    sha,save,retain=(helpers[k] for k in ["sha","save","retain"])
    scoring=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))
    dataset=args.dataset_root.resolve();binary=args.binary.resolve();output=args.output_dir.resolve()
    output.mkdir(parents=True,exist_ok=False)
    sources=list((repo/"src").rglob("*.rs"))+[repo/p for p in [
        "scripts/run_ppc_datum_alignment.py","scripts/run_ppc_gauge_probe.py","scripts/run_ppc_rtk_float.py","scripts/run_native_ppc_comparison.py",
        "examples/ppc_rtk_float.rs","examples/ppc_rtk_gauge.rs","tests/rtk_ambiguity.rs","tests/rtk_multisystem.rs",
        "scripts/upstream_rtk_datum.cpp","tests/fixtures/upstream_rtk_datum.csv","Cargo.toml","Cargo.lock","rust-toolchain.toml"]]
    inputs={}
    for case in prior["cases"]:
        city,run=case.split("-")
        for name in ["rover.obs","base.obs","base.nav","reference.csv"]:
            p=dataset/city/run/name
            expected=next(h for k,h in prior["input_sha256"].items() if k.endswith(f"/{city}/{run}/{name}"))
            assert sha(p)==expected;inputs[str(p)]=expected
    profiles={}
    for name,mode in prior["profiles"].items():profiles[name]=mode;profiles[name+"_aligned"]=mode+"-align"
    manifest=dict(base_commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=repo,text=True).strip(),
        source_sha256={str(p.relative_to(repo)):sha(p) for p in sources},input_sha256=inputs,binary_sha256=sha(binary),
        prior_manifest_sha256=sha(previous/"predeclared-inputs.json"),prior_report_sha256=sha(previous/"report.json"),
        cases=prior["cases"],profiles=profiles,base_ecef_m=prior["base_ecef_m"],
        changed_setting="align_ambiguity_datums false -> true only; means per system, covariance untouched",
        thresholds_changed=False,reference_used_for_configuration=False,default_changed=False,
        native_whole_solver_rerun=False,practical_PPC_parity_accepted=False)
    save(output/"predeclared-inputs.json",manifest);report=dict(practical_PPC_parity_accepted=False,cases={})
    for case in manifest["cases"]:
        city,run=case.split("-");references={}
        with (dataset/city/run/"reference.csv").open() as stream:
            for raw in csv.DictReader(stream):
                r={k.strip():v.strip() for k,v in raw.items()};references[(int(r["GPS Week"]),Decimal(r["GPS TOW (s)"]))]=r
        report["cases"][case]={};ledger=None
        for profile,mode in profiles.items():
            assert all(sha(repo/p)==h for p,h in manifest["source_sha256"].items())
            assert sha(binary)==manifest["binary_sha256"]
            command=[str(binary),str(dataset/city/run),mode,",".join(map(str,manifest["base_ecef_m"][city]))]
            result=subprocess.run(command,capture_output=True,cwd=repo)
            retain(output,case+"."+profile,command,result.stdout,result.stderr,result.returncode)
            assert result.returncode==0,result.stderr.decode()
            rows,metrics=scoring["inspect"](result.stdout)
            keys=[(r["kind"],int(r["week"]),float(r["tow"]).hex(),r["side"]) for r in rows]
            assert ledger is None or ledger==keys;ledger=keys
            if profile.endswith("_aligned"):metrics["datum_alignment"]=alignment_summary(rows)
            else:
                with gzip.open(previous/(case+"."+profile+".stdout.csv.gz"),"rt") as stream:assert rows==list(csv.DictReader(stream))
                metrics["all_columns_identical_to_prior_control"]=True
            metrics["reference_scores"]=scoring["score"](rows,references,city)
            report["cases"][case][profile]=metrics;save(output/"report.json",report)
            print(case,profile,metrics["status_counts"].get("FLOAT",0),flush=True)
    assert all(sha(Path(p))==h for p,h in inputs.items())


if __name__=="__main__":main()
