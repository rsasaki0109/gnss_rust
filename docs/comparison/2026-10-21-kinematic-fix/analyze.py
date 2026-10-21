"""Read-only status/continuity comparison with the previously frozen native run."""
import collections
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
from pathlib import Path
import runpy


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def rows(path):return list(csv.DictReader(io.StringIO(gzip.decompress(path.read_bytes()).decode())))
def key(row):return int(row["week"]),Decimal(row["tow"])
def time(row):return Decimal(row["week"])*604800+Decimal(row["tow"])


def continuity(data):
    previous=None;length=0;runs=[]
    for r in data:
        if r["kind"]!="pair":continue
        if r["status"]=="FIX":
            if previous is not None and time(r)-previous==1:length+=1
            else:
                if length:runs.append(length)
                length=1
            previous=time(r)
        else:
            if length:runs.append(length)
            length=0;previous=None
    if length:runs.append(length)
    return dict(FIX_epochs=sum(runs),one_Hz_FIX_runs=runs,longest_FIX_run_epochs=max(runs,default=0))


def analyze(archive,dataset=None):
    repo=archive.parents[2];native=repo/"docs/comparison/2026-10-21-gps-galileo-float/native-comparison"
    manifest=json.loads((archive/"predeclared-inputs.json").read_text());report=json.loads((archive/"report.json").read_text())
    score=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))["score"]
    result=dict(native_rerun=False,native_commit="72f3b7c2c4c088dc575286a6f0edf4e407bc499a",native_source="2026-10-21-gps-galileo-float/native-comparison/*.spp.stdout.csv.gz",practical_PPC_parity_accepted=False,
        native_evidence_sha256={},profiles={},cases={})
    for name in ["predeclared-inputs.json","report.json","build.json","helper-checks.json"]:
        path=native/name;result["native_evidence_sha256"][str(path.relative_to(repo))]=sha(path)
    for case in manifest["cases"]:
        path=native/(case+".spp.stdout.csv.gz");result["native_evidence_sha256"][str(path.relative_to(repo))]=sha(path)
        nr=rows(path);native_pairs={key(r):r for r in nr if r["kind"]=="pair"};native_valid={k:r for k,r in native_pairs.items() if r["status"] in {"FLOAT","FIX"}}
        references=None;city,iteration=case.split("-")
        if dataset is not None:
            with (dataset/city/iteration/"reference.csv").open() as f:
                reference_rows=[{k.strip():v.strip() for k,v in r.items()} for r in csv.DictReader(f)]
            references={(int(r["GPS Week"]),Decimal(r["GPS TOW (s)"])):r for r in reference_rows}
        cr=dict(native_status_counts=dict(collections.Counter(r["status"] for r in native_pairs.values())),native_continuity=continuity(nr),profiles={});result["cases"][case]=cr
        for profile in manifest["profiles"]:
            rr=rows(archive/(case+"."+profile+".stdout.csv.gz"));rp={key(r):r for r in rr if r["kind"]=="pair"};assert rp.keys()==native_pairs.keys()
            rust_valid={k:r for k,r in rp.items() if r["status"] in {"FLOAT","FIX"}};common=rust_valid.keys()&native_valid.keys()
            rust_fixed={k for k,r in rp.items() if r["status"]=="FIX"};native_fixed={k for k,r in native_pairs.items() if r["status"]=="FIX"}
            metrics=dict(status_transition_counts=dict(collections.Counter(native_pairs[k]["status"]+"->"+r["status"] for k,r in rp.items())),
                common_RTK_epochs=len(common),Rust_only_RTK_epochs=len(rust_valid.keys()-native_valid.keys()),native_only_RTK_epochs=len(native_valid.keys()-rust_valid.keys()),
                common_FIX_epochs=len(rust_fixed&native_fixed),Rust_only_FIX_epochs=len(rust_fixed-native_fixed),native_only_FIX_epochs=len(native_fixed-rust_fixed),
                Rust_continuity=continuity(rr))
            if references is not None:
                def scores(data):return score([{**r,"status":"FLOAT"} for r in data],references,city)
                metrics["common_RTK_reference_scores"]={"Rust":scores([rust_valid[k] for k in sorted(common)]),"native":scores([native_valid[k] for k in sorted(common)])}
            cr["profiles"][profile]=metrics
            totals=result["profiles"].setdefault(profile,dict(paired_epochs=0,FLOAT_epochs=0,FIX_epochs=0,INVALID_epochs=0,motion_checks=0,decision_counts={},native_RTK_epochs=0,native_FIX_epochs=0,common_FIX_epochs=0))
            run=report["cases"][case][profile];totals["paired_epochs"]+=run["paired_epochs"]
            for s in ["FLOAT","FIX","INVALID"]:totals[s+"_epochs"]+=run["status_counts"].get(s,0)
            totals["motion_checks"]+=run["motion_checks"];totals["decision_counts"]=dict(collections.Counter(totals["decision_counts"])+collections.Counter(run["decision_counts"]))
            totals["native_RTK_epochs"]+=len(native_valid);totals["native_FIX_epochs"]+=len(native_fixed);totals["common_FIX_epochs"]+=metrics["common_FIX_epochs"]
    return result


if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--dataset-root",type=Path);args=parser.parse_args()
    archive=Path(__file__).resolve().parent
    print(json.dumps(analyze(archive,args.dataset_root),indent=2))
