"""Read-only common-epoch comparison; native captures are reused unchanged."""
import collections
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
from pathlib import Path
import runpy


def rows(path):return list(csv.DictReader(io.StringIO(gzip.decompress(path.read_bytes()).decode())))
def key(r):return int(r["week"]),Decimal(r["tow"])
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def analyze(archive,dataset=None):
    repo=archive.parents[2];m=json.loads((archive/"predeclared-inputs.json").read_text());report=json.loads((archive/"report.json").read_text());native=repo/"docs/comparison/2026-10-21-gps-galileo-float/native-comparison"
    score=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))["score"]
    result=dict(native_whole_rerun=False,Rust_FIX_policy_evaluated=False,practical_PPC_parity_accepted=False,native_evidence_sha256={},profiles={},cases={})
    for name in ["predeclared-inputs.json","report.json","build.json"]:p=native/name;result["native_evidence_sha256"][str(p.relative_to(repo))]=sha(p)
    for case in m["cases"]:
        city,iteration=case.split("-");references=None
        if dataset:
            with (dataset/city/iteration/"reference.csv").open() as f:refs=[{k.strip():v.strip() for k,v in r.items()} for r in csv.DictReader(f)]
            references={(int(r["GPS Week"]),Decimal(r["GPS TOW (s)"])):r for r in refs}
        def scores(data):return score([{**r,"status":"FLOAT"} for r in data],references,city)
        path=native/(case+".spp.stdout.csv.gz");result["native_evidence_sha256"][str(path.relative_to(repo))]=sha(path);nr={key(r):r for r in rows(path) if r["kind"]=="pair"};nv={k:r for k,r in nr.items() if r["status"] in {"FLOAT","FIX"}}
        cr=dict(native_status_counts=dict(collections.Counter(r["status"] for r in nr.values())),profiles={},GPS_vs_mixed_seed={});result["cases"][case]=cr;captured={}
        for profile in m["profiles"]:
            rr=rows(archive/(case+"."+profile+".stdout.csv.gz"));rp={key(r):r for r in rr if r["kind"]=="pair"};assert rp.keys()==nr.keys();rv={k:r for k,r in rp.items() if r["status"]=="FLOAT"};common=sorted(rv.keys()&nv.keys());captured[profile]=rp
            metrics=dict(status_transitions_native_to_Rust=dict(collections.Counter(nr[k]["status"]+"->"+r["status"] for k,r in rp.items())),common_RTK_epochs=len(common),Rust_only_RTK_epochs=len(rv.keys()-nv.keys()),native_only_RTK_epochs=len(nv.keys()-rv.keys()))
            if references is not None:metrics["common_RTK_reference_scores"]={"Rust":scores([rv[k] for k in common]),"native":scores([nv[k] for k in common])}
            cr["profiles"][profile]=metrics
            t=result["profiles"].setdefault(profile,dict(paired_epochs=0,FLOAT_epochs=0,INVALID_epochs=0,accepted_SPP_seed_epochs=0,accepted_seed_fallbacks=0))
            run=report["cases"][case][profile];t["paired_epochs"]+=run["paired_epochs"];t["FLOAT_epochs"]+=len(rv);t["INVALID_epochs"]+=run["status_counts"].get("INVALID",0)
            s=run["SPP_seed"]["accepted_seed_sources"].get("SinglePoint",0);t["accepted_SPP_seed_epochs"]+=s;t["accepted_seed_fallbacks"]+=len(rv)-s
        for profile,c in m["profiles"].items():
            if not c["mixed_seed"]:continue
            old=c["control_profile"];a=captured[old];b=captured[profile];av={k:r for k,r in a.items() if r["status"]=="FLOAT"};bv={k:r for k,r in b.items() if r["status"]=="FLOAT"};common=sorted(av.keys()&bv.keys());seeds=[k for k in common if av[k]["seed_source"]==bv[k]["seed_source"]=="SinglePoint"]
            metrics=dict(pair_status_transitions=dict(collections.Counter(a[k]["status"]+"->"+b[k]["status"] for k in a)),common_FLOAT_epochs=len(common),GPS_seed_only_FLOAT_epochs=len(av.keys()-bv.keys()),mixed_seed_only_FLOAT_epochs=len(bv.keys()-av.keys()),common_successful_seed_epochs=len(seeds))
            if references is not None:
                metrics["common_FLOAT_reference_scores"]={"GPS_seed":scores([av[k] for k in common]),"mixed_seed":scores([bv[k] for k in common])}
                metrics["common_successful_seed_reference_scores"]={}
                for name,data in [("GPS_seed",av),("mixed_seed",bv)]:metrics["common_successful_seed_reference_scores"][name]=scores([{**data[k],"x_m":data[k]["seed_x_m"],"y_m":data[k]["seed_y_m"],"z_m":data[k]["seed_z_m"]} for k in seeds])
            cr["GPS_vs_mixed_seed"][profile]=metrics
    return result


if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--dataset-root",type=Path);args=p.parse_args();print(json.dumps(analyze(Path(__file__).resolve().parent,args.dataset_root),indent=2))
