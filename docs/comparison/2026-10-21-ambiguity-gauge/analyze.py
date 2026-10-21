"""Post-run coverage and error diagnostics; no solver configuration changes."""
import argparse
import csv
from decimal import Decimal
import gzip
import json
import math
from pathlib import Path
import runpy


def rows(path):
    with gzip.open(path,"rt") as stream:return list(csv.DictReader(stream))


def compare(a,b,pairs,references,city,score):
    common=sorted(a.keys()&b.keys())
    scores=lambda data:score([{**data[t],"status":"FLOAT"} for t in common],references,city)
    differences=[math.dist([float(a[t][k]) for k in ["x_m","y_m","z_m"]],
                           [float(b[t][k]) for k in ["x_m","y_m","z_m"]]) for t in common]
    return dict(a_epochs=len(a),b_epochs=len(b),common_epochs=len(common),
        a_only=len(a.keys()-b.keys()),b_only=len(b.keys()-a.keys()),neither=pairs-len(a.keys()|b.keys()),
        common_reference_scores=dict(a=scores(a),b=scores(b)),
        common_position_difference_m=dict(epochs=len(differences),maximum=max(differences),
            rms=math.sqrt(sum(d*d for d in differences)/len(differences))) if differences else None)


def diagnostics(archive,dataset):
    repo=archive.parents[2];run=archive/"ppc-comparison"
    manifest=json.loads((run/"predeclared-inputs.json").read_text())
    score=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))["score"]
    previous=repo/"docs/comparison/2026-10-21-gps-galileo-float/native-comparison"
    key=lambda r:(int(r["week"]),float(r["tow"]).hex())
    result=dict(post_run_only=True,configuration_selected_using_reference=False,
        native_outputs_reused_from_prior_immutable_archive=True,cases={})
    for case in manifest["cases"]:
        city,iteration=case.split("-");references={}
        with (dataset/city/iteration/"reference.csv").open() as stream:
            for raw in csv.DictReader(stream):
                r={k.strip():v.strip() for k,v in raw.items()};references[(int(r["GPS Week"]),Decimal(r["GPS TOW (s)"]))]=r
        result["cases"][case]={};solutions={};pairs=0
        for profile in manifest["profiles"]:
            data=rows(run/(case+"."+profile+".stdout.csv.gz"))
            solutions[profile]={key(r):r for r in data if r["status"]=="FLOAT"}
            pairs=sum(r["kind"]=="pair" for r in data)
        for profile in manifest["profiles"]:
            if profile.endswith("_aligned"):continue
            result["cases"][case][profile]=compare(solutions[profile],solutions[profile+"_aligned"],pairs,references,city,score)
        native={key(r):r for r in rows(previous/(case+".spp.stdout.csv.gz")) if r["status"] in {"FLOAT","FIX"}}
        result["cases"][case]["native_default_SPP_vs_aligned_Galileo"]={}
        for profile in ["galileo_aligned","galileo_slips_aligned"]:
            result["cases"][case]["native_default_SPP_vs_aligned_Galileo"][profile]=compare(native,solutions[profile],pairs,references,city,score)
        special={}
        if case=="nagoya-run2":
            for profile,data in solutions.items():
                special[profile]={}
                for tow in [556685.0,556760.0]:
                    t=(2323,tow.hex());r=data.get(t)
                    special[profile][str(tow)]={"status":"FLOAT","row":r,"reference_score":score([r],references,city)} if r else {"status":"NO_FLOAT"}
        result["cases"][case]["prior_GPS_outlier_times"]=special
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root",type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(diagnostics(Path(__file__).resolve().parent,args.dataset_root),indent=2))
