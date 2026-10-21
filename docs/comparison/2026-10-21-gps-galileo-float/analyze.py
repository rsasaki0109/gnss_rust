"""Post-run shared-epoch diagnostics, without changing any solver decisions."""
import argparse
import csv
from decimal import Decimal
import gzip
import json
import math
from pathlib import Path
import runpy


def load(path):
    with gzip.open(path,"rt") as stream:return list(csv.DictReader(stream))


def diagnostics(archive,dataset):
    scoring=runpy.run_path(str(archive.parents[2]/"scripts/run_ppc_rtk_float.py"))
    report=json.loads((archive/"report.json").read_text())
    result=dict(post_run_only=True,reference_hypotheses_unchanged=True,cases={})
    key=lambda r:(int(r["week"]),float(r["tow"]).hex())
    for case in report["cases"]:
        city,run=case.split("-");references={}
        with (dataset/city/run/"reference.csv").open() as stream:
            for raw in csv.DictReader(stream):
                r={k.strip():v.strip() for k,v in raw.items()};references[(int(r["GPS Week"]),Decimal(r["GPS TOW (s)"]))]=r
        native_rows=load(archive/"native-comparison"/(case+".spp.stdout.csv.gz"))
        native={key(r):r for r in native_rows if r["status"] in {"FLOAT","FIX"}}
        pairs=sum(r["kind"]=="pair" for r in native_rows)
        result["cases"][case]={}
        for profile in ["galileo","galileo_slips"]:
            rust_rows=load(archive/(case+"."+profile+".stdout.csv.gz"))
            rust={key(r):r for r in rust_rows if r["status"]=="FLOAT"}
            common=native.keys()&rust.keys()
            delta=[math.dist([float(native[t][k]) for k in ["x_m","y_m","z_m"]],
                             [float(rust[t][k]) for k in ["x_m","y_m","z_m"]]) for t in sorted(common)]
            score=lambda rows:scoring["score"]([{**r,"status":"FLOAT"} for r in rows],references,city)
            result["cases"][case][profile]=dict(
                native_default_SPP_RTK_epochs=len(native),rust_FLOAT_epochs=len(rust),
                common_RTK_epochs=len(common),native_RTK_only=len(native.keys()-rust.keys()),
                rust_FLOAT_only=len(rust.keys()-native.keys()),neither_RTK=pairs-len(native.keys()|rust.keys()),
                common_epochs_reference_scores=dict(native=score(native[t] for t in sorted(common)),
                                                   rust=score(rust[t] for t in sorted(common))),
                common_solver_position_difference_m=dict(epochs=len(delta),rms=math.sqrt(sum(d*d for d in delta)/len(delta)),maximum=max(delta)) if delta else None,
            )
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root",type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(diagnostics(Path(__file__).resolve().parent,args.dataset_root),indent=2))
