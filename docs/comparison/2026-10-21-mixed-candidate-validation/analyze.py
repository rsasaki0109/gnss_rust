"""Apply only the existing four-pair eligibility bound to recorded diagnostics."""
import argparse
import csv
from decimal import Decimal
import gzip
import json
from pathlib import Path
import runpy


def diagnostics(archive,dataset):
    repo=archive.parents[2]
    scoring=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))
    runner=runpy.run_path(str(repo/"scripts/run_ppc_mixed_candidates.py"))
    manifest=json.loads((archive/"predeclared-inputs.json").read_text())
    result=dict(post_run_only=True,minimum_pairs_from_existing_policy=4,
        lock_confirmation_motion_history_not_evaluated=True,FIX_declared=0,cases={})
    for case in manifest["cases"]:
        city,iteration=case.split("-");references={};result["cases"][case]={}
        with (dataset/city/iteration/"reference.csv").open() as stream:
            for raw in csv.DictReader(stream):
                r={k.strip():v.strip() for k,v in raw.items()};references[(int(r["GPS Week"]),Decimal(r["GPS TOW (s)"]))]=r
        for profile in manifest["diagnostic_modes"]:
            with gzip.open(archive/(case+"."+profile+".diagnostic.stdout.csv.gz"),"rt") as stream:rows=list(csv.DictReader(stream))
            selected=runner["selections"](rows)["ratio_and_all_reference_gates"]
            eligible=[r for r in selected if int(r["ambiguity_pairs"])>=4]
            scored=[{**r,"x_m":r["conditioned_x_m"],"y_m":r["conditioned_y_m"],"z_m":r["conditioned_z_m"]} for r in eligible]
            result["cases"][case][profile]=dict(
                ratio_and_all_reference_gates=len(selected),below_existing_minimum_pairs=len(selected)-len(eligible),
                additionally_meets_existing_minimum_pairs=len(eligible),candidate_reference_scores=scoring["score"](scored,references,city))
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root",type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(diagnostics(Path(__file__).resolve().parent,args.dataset_root),indent=2))
