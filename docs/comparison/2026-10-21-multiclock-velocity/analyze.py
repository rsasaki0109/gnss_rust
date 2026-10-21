"""Post-run observability diagnostics; never choose thresholds using reference errors."""
import argparse
import csv
from decimal import Decimal
import gzip
import json
import math
from pathlib import Path
import runpy


def diagnostics(archive,dataset):
    repo=archive.parents[2];rover=runpy.run_path(str(repo/"scripts/run_ppc_rover_velocity.py"));targets=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))["reference_targets"]
    manifest=json.loads((archive/"predeclared-inputs.json").read_text());result=dict(post_run_only=True,thresholds_selected_using_reference=False,cases={})
    for case in manifest["cases"]:
        city,iteration=case.split("-");references=rover["reference_velocities"](dataset/city/iteration/"reference.csv",city,targets);result["cases"][case]={}
        for profile in manifest["profiles"]:
            with gzip.open(archive/(case+"."+profile+".mixed.stdout.csv.gz"),"rt") as stream:rows=list(csv.DictReader(stream))
            valid=[r for r in rows if r["velocity_status"]=="VELOCITY"];by_dof={};worst=[]
            dof=lambda r:int(r["velocity_prepared_count"])-3-len(r["velocity_clock_systems"].split(";"))
            for freedom in sorted({dof(r) for r in valid}):
                group=[r for r in valid if dof(r)==freedom]
                by_dof[str(freedom)]=dict(epochs=len(group),reference_scores=rover["scores"](group,references))
            for row in valid:
                key=(int(row["week"]),Decimal(row["tow"]))
                if key in references:
                    error=math.dist([float(row[k]) for k in ["vx_mps","vy_mps","vz_mps"]],references[key][0]["reference_is_antenna"])
                    worst.append(dict(error_3d_mps=error,row=row))
            result["cases"][case][profile]=dict(by_degrees_of_freedom=by_dof,worst_reference_is_antenna=sorted(worst,key=lambda r:r["error_3d_mps"],reverse=True)[:5])
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--dataset-root",type=Path,required=True);args=parser.parse_args()
    print(json.dumps(diagnostics(Path(__file__).resolve().parent,args.dataset_root),indent=2))
