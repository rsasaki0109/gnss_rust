"""Post-run velocity diagnostics grouped by observed satellite count; no gate changes."""
import argparse
import csv
from decimal import Decimal
import gzip
import json
import math
from pathlib import Path
import runpy


def diagnostics(archive,dataset):
    repo=archive.parents[2];runner=runpy.run_path(str(repo/"scripts/run_ppc_rover_velocity.py"))
    targets=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))["reference_targets"]
    manifest=json.loads((archive/"predeclared-inputs.json").read_text())
    result=dict(post_run_only=True,configuration_selected_using_reference=False,cases={})
    for case in manifest["cases"]:
        city,iteration=case.split("-");references=runner["reference_velocities"](dataset/city/iteration/"reference.csv",city,targets)
        result["cases"][case]={}
        for profile in manifest["profiles"]:
            with gzip.open(archive/(case+"."+profile+".stdout.csv.gz"),"rt") as stream:rows=list(csv.DictReader(stream))
            valid=[r for r in rows if r["velocity_status"]=="VELOCITY"];groups={};worst=[]
            for count in sorted({int(r["velocity_prepared_count"]) for r in valid}):
                selected=[r for r in valid if int(r["velocity_prepared_count"])==count]
                groups[str(count)]=dict(epochs=len(selected),reference_scores=runner["scores"](selected,references))
            for r in valid:
                key=(int(r["week"]),Decimal(r["tow"]))
                if key not in references:continue
                velocity=references[key][0]["reference_is_antenna"]
                error=math.dist([float(r[k]) for k in ["vx_mps","vy_mps","vz_mps"]],velocity)
                worst.append(dict(error_3d_mps=error,row=r))
            candidates=[r for r in rows if r["conservative_rejection"]=="PASS_LOCAL_GATES" and int(r["ambiguity_pairs"])>=4]
            result["cases"][case][profile]=dict(by_GPS_Doppler_count=groups,worst_reference_is_antenna=sorted(worst,key=lambda x:x["error_3d_mps"],reverse=True)[:5],
                ratio_all_reference_four_DD_snapshots=len(candidates),with_velocity=sum(r["velocity_status"]=="VELOCITY" for r in candidates))
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--dataset-root",type=Path,required=True)
    args=parser.parse_args();print(json.dumps(diagnostics(Path(__file__).resolve().parent,args.dataset_root),indent=2))
