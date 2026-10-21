"""Read-only worst-point diagnostics; truth never feeds positioning or gates."""
import csv
from decimal import Decimal
import gzip
import json
import math
from pathlib import Path
import runpy


def diagnose(archive,dataset):
    repo=archive.parents[2];targets=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))["reference_targets"]
    m=json.loads((archive/"predeclared-inputs.json").read_text());result=dict(scope="five worst antenna-assumption QC-gated FLOAT points per case; selected after solving, diagnostics only",cases={})
    def key(r):return int(r["week"]),Decimal(r["tow"])
    for case in m["cases"]:
        city,iteration=case.split("-")
        with (dataset/city/iteration/"reference.csv").open() as f:refs={(int(r["GPS Week"].strip()),Decimal(r["GPS TOW (s)"].strip())):{k.strip():v.strip() for k,v in r.items()} for r in ({k.strip():v for k,v in row.items()} for row in csv.DictReader(f))}
        name="galileo_spp_slips_gate_ggspp_nativevar"
        def load(profile):
            with gzip.open(archive/(case+"."+profile+".stdout.csv.gz"),"rt") as f:return {key(r):r for r in csv.DictReader(f) if r["kind"]=="pair"}
        old=load(name);new=load(name+"_sppqc")
        def errors(r,prefix=""):
            if r["status"]!="FLOAT" or key(r) not in refs:return None
            positions,_,_=targets(refs[key(r)],city);xyz=[float(r[prefix+k]) for k in ["x_m","y_m","z_m"]]
            return {n:math.sqrt(sum((a-b)**2 for a,b in zip(xyz,p))) for n,p in positions.items()}
        chosen=sorted((r for r in new.values() if r["status"]=="FLOAT"),key=lambda r:errors(r)["reference_is_antenna"],reverse=True)[:5]
        result["cases"][case]=[dict(QC_row=r,position_errors_m=errors(r),seed_position_errors_m=errors(r,"seed_") if r["seed_source"]=="SinglePoint" else None,without_QC_row=old[key(r)],without_QC_position_errors_m=errors(old[key(r)])) for r in chosen]
    return result


if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--dataset-root",type=Path,required=True);args=p.parse_args();print(json.dumps(diagnose(Path(__file__).resolve().parent,args.dataset_root),indent=2))
