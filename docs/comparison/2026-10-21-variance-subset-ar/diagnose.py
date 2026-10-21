"""Read-only preselection observability and changes, independent of reference truth."""
import collections
import csv
import gzip
import io
import json
from pathlib import Path


def rows(path):return list(csv.DictReader(io.StringIO(gzip.decompress(path.read_bytes()).decode())))


def diagnose(archive):
    m=json.loads((archive/"predeclared-inputs.json").read_text());result=dict(reference_truth_read=False,cases={})
    for case in m["cases"]:
        result["cases"][case]={}
        for profile,c in m["profiles"].items():
            if not c["relative_variance"]:continue
            enabled=rows(archive/(case+"."+profile+".stdout.csv.gz"));disabled=rows(archive/(case+"."+c["control_profile"]+".stdout.csv.gz"))
            changed=collections.Counter();decisions=collections.Counter();applied=[];admitted=0
            for r,b in zip(enabled,disabled):
                assert (r["kind"],r["week"],r["tow"],r["side"])==(b["kind"],b["week"],b["tow"],b["side"])
                for key in list(b)[:53]:
                    if r[key]!=b[key]:changed[key]+=1
                if r["variance_applied"]=="true":
                    values=[float(v) for v in r["variance_dd_cycles2"].split(";")];median=float(r["variance_median_cycles2"])
                    motion=r["velocity_status"]=="VELOCITY" and float(r["velocity_max_std_mps"])<=c["motion_std_mps"];admitted+=motion;decisions[r["decision"]]+=1
                    applied.append(dict(week=int(r["week"]),tow=r["tow"],decision=r["decision"],prior_decision=b["decision"],motion_admitted=motion,maximum_to_upper_median=max(values)/median if median else None,
                        prior_AR_satellites=b["ar_satellites"],AR_satellites=r["ar_satellites"],prior_ratio=b["ratio"],ratio=r["ratio"]))
            result["cases"][case][profile]=dict(changed_original_53_column_counts=dict(changed),applied_count=len(applied),applied_with_admitted_motion=admitted,applied_decisions=dict(decisions),applied=applied)
    return result


if __name__=="__main__":print(json.dumps(diagnose(Path(__file__).resolve().parent),indent=2))
