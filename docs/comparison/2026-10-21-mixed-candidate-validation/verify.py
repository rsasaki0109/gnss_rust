"""Recheck mixed-candidate evidence without executing GNSS solvers."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import runpy


def sha(data):return hashlib.sha256(data).hexdigest()
def read(path):return gzip.decompress(path.read_bytes())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root",type=Path)
    args=parser.parse_args();archive=Path(__file__).resolve().parent;repo=archive.parents[2]
    previous=repo/"docs/comparison/2026-10-21-gps-galileo-float"
    runner=runpy.run_path(str(repo/"scripts/run_ppc_mixed_candidates.py"))
    scoring=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))
    manifest=json.loads((archive/"predeclared-inputs.json").read_text());report=json.loads((archive/"report.json").read_text())
    assert not report["practical_PPC_parity_accepted"] and not manifest["FIX_policy_enabled"] and not manifest["thresholds_changed"]
    assert sha((previous/"predeclared-inputs.json").read_bytes())==manifest["prior_manifest_sha256"]
    for path,h in manifest["source_sha256"].items():assert sha((repo/path).read_bytes())==h,path
    totals={p:{} for p in manifest["diagnostic_modes"]}
    for case in manifest["cases"]:
        references=None;city,iteration=case.split("-")
        if args.dataset_root:
            root=args.dataset_root/city/iteration;references={}
            for name in ["rover.obs","base.obs","base.nav","reference.csv"]:
                expected=next(h for k,h in manifest["input_sha256"].items() if k.endswith(f"/{city}/{iteration}/{name}"))
                assert sha((root/name).read_bytes())==expected
            with (root/"reference.csv").open() as stream:
                for raw in csv.DictReader(stream):
                    r={k.strip():v.strip() for k,v in raw.items()};references[(int(r["GPS Week"]),Decimal(r["GPS TOW (s)"]))]=r
        for profile in manifest["controls"]:
            def execution(suffix):
                stem=case+"."+profile+"."+suffix
                record=json.loads((archive/(stem+".execution.json")).read_text())
                data=read(archive/(stem+".stdout.csv.gz"));stderr=read(archive/(stem+".stderr.log.gz"))
                assert record["returncode"]==0 and sha(data)==record["csv_sha256"] and sha(stderr)==record["stderr_sha256"]
                return data,stderr
            data,_=execution("control");control,metrics=scoring["inspect"](data)
            assert control==list(csv.DictReader(io.StringIO(read(previous/(case+"."+profile+".stdout.csv.gz")).decode())))
            metrics["all_columns_identical_to_prior_control"]=True
            assert metrics==report["cases"][case][profile]["control"]
            if profile in manifest["diagnostic_modes"]:
                data,stderr=execution("diagnostic");actual=runner["candidate_summary"](data,control)
                assert b"FIX_policy_evaluated=false" in stderr and b"no_feedback_to_FLOAT=true" in stderr
                expected=report["cases"][case][profile]["diagnostic"]
                assert actual=={k:v for k,v in expected.items() if k!="candidate_reference_scores"}
                if references is not None:
                    decoded=list(csv.DictReader(io.StringIO(data.decode())))
                    assert runner["candidate_scores"](decoded,references,city,scoring["score"])==expected["candidate_reference_scores"]
                for k,v in actual["counts"].items():totals[profile][k]=totals[profile].get(k,0)+v
    qualification=json.loads((archive/"qualification.json").read_text())
    assert qualification["predeclared_inputs_sha256"]==sha((archive/"predeclared-inputs.json").read_bytes())
    assert {c["name"] for c in qualification["checks"]}=={"rust","examples","fmt","clippy","release"}
    for c in qualification["checks"]:
        data=read(archive/c["log"]);assert c["returncode"]==0 and sha(data)==c["uncompressed_sha256"]
        if c["name"]=="rust":assert sum(map(int,re.findall(rb"test result: ok\. (\d+) passed;",data)))==c["passed_tests"]==330
    if args.dataset_root:
        analyze=runpy.run_path(str(archive/"analyze.py"))
        assert analyze["diagnostics"](archive,args.dataset_root)==json.loads((archive/"eligibility-diagnostics.json").read_text())
    for path,h in json.loads((archive/"evidence-sha256.json").read_text()).items():assert sha((archive/path).read_bytes())==h,path
    print("Verified 18 unchanged FLOAT controls + 12 read-only candidate runs; 330 tests; no FIX declared.")
    print(totals)


if __name__=="__main__":main()
