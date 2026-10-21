"""Verify mixed-clock SPP seed evidence without rerunning positioning solvers."""
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
import subprocess


def sha(data):return hashlib.sha256(data).hexdigest()
def read(path):return gzip.decompress(path.read_bytes())


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--dataset-root",type=Path);args=p.parse_args();archive=Path(__file__).resolve().parent;repo=archive.parents[2]
    m=json.loads((archive/"predeclared-inputs.json").read_text());report=json.loads((archive/"report.json").read_text());previous=repo/"docs/comparison/2026-10-21-gps-galileo-float"
    assert not m["SPP_covariance_fused"] and not m["existing_thresholds_changed"] and not m["default_seed_scope_changed"] and not m["practical_PPC_parity_accepted"] and not report["practical_PPC_parity_accepted"]
    assert m["independent_receiver_clocks"] and m["legacy_primary_group_delay"] and not m["phase_required_for_SPP"] and not m["native_whole_rerun"] and not m["held_out_validation"]
    assert len(m["cases"])==len(m["profiles"])==6 and sha((previous/"predeclared-inputs.json").read_bytes())==m["prior_manifest_sha256"]
    for path,h in m["source_sha256"].items():assert sha((repo/path).read_bytes())==h,path
    runner=runpy.run_path(str(repo/"scripts/run_ppc_multiclock_spp.py"));scoring=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"));mixed=runpy.run_path(str(repo/"scripts/run_ppc_galileo.py"))["mixed_summary"]
    for case in m["cases"]:
        city,iteration=case.split("-");references=None;ledger=None
        if args.dataset_root:
            root=args.dataset_root/city/iteration
            for name in ["rover.obs","base.obs","base.nav","reference.csv"]:
                expected=next(h for path,h in m["input_sha256"].items() if path.endswith(f"/{city}/{iteration}/{name}"));assert sha((root/name).read_bytes())==expected
            with (root/"reference.csv").open() as f:refs=[{k.strip():v.strip() for k,v in r.items()} for r in csv.DictReader(f)]
            references={(int(r["GPS Week"]),Decimal(r["GPS TOW (s)"])):r for r in refs}
        for profile,c in m["profiles"].items():
            stem=case+"."+profile;e=json.loads((archive/(stem+".execution.json")).read_text());data=read(archive/(stem+".stdout.csv.gz"));stderr=read(archive/(stem+".stderr.log.gz"))
            assert e["returncode"]==0 and sha(data)==e["csv_sha256"] and sha(stderr)==e["stderr_sha256"]
            assert b"no_SPP_covariance_fusion=true" in stderr and b"FIX_policy_evaluated=false" in stderr
            rows,actual=scoring["inspect"](data);keys=[(r["kind"],r["week"],r["tow"],r["side"]) for r in rows];assert ledger is None or ledger==keys;ledger=keys
            if not c["mixed_seed"]:
                with gzip.open(previous/(case+"."+profile+".stdout.csv.gz"),"rt") as f:assert rows==list(csv.DictReader(f))
                actual["all_columns_identical_to_prior_GPS_seed_control"]=True
            else:assert b"independent_clocks=true legacy_primary_BGD=true" in stderr
            actual["mixed_system"]=mixed(rows);actual["SPP_seed"]=runner["seed_summary"](rows,c["mixed_seed"])
            expected=report["cases"][case][profile];assert actual=={k:v for k,v in expected.items() if k not in {"reference_scores","accepted_SPP_seed_reference_scores"}}
            if references is not None:
                assert scoring["score"](rows,references,city)==expected["reference_scores"]
                seeds=[{**r,"x_m":r["seed_x_m"],"y_m":r["seed_y_m"],"z_m":r["seed_z_m"]} for r in rows if r["status"]=="FLOAT" and r["seed_source"]=="SinglePoint"]
                assert scoring["score"](seeds,references,city)==expected["accepted_SPP_seed_reference_scores"]
    a=json.loads((archive/"analysis.json").read_text());actual=runpy.run_path(str(archive/"analyze.py"))["analyze"](archive,args.dataset_root)
    if args.dataset_root is None:
        for case in a["cases"].values():
            for v in case["profiles"].values():v.pop("common_RTK_reference_scores",None)
            for v in case["GPS_vs_mixed_seed"].values():
                v.pop("common_FLOAT_reference_scores",None);v.pop("common_successful_seed_reference_scores",None)
    assert a==actual
    transfer=json.loads((archive/"clock-QR-source-transfer.json").read_text());source=subprocess.check_output(["git","show",transfer["base_commit"]+":"+transfer["source_file"]],cwd=repo,text=True)
    assert sha(source.encode())==transfer["source_sha256"] and sha((repo/transfer["new_file"]).read_bytes())==transfer["new_sha256"]
    begin=source.index("struct VelocityLinearFit {");end=source.index("#[derive(Debug, Clone, Copy, PartialEq)]\npub struct VelocityEpochConfig",begin)
    body=source[begin:end].replace("struct VelocityLinearFit {","pub(crate) struct ClockLinearFit {").replace("VelocityLinearFit","ClockLinearFit").replace("fn least_squares_with_clocks(","pub(crate) fn least_squares_with_clocks(")
    for field in ["mean","covariance","design"]:body=body.replace("    "+field+":","    pub "+field+":",1)
    new=(repo/transfer["new_file"]).read_text();new=new[new.index("pub(crate) struct ClockLinearFit {"):];assert re.sub(r"\s+","",body)==re.sub(r"\s+","",new)
    q=json.loads((archive/"qualification.json").read_text());assert q["predeclared_inputs_sha256"]==sha((archive/"predeclared-inputs.json").read_bytes()) and {c["name"] for c in q["checks"]}=={"rust","examples","fmt","clippy","release"}
    for c in q["checks"]:
        data=read(archive/c["log"]);assert c["returncode"]==0 and sha(data)==c["uncompressed_sha256"]
        if c["name"]=="rust":assert sum(map(int,re.findall(rb"test result: ok\. (\d+) passed;",data)))==c["passed_tests"]==371
    sealed=json.loads((archive/"evidence-sha256.json").read_text());assert set(sealed)=={str(p.relative_to(archive)) for p in archive.rglob("*") if p.is_file() and p.name!="evidence-sha256.json"}
    for path,h in sealed.items():assert sha((archive/path).read_bytes())==h,path
    print("Verified 36 SPP-seed runs, 18 exact GPS controls, unchanged clock QR body, all-origin comparisons and 371 Rust tests")


if __name__=="__main__":main()
