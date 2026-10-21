"""Recheck frozen datum-alignment evidence without running a GNSS solver."""
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
def rows(data):return list(csv.DictReader(io.StringIO(data.decode())))
def key(row):return (row["kind"],row["week"],row["tow"],row["side"])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root",type=Path)
    args=parser.parse_args();archive=Path(__file__).resolve().parent;repo=archive.parents[2]
    previous=repo/"docs/comparison/2026-10-21-gps-galileo-float"
    scoring=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))
    alignment=runpy.run_path(str(repo/"scripts/run_ppc_datum_alignment.py"))
    probe=runpy.run_path(str(repo/"scripts/run_ppc_gauge_probe.py"))
    run=archive/"ppc-comparison";manifest=json.loads((run/"predeclared-inputs.json").read_text())
    report=json.loads((run/"report.json").read_text())
    assert not report["practical_PPC_parity_accepted"] and not manifest["default_changed"] and not manifest["thresholds_changed"]
    assert sha((previous/"predeclared-inputs.json").read_bytes())==manifest["prior_manifest_sha256"]
    assert sha((previous/"report.json").read_bytes())==manifest["prior_report_sha256"]
    for path,h in manifest["source_sha256"].items():assert sha((repo/path).read_bytes())==h,path
    diagnostic=archive/"read-only-probe";pm=json.loads((diagnostic/"predeclared-inputs.json").read_text());pr=json.loads((diagnostic/"report.json").read_text())
    assert pm["read_only_diagnostic"] and not pm["reference_read"] and not pm["solver_configuration_changed"]
    for path,h in pm["source_sha256"].items():
        data=subprocess.check_output(["git","show",pm["base_commit"]+":"+path],cwd=repo) if path.startswith("src/") else (repo/path).read_bytes()
        assert sha(data)==h,path
    counts={p:0 for p in manifest["profiles"]}
    for case in manifest["cases"]:
        references=None;city,iteration=case.split("-")
        if args.dataset_root:
            root=args.dataset_root/city/iteration;references={}
            for name in ["rover.obs","base.obs","base.nav","reference.csv"]:
                path=root/name;expected=next(h for k,h in manifest["input_sha256"].items() if k.endswith(f"/{city}/{iteration}/{name}"))
                assert sha(path.read_bytes())==expected
            with (root/"reference.csv").open() as stream:
                for raw in csv.DictReader(stream):
                    r={k.strip():v.strip() for k,v in raw.items()};references[(int(r["GPS Week"]),Decimal(r["GPS TOW (s)"]))]=r
        ledger=None
        for profile in manifest["profiles"]:
            stem=case+"."+profile;execution=json.loads((run/(stem+".execution.json")).read_text())
            data=read(run/(stem+".stdout.csv.gz"));stderr=read(run/(stem+".stderr.log.gz"))
            assert execution["returncode"]==0 and sha(data)==execution["csv_sha256"] and sha(stderr)==execution["stderr_sha256"]
            decoded,actual=scoring["inspect"](data)
            assert ledger is None or list(map(key,decoded))==ledger;ledger=list(map(key,decoded))
            if profile.endswith("_aligned"):
                actual["datum_alignment"]=alignment["alignment_summary"](decoded)
                assert b"ambiguity_datum_alignment=per_system_retained_mean" in stderr
            else:
                assert decoded==rows(read(previous/(stem+".stdout.csv.gz")))
                actual["all_columns_identical_to_prior_control"]=True
            expected=report["cases"][case][profile]
            assert actual=={k:v for k,v in expected.items() if k!="reference_scores"}
            if references is not None:assert scoring["score"](decoded,references,city)==expected["reference_scores"]
            counts[profile]+=actual["status_counts"].get("FLOAT",0)
        for profile in pm["profiles"]:
            stem=case+"."+profile;execution=json.loads((diagnostic/(stem+".execution.json")).read_text())
            data=read(diagnostic/(stem+".stdout.csv.gz"));stderr=read(diagnostic/(stem+".stderr.log.gz"))
            assert execution["returncode"]==0 and sha(data)==execution["csv_sha256"] and sha(stderr)==execution["stderr_sha256"]
            assert probe["summarize"](data,rows(read(previous/(stem+".stdout.csv.gz"))))==pr[case][profile]
    qualification=json.loads((archive/"qualification.json").read_text())
    assert {c["name"] for c in qualification["checks"]}=={"rust","examples","fmt","clippy","release"}
    for c in qualification["checks"]:
        data=read(archive/c["log"]);assert c["returncode"]==0 and sha(data)==c["uncompressed_sha256"]
        if c["name"]=="rust":assert sum(map(int,re.findall(rb"test result: ok\. (\d+) passed;",data)))==c["passed_tests"]==324
    component=json.loads((archive/"component-generation.json").read_text())
    assert component["compile_returncode"]==component["run_returncode"]==0 and component["regenerated_bytes_identical"]
    assert sha((repo/"scripts/upstream_rtk_datum.cpp").read_bytes())==component["generator_sha256"]
    assert sha((repo/"tests/fixtures/upstream_rtk_datum.csv").read_bytes())==component["fixture_sha256"]
    assert sha(read(archive/"component-build.log.gz"))==component["compile_log_sha256"]
    if args.dataset_root:
        analysis=runpy.run_path(str(archive/"analyze.py"))
        assert analysis["diagnostics"](archive,args.dataset_root)==json.loads((archive/"shared-epoch-diagnostics.json").read_text())
    for path,h in json.loads((archive/"evidence-sha256.json").read_text()).items():assert sha((archive/path).read_bytes())==h,path
    print("Verified 84 PPC runs, 12 read-only probes, native updateBias cases and 324 tests.")
    print("FLOAT totals:",counts)


if __name__=="__main__":main()
