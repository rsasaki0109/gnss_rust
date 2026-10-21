"""Verify 42 Rust and 12 whole-native PPC runs without executing solvers."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
import runpy


def sha(data):return hashlib.sha256(data).hexdigest()
def read(path):return gzip.decompress(path.read_bytes())
def rows(data):return list(csv.DictReader(io.StringIO(data.decode())))
def keys(data):return [(r["kind"],int(r["week"]),float(r["tow"]).hex(),r["side"]) for r in data]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root",type=Path)
    args=parser.parse_args();archive=Path(__file__).resolve().parent;repo=archive.parents[2]
    scoring=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))
    runner=runpy.run_path(str(repo/"scripts/run_ppc_galileo.py"))
    slip=runpy.run_path(str(repo/"scripts/run_ppc_slip.py"))
    prefit=runpy.run_path(str(repo/"scripts/run_ppc_prefit.py"))
    native=runpy.run_path(str(repo/"scripts/run_native_ppc_comparison.py"))
    mixed=runpy.run_path(str(repo/"scripts/run_native_galileo.py"))
    manifest=json.loads((archive/"predeclared-inputs.json").read_text());report=json.loads((archive/"report.json").read_text())
    assert not report["practical_PPC_parity_accepted"] and not manifest["thresholds_changed"]
    assert set(report["cases"])==set(manifest["cases"])
    for p,h in manifest["source_sha256"].items():assert sha((repo/p).read_bytes())==h,p
    child=archive/"native-comparison";cm=json.loads((child/"predeclared-inputs.json").read_text());cr=json.loads((child/"report.json").read_text())
    assert not cr["practical_PPC_parity_accepted"] and cm["native_AR_QC_defaults_preserved"]
    assert cm["parent_manifest_sha256"]==sha((archive/"predeclared-inputs.json").read_bytes())
    assert cm["parent_report_sha256"]==sha((archive/"report.json").read_bytes())
    for p,h in cm["source_sha256"].items():assert sha((repo/p).read_bytes())==h,p
    source=(child/"native_ppc_mixed.cpp").read_text()
    assert source==mixed["variant"]((repo/"scripts/native_ppc_rtk.cpp").read_text())
    assert sha(source.encode())==cm["adapter_sha256"] and "processor.processRTKEpoch" in source
    build=json.loads((child/"build.json").read_text());assert build["returncode"]==0
    assert sha(read(child/"build.stderr.log.gz"))==build["stderr_sha256"]
    rust_totals={p:0 for p in manifest["profiles"]};native_totals={p:0 for p in cm["seed_policies"]};input_rows=0
    for case in manifest["cases"]:
        city,run=case.split("-");references=None
        if args.dataset_root:
            root=args.dataset_root/city/run;references={}
            for name in ["rover.obs","base.obs","base.nav","reference.csv"]:
                k=next(k for k in manifest["input_sha256"] if k.endswith(f"/{city}/{run}/{name}"))
                assert sha((root/name).read_bytes())==manifest["input_sha256"][k]
            with (root/"reference.csv").open() as stream:
                for raw in csv.DictReader(stream):
                    r={k.strip():v.strip() for k,v in raw.items()};references[(int(r["GPS Week"]),Decimal(r["GPS TOW (s)"]))]=r
        ledger=None
        for profile in manifest["profiles"]:
            stem=case+"."+profile;execution=json.loads((archive/(stem+".execution.json")).read_text())
            data=read(archive/(stem+".stdout.csv.gz"));stderr=read(archive/(stem+".stderr.log.gz"))
            assert execution["returncode"]==0 and sha(data)==execution["csv_sha256"] and sha(stderr)==execution["stderr_sha256"]
            decoded,actual=scoring["inspect"](data)
            assert ledger is None or keys(decoded)==ledger;ledger=keys(decoded)
            if profile in {"gps","gps_slips"}:
                previous="2026-10-21-native-rtk-common-inputs" if profile=="gps" else "2026-10-21-ppc-slip-detection"
                suffix=".rust.stdout.csv.gz" if profile=="gps" else ".both.stdout.csv.gz"
                assert decoded==rows(read(archive.parent/previous/(case+suffix)))
                actual["all_columns_identical_to_prior_control"]=True
            else:
                actual["mixed_system"]=runner["mixed_summary"](decoded)
                actual["slip_evidence"]=slip["slip_summary"](decoded)
                if profile.endswith("_gate"):actual["prefit"]=prefit["prefit_summary"](decoded)
            expected=report["cases"][case][profile]
            assert actual=={k:v for k,v in expected.items() if k!="reference_scores"}
            if references is not None:assert scoring["score"](decoded,references,city)==expected["reference_scores"]
            rust_totals[profile]+=actual["status_counts"].get("FLOAT",0)
        comparisons=[]
        for policy in cm["seed_policies"]:
            stem=case+"."+policy;execution=json.loads((child/(stem+".execution.json")).read_text())
            data=read(child/(stem+".stdout.csv.gz"));stderr=read(child/(stem+".stderr.log.gz"))
            assert execution["returncode"]==0 and sha(data)==execution["csv_sha256"] and sha(stderr)==execution["stderr_sha256"]
            decoded,rtk,actual=native["native_summary"](data)
            assert keys(decoded)==ledger and actual["valid_outputs_with_mismatched_time"]==0
            assert b"input_drained=1" in stderr and b"native_whole_RTKProcessor=true" in stderr
            expected=cr["cases"][case][policy]["native"]
            assert actual=={k:v for k,v in expected.items() if k not in {"reference_scores","FIX_reference_scores"}}
            if references is not None:
                assert scoring["score"](rtk,references,city)==expected["reference_scores"]
                assert scoring["score"]([{**r,"status":"FLOAT"} for r in decoded if r["status"]=="FIX"],references,city)==expected["FIX_reference_scores"]
            native_totals[policy]+=actual["valid_RTK_positions"]
            comparisons.append(cr["cases"][case][policy]["input_comparison"])
        assert comparisons[0]==comparisons[1]
        native_ledger=read(child/(case+".native.input-ledger.csv.gz"));rust_ledger=read(child/(case+".rust.input-ledger.csv.gz"))
        comparison=comparisons[0]
        assert sha(native_ledger)==comparison["native_ledger_sha256"] and sha(rust_ledger)==comparison["rust_ledger_sha256"]
        from itertools import zip_longest
        digest=hashlib.sha256();count=0
        for count,(a,b) in enumerate(zip_longest(rows(native_ledger),rows(rust_ledger)),1):
            assert a is not None and b is not None
            value=mixed["ledger_key"](a,native["ledger_key"])
            assert value==mixed["ledger_key"](b,native["ledger_key"])
            digest.update((json.dumps(value,separators=(",",":"))+"\n").encode())
        assert count==comparison["rows"] and digest.hexdigest()==comparison["normalized_sha256"]
        input_rows+=count
    qualification=json.loads((archive/"qualification.json").read_text())
    assert sha((archive/"predeclared-inputs.json").read_bytes())==qualification["predeclared_inputs_sha256"]
    assert {c["name"] for c in qualification["checks"]}=={"rust","examples","fmt","clippy","release"}
    for check in qualification["checks"]:
        data=read(archive/check["log"]);assert check["returncode"]==0 and sha(data)==check["uncompressed_sha256"]
        if check["name"]=="rust":assert sum(map(int,re.findall(rb"test result: ok\. (\d+) passed;",data)))==check["passed_tests"]==316
    helper=json.loads((child/"helper-checks.json").read_text())
    assert {c["name"] for c in helper["checks"]}=={"fmt","clippy","ledger-release"}
    for c in helper["checks"]:assert c["returncode"]==0 and sha(read(child/c["log"]))==c["uncompressed_sha256"]
    component=json.loads((archive/"component-generation.json").read_text())
    assert (component["geometry_cases"],component["moving_epochs"],component["startup_cases"])==(24,32,3)
    for c in component["records"]:
        assert c["returncode"]==0 and c["regenerated_bytes_identical"]
        for f,h in c["fixtures_sha256"].items():assert sha((repo/"tests/fixtures"/f).read_bytes())==h
    smoke=json.loads((child/"synthetic-smoke/report.json").read_text());truth={}
    with (repo/"tests/fixtures/synthetic_rtk_mixed_truth.csv").open() as stream:
        for r in csv.DictReader(stream):truth[(int(r["week"]),float(r["tow"]))]=[float(r[k]) for k in ["x_m","y_m","z_m"]]
    for policy,c in smoke["policies"].items():
        assert c["returncode"]==0 and c["input_comparison"]["rows"]==512
        decoded,_,actual=native["native_summary"](read(child/("synthetic-smoke/"+policy+".output.csv.gz")))
        errors=[math.dist([float(r[k]) for k in ["x_m","y_m","z_m"]],truth[(int(r["week"]),float(r["tow"]))]) for r in decoded if r["status"] in {"FLOAT","FIX"}]
        actual["maximum_3d_truth_error_m"]=max(errors)
        assert actual==c["summary"] and actual["paired_epochs"]==32 and actual["valid_RTK_positions"]==32
    diagnostic=archive/"input-diagnostics";records=json.loads((diagnostic/"executions.json").read_text())
    for c in records["inventory"]:assert c["returncode"]==0 and sha(read(diagnostic/(c["case"]+".inventory.csv.gz")))==c["stdout_sha256"]
    assert records["spp_probe"]["returncode"]==0 and sha(read(diagnostic/"nagoya-run2.spp-worst.csv.gz"))==records["spp_probe"]["stdout_sha256"]
    if args.dataset_root:
        analyze=runpy.run_path(str(archive/"analyze.py"))
        assert analyze["diagnostics"](archive,args.dataset_root)==json.loads((archive/"shared-epoch-diagnostics.json").read_text())
    for p,h in json.loads((archive/"evidence-sha256.json").read_text()).items():assert sha((archive/p).read_bytes())==h,p
    print("Verified 42 Rust + 12 whole-native PPC runs; input rows:",input_rows,"; 316 tests.")
    print("Rust FLOAT:",rust_totals,"; native RTK:",native_totals)


if __name__=="__main__":main()
