"""Verify explicit-gauge product and raw PPC admission evidence without positioning."""
import argparse
import csv
from datetime import date
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import runpy


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--verify-local-inputs",action="store_true");args=parser.parse_args()
    archive=Path(__file__).resolve().parent;repo=archive.parents[2]
    m=json.loads((archive/"predeclared-inputs.json").read_text());report=json.loads((archive/"report.json").read_text())
    assert not any(m[k] for k in ["phase_OSB_generated","positioning_executed","PPP_executed","truth_used","default_settings_changed","practical_parity_accepted"])
    assert not any(report[k] for k in ["positioning_executed","PPP_executed","practical_parity_accepted"])
    assert "unverified" in m["clock_gauge"]
    assert len(m["product_identities"])==len(report["products"])==9 and len(m["case_dates"])==len(report["cases"])==6
    for path,h in m["source_sha256"].items():assert sha(repo/path)==h,path
    if args.verify_local_inputs:
        for path,h in m["input_sha256"].items():assert sha(Path(path))==h,path
    inspect=runpy.run_path(str(repo/"scripts/run_ppc_dsb_clock.py"))["inspect"]
    total_epochs=0
    for kind,items in [("product",report["products"]),("rover",report["cases"])]:
        for stem,expected in items.items():
            e=json.loads((archive/(stem+".execution.json")).read_text());data=gzip.decompress((archive/(stem+".stdout.csv.gz")).read_bytes());err=gzip.decompress((archive/(stem+".stderr.log.gz")).read_bytes())
            assert e["returncode"]==0 and hashlib.sha256(data).hexdigest()==e["csv_sha256"] and hashlib.sha256(err).hexdigest()==e["stderr_sha256"]
            assert b"phase_OSB_generated=false" in err and b"clock_compatibility_verified=false" in err and b"PPP_executed=false" in err and b"positioning_executed=false" in err
            assert inspect(data,kind)==expected
            if kind=="rover":total_epochs+=expected["counts"]["raw_rover_epochs"]
            else:
                assert expected["rows"]==128
                rows=list(csv.DictReader(io.StringIO(data.decode())))
                for r in rows:
                    if r["status"]!="AVAILABLE":continue
                    terms=r["source_terms"].split(";");value=0.0;bound=0.0
                    for term in terms:
                        pair,coefficient,raw=term.split("*");coefficient=float(coefficient)
                        assert raw.endswith("ns") or raw.endswith("m")
                        scalar=float(raw[:-2] if raw.endswith("ns") else raw[:-1])*(299792458e-9 if raw.endswith("ns") else 1.0)
                        value+=coefficient*scalar
                    assert abs(value-float(r["code_add_m"]))<1e-12
    assert total_epochs==58634
    native=json.loads((archive/"native-generator.json").read_text())
    assert native["returncode"]==0 and native["cases"]==24 and native["code_corrections"]==96 and not native["native_DSB_clock_conversion_called"]
    assert native["source_sha256"]==sha(repo/"scripts/upstream_dsb_clock.cpp") and native["build_log_sha256"]==sha(archive/"native-generator-build.log.gz")
    assert len(native["repeat_runs"])==2 and all(r["returncode"]==0 and r["fixture_sha256"]==sha(repo/"tests/fixtures/upstream_dsb_clock.csv") for r in native["repeat_runs"])
    for path,h in native["earlier_native_source_proof"].items():assert sha(repo/path)==h,path
    for path,h in native["header_sha256"].items():
        if Path(path).exists():assert sha(Path(path))==h,path
    q=json.loads((archive/"qualification.json").read_text());assert q["predeclared_inputs_sha256"]==sha(archive/"predeclared-inputs.json") and {c["name"] for c in q["checks"]}=={"rust","examples","fmt","clippy","release"}
    for c in q["checks"]:
        data=gzip.decompress((archive/c["log"]).read_bytes());assert c["returncode"]==0 and hashlib.sha256(data).hexdigest()==c["uncompressed_sha256"]
        if c["name"]=="rust":assert sum(map(int,re.findall(rb"test result: ok\. (\d+) passed;",data)))==c["passed_tests"]==390
    sealed=json.loads((archive/"evidence-sha256.json").read_text());assert set(sealed)=={str(p.relative_to(archive)) for p in archive.rglob("*") if p.is_file() and p.name!="evidence-sha256.json"}
    for path,h in sealed.items():assert sha(archive/path)==h,path
    print("Verified 1,152 product requests, all 58,634 rover epochs, 96 independent gauge corrections and 390 Rust tests; no PPP executed")


if __name__=="__main__":main()
