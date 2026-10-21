"""Read-only verification of frozen prefit admission experiments."""
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


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root",type=Path)
    args=parser.parse_args()
    archive=Path(__file__).resolve().parent
    repo=archive.parents[2]
    baseline=archive.parent/"2026-10-21-native-rtk-common-inputs"
    spp_baseline=archive.parent/"2026-10-21-rust-spp-seed"
    helpers=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))
    prefit=runpy.run_path(str(repo/"scripts/run_ppc_prefit.py"))
    manifest=json.loads((archive/"predeclared-inputs.json").read_text())
    report=json.loads((archive/"report.json").read_text())
    assert not report["practical_PPC_parity_accepted"] and not manifest["other_thresholds_changed"]
    assert manifest["max_prefit_residual_m"]==3.0
    assert set(report["cases"])==set(manifest["cases"])
    for name,expected in manifest["source_sha256"].items():assert sha((repo/name).read_bytes())==expected,name
    totals={p:0 for p in manifest["profiles"]}
    worst=[]
    for case in manifest["cases"]:
        city,run=case.split("-")
        references=None
        if args.dataset_root:
            root=args.dataset_root/city/run
            for name in ["rover.obs","base.obs","base.nav","reference.csv"]:
                key=next(k for k in manifest["input_sha256"] if k.endswith(f"/{city}/{run}/{name}"))
                assert sha((root/name).read_bytes())==manifest["input_sha256"][key]
            references={}
            with (root/"reference.csv").open() as stream:
                for raw in csv.DictReader(stream):
                    r={k.strip():v.strip() for k,v in raw.items()}
                    references[(int(r["GPS Week"]),Decimal(r["GPS TOW (s)"]))]=r
        ledger=None
        for profile in manifest["profiles"]:
            stem=case+"."+profile
            execution=json.loads((archive/(stem+".execution.json")).read_text())
            data=gzip.decompress((archive/(stem+".stdout.csv.gz")).read_bytes())
            stderr=gzip.decompress((archive/(stem+".stderr.log.gz")).read_bytes())
            assert execution["returncode"]==0 and sha(data)==execution["csv_sha256"] and sha(stderr)==execution["stderr_sha256"]
            rows,actual=helpers["inspect"](data)
            keys=[(r["kind"],int(r["week"]),float(r["tow"]).hex(),r["side"]) for r in rows]
            assert ledger is None or ledger==keys
            ledger=keys
            if profile in {"none","spp"}:
                path=(baseline/(case+".rust.stdout.csv.gz") if profile=="none" else spp_baseline/(case+".spp.stdout.csv.gz"))
                with gzip.open(path,"rt") as stream:assert rows==list(csv.DictReader(stream))
                actual["all_columns_identical_to_prior_control"]=True
            else:actual["prefit"]=prefit["prefit_summary"](rows)
            expected=report["cases"][case][profile]
            assert actual=={k:v for k,v in expected.items() if k!="reference_scores"}
            totals[profile]+=actual["status_counts"].get("FLOAT",0)
            if references is not None:
                assert helpers["score"](rows,references,city)==expected["reference_scores"]
                if case=="nagoya-run2" and profile=="spp_gate":
                    for r in rows:
                        if r["status"]!="FLOAT":continue
                        t=references[(int(r["week"]),Decimal(r["tow"]))]
                        xyz=[float(t[k]) for k in ["ECEF X (m)","ECEF Y (m)","ECEF Z (m)"]]
                        worst.append(dict(error_3d_m=math.dist(xyz,[float(r[k]) for k in ["x_m","y_m","z_m"]]),
                                          seed_error_3d_m=math.dist(xyz,[float(r[k]) for k in ["seed_x_m","seed_y_m","seed_z_m"]]),output_row=r))
    if args.dataset_root:
        expected=json.loads((archive/"worst-spp-gated-nagoya2.json").read_text())
        assert expected["post_run_only"] and expected["reference_hypothesis"]=="antenna"
        assert sorted(worst,key=lambda r:r["error_3d_m"],reverse=True)[:10]==expected["positions"]
    native=json.loads((archive/"native-prefit-generation.json").read_text())
    assert native["returncode"]==0 and native["cases"]==10 and native["generated_fixture_exactly_matches"]
    assert sha((repo/"tests/fixtures/upstream_rtk_prefit.csv").read_bytes())==native["fixture_sha256"]
    assert sha((repo/"scripts/upstream_rtk_prefit.cpp").read_bytes())==native["source_sha256"]
    q=json.loads((archive/"qualification.json").read_text())
    assert sha((archive/"predeclared-inputs.json").read_bytes())==q["predeclared_inputs_sha256"]
    assert {c["name"] for c in q["checks"]}=={"rust","examples","fmt","clippy","release"}
    for check in q["checks"]:
        data=gzip.decompress((archive/check["log"]).read_bytes())
        assert check["returncode"]==0 and sha(data)==check["uncompressed_sha256"]
        if check["name"]=="rust":
            assert sum(map(int,re.findall(rb"test result: ok\. (\d+) passed;",data)))==check["passed_tests"]==302
    for name,expected in json.loads((archive/"evidence-sha256.json").read_text()).items():assert sha((archive/name).read_bytes())==expected,name
    print("Verified 36 runs, unchanged RTK/SPP control CSVs, 10 native classifications, 302 tests; FLOAT:",totals)


if __name__=="__main__":main()
