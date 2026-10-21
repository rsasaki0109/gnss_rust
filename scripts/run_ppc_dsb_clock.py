"""Audit explicit GPS IF DSB gauges and all six PPC raw dual-frequency streams."""
import argparse
import collections
import csv
from datetime import date
import gzip
import hashlib
import io
import json
from pathlib import Path
import runpy
import subprocess


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def inspect(data,kind):
    rows=list(csv.DictReader(io.StringIO(data.decode())))
    if kind=="product":
        assert len(rows)==128 and len({(r["satellite"],r["tracking"]) for r in rows})==128
        statuses=collections.Counter(r["status"] for r in rows)
        for r in rows:
            if r["status"]=="AVAILABLE":
                assert r["code_add_m"] and r["sigma_upper_bound_m"] and r["source_terms"] and not r["error"]
            else:assert not r["code_add_m"] and not r["sigma_upper_bound_m"] and r["error"]
        return dict(rows=len(rows),status_counts=dict(statuses))
    assert len({(r["week"],r["tow"]) for r in rows})==len(rows)
    stats=collections.Counter();epochs=collections.Counter()
    for r in rows:
        stats["raw_rover_epochs"]+=1;epochs[r["kind"]]+=1
        for pair in ["1c2w","1c2l"]:
            values={k:int(r[k+"_"+pair]) for k in ["complete","unflagged","dsb_available","dsb_missing","dsb_expired","dsb_ambiguous","dsb_other"]}
            assert 0<=values["unflagged"]<=values["complete"]<=int(r["gps_satellites"])
            assert values["complete"]==sum(values[k] for k in ["dsb_available","dsb_missing","dsb_expired","dsb_ambiguous","dsb_other"])
            for k,v in values.items():stats[k+"_"+pair]+=v
            stats["epochs_with_at_least_four_complete_"+pair]+=values["complete"]>=4
            stats["epochs_with_at_least_four_unflagged_"+pair]+=values["unflagged"]>=4
            stats["epochs_with_at_least_four_DSB_available_"+pair]+=values["dsb_available"]>=4
        stats["missing_LLI_pair_rows"]+=int(r["missing_lli_pair_rows"])
    return dict(counts=dict(stats),epoch_kinds=dict(epochs),scope="all original rover epochs; pair counts can overlap, DSB is code only")


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--dataset-root",type=Path,required=True);p.add_argument("--product-root",type=Path,required=True);p.add_argument("--binary",type=Path,required=True);p.add_argument("--output-dir",type=Path,required=True);a=p.parse_args()
    repo=Path(__file__).resolve().parents[1];dataset=a.dataset_root.resolve();product=a.product_root.resolve();binary=a.binary.resolve();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    acquisition=json.loads((product/"acquisition.json").read_text());items=[p for p in acquisition["products"] if p["kind"]=="dsb"]
    assert len(items)==9 and all(p["acquired"] for p in items)
    sources=list((repo/"src").rglob("*.rs"))+[repo/f for f in ["scripts/run_ppc_dsb_clock.py","scripts/upstream_dsb_clock.cpp","examples/ppc_dsb_clock.rs","tests/dsb_clock.rs","tests/fixtures/upstream_dsb_clock.csv","docs/dsb-clock-gauge.md","Cargo.toml","Cargo.lock","rust-toolchain.toml"]]
    dates={"nagoya-run1":"2024-08-03","nagoya-run2":"2024-07-20","nagoya-run3":"2024-08-03","tokyo-run1":"2024-07-23","tokyo-run2":"2024-07-23","tokyo-run3":"2024-07-23"}
    inputs={str(product/"acquisition.json"):sha(product/"acquisition.json")}
    for item in items:
        path=product/item["expanded_filename"];assert sha(path)==item["expanded_sha256"];inputs[str(path)]=sha(path)
    for case in dates:
        city,iteration=case.split("-");path=dataset/city/iteration/"rover.obs";inputs[str(path)]=sha(path)
    manifest=dict(base_commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=repo,text=True).strip(),source_sha256={str(s.relative_to(repo)):sha(s) for s in sources},input_sha256=inputs,binary_sha256=sha(binary),product_identities=items,case_dates=dates,clock_gauge="explicit IF(C1W,C2W); compatibility unverified",phase_OSB_generated=False,positioning_executed=False,PPP_executed=False,truth_used=False,default_settings_changed=False,practical_parity_accepted=False)
    helpers=runpy.run_path(str(repo/"scripts/run_native_ppc_comparison.py"));save,retain=(helpers[k] for k in ["save","retain"])
    save(out/"predeclared-inputs.json",manifest);report=dict(products={},cases={},positioning_executed=False,PPP_executed=False,practical_parity_accepted=False)
    def execute(stem,command,kind):
        assert all(sha(repo/path)==h for path,h in manifest["source_sha256"].items()) and sha(binary)==manifest["binary_sha256"]
        r=subprocess.run(command,capture_output=True,cwd=repo);retain(out,stem,command,r.stdout,r.stderr,r.returncode);assert r.returncode==0,r.stderr.decode();metrics=inspect(r.stdout,kind);print(stem,metrics,flush=True);return metrics
    for item in items:
        day=date.fromisoformat(item["date"]);days=(day-date(1980,1,6)).days;week,tow=days//7,(days%7)*86400+43200
        report["products"][item["date"]]=execute(item["date"],[str(binary),str(product/item["expanded_filename"]),"--at",str(week),str(tow)],"product");save(out/"report.json",report)
    for case,day in dates.items():
        city,iteration=case.split("-");item=next(p for p in items if p["date"]==day)
        report["cases"][case]=execute(case,[str(binary),str(product/item["expanded_filename"]),"--rover",str(dataset/city/iteration/"rover.obs")],"rover");save(out/"report.json",report)
    assert all(sha(Path(path))==h for path,h in inputs.items())


if __name__=="__main__":main()
