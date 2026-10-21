"""Freeze GPS/Galileo within-system DD FLOAT comparisons on all PPC runs."""
import argparse
import collections
import csv
from decimal import Decimal
import gzip
import json
from pathlib import Path
import runpy
import subprocess


def save(path,value):path.write_text(json.dumps(value,indent=2)+"\n")


def mixed_summary(rows):
    counts=collections.Counter();tracking=collections.Counter();seeds=collections.Counter();failures=collections.Counter()
    for r in rows:
        if r["status"]!="FLOAT":continue
        refs=r["reference_satellites"].split(";");used=r["used_satellites"].split(";")
        choices=r["tracking_pairs"].split(";")
        assert len(used)==int(r["satellites"]) and len(choices)==len(used)
        assert set(refs)<=set(used) and len({s[0] for s in refs})==len(refs)
        assert set(r["geometry_clock_systems"].split(";"))=={s[0] for s in refs}
        assert int(r["dd_phase_rows"])==len(used)-len(refs)>=3
        parsed=[s.split(":") for s in choices]
        assert {p[0] for p in parsed}==set(used)
        for sat,ro,ba in parsed:
            assert sat[0] in {"G","E"} and ro in {"1C","1X"} and ba in {"1C","1X"}
            if sat[0]=="G":assert ro==ba=="1C"
            tracking[sat[0]+":"+ro+":"+ba]+=1
        counts["two_system_epochs"]+=len(refs)==2
        counts["single_system_epochs"]+=len(refs)==1
        counts["GPS_satellite_rows"]+=sum(s[0]=="G" for s in used)
        counts["Galileo_satellite_rows"]+=sum(s[0]=="E" for s in used)
        seeds[r["seed_source"]]+=1
        if r["spp_failure"]:failures[r["spp_failure"]]+=1
    return dict(counts=dict(counts),tracking_pairs=dict(tracking),accepted_seed_sources=dict(seeds),accepted_SPP_failures=dict(failures))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root",type=Path,required=True)
    parser.add_argument("--binary",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args();repo=Path(__file__).resolve().parents[1]
    baseline=repo/"docs/comparison/2026-10-21-native-rtk-common-inputs"
    slip_baseline=repo/"docs/comparison/2026-10-21-ppc-slip-detection"
    control=json.loads((baseline/"predeclared-inputs.json").read_text())
    scoring=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"));sha=scoring["sha"]
    slip=runpy.run_path(str(repo/"scripts/run_ppc_slip.py"))
    prefit=runpy.run_path(str(repo/"scripts/run_ppc_prefit.py"))
    capture=runpy.run_path(str(repo/"scripts/run_native_ppc_comparison.py"))["retain"]
    output,dataset,binary=args.output_dir.resolve(),args.dataset_root.resolve(),args.binary.resolve()
    output.mkdir(parents=True,exist_ok=False)
    sources=list((repo/"src").rglob("*.rs"))+[repo/n for n in [
        "scripts/run_ppc_galileo.py","scripts/run_ppc_rtk_float.py","scripts/run_ppc_slip.py","scripts/run_ppc_prefit.py","scripts/run_native_ppc_comparison.py",
        "examples/ppc_rtk_float.rs","examples/ppc_signal_inventory.rs","examples/ppc_spp_diagnostics.rs",
        "tests/rtk_multisystem.rs","tests/geometry.rs","scripts/upstream_rtk_mixed.cpp","scripts/upstream_rtk_mixed_startup.cpp",
        "scripts/upstream_multiclock_geometry.cpp","scripts/upstream_navigation.cpp","scripts/upstream_rtk_float.cpp",
        "tests/fixtures/synthetic_rtk_mixed.nav","tests/fixtures/synthetic_rtk_mixed_rover.obs","tests/fixtures/synthetic_rtk_mixed_base.obs",
        "tests/fixtures/synthetic_rtk_mixed_truth.csv","tests/fixtures/upstream_rtk_mixed_startup.csv","tests/fixtures/upstream_multiclock_geometry.csv",
        "Cargo.toml","Cargo.lock","rust-toolchain.toml",
    ]]
    inputs={}
    for case in control["cases"]:
        city,run=case.split("-")
        for name in ["rover.obs","base.obs","base.nav","reference.csv"]:
            p=dataset/city/run/name
            k=next(k for k in control["input_sha256"] if k.endswith(f"/{city}/{run}/{name}"))
            assert sha(p)==control["input_sha256"][k];inputs[str(p)]=sha(p)
    manifest=dict(source_sha256={str(p.relative_to(repo)):sha(p) for p in sources},input_sha256=inputs,
        binary_sha256=sha(binary),cases=control["cases"],base_ecef_m=control["base_ecef_m"],
        profiles={"gps":"kinematic","gps_slips":"kinematic-slips","galileo":"kinematic-galileo",
                  "galileo_slips":"kinematic-galileo-slips","galileo_spp":"kinematic-galileo-spp",
                  "galileo_spp_slips":"kinematic-galileo-spp-slips","galileo_spp_slips_gate":"kinematic-galileo-spp-slips-prefit3"},
        Galileo_tracking="complete E1 group: prefer 1C, otherwise 1X independently per receiver; tracking changes restart arc",
        clocks="within-system DD reference and R blocks; no cross-system clock alignment",SPP_seed="GPS 1C only, unchanged",
        thresholds_changed=False,reference_used_to_choose_configuration=False,
        mixed_constellation_FIX_enabled=False,native_whole_solver_rerun_here=False,
        timing_memory_acceptance=False,practical_PPC_parity_accepted=False)
    save(output/"predeclared-inputs.json",manifest)
    report=dict(practical_PPC_parity_accepted=False,cases={})
    for case in control["cases"]:
        city,run=case.split("-");root=dataset/city/run;references={}
        with (root/"reference.csv").open() as stream:
            for raw in csv.DictReader(stream):
                r={k.strip():v.strip() for k,v in raw.items()};references[(int(r["GPS Week"]),Decimal(r["GPS TOW (s)"]))]=r
        report["cases"][case]={};ledger=None
        for profile,mode in manifest["profiles"].items():
            assert all(sha(repo/p)==h for p,h in manifest["source_sha256"].items())
            assert sha(binary)==manifest["binary_sha256"]
            command=[str(binary),str(root),mode,",".join(map(str,control["base_ecef_m"][city]))]
            result=subprocess.run(command,capture_output=True,cwd=repo)
            capture(output,case+"."+profile,command,result.stdout,result.stderr,result.returncode)
            assert result.returncode==0,result.stderr.decode()
            rows,metrics=scoring["inspect"](result.stdout)
            keys=[(r["kind"],int(r["week"]),float(r["tow"]).hex(),r["side"]) for r in rows]
            assert ledger is None or ledger==keys;ledger=keys
            if profile in {"gps","gps_slips"}:
                path=baseline/(case+".rust.stdout.csv.gz") if profile=="gps" else slip_baseline/(case+".both.stdout.csv.gz")
                with gzip.open(path,"rt") as stream:assert rows==list(csv.DictReader(stream))
                metrics["all_columns_identical_to_prior_control"]=True
            else:
                metrics["mixed_system"]=mixed_summary(rows)
                metrics["slip_evidence"]=slip["slip_summary"](rows)
                if profile.endswith("_gate"):metrics["prefit"]=prefit["prefit_summary"](rows)
            metrics["reference_scores"]=scoring["score"](rows,references,city)
            report["cases"][case][profile]=metrics;save(output/"report.json",report)
            print(case,profile,metrics["status_counts"].get("FLOAT",0),flush=True)
    assert all(sha(Path(p))==h for p,h in inputs.items())
    assert all(sha(repo/p)==h for p,h in manifest["source_sha256"].items())


if __name__=="__main__":main()
