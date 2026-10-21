"""Record retained SD ambiguity datum differences without changing the filter."""
import argparse
import collections
import csv
import gzip
import io
import json
from pathlib import Path
import runpy
import subprocess


def summarize(data, baseline):
    rows = list(csv.DictReader(io.StringIO(data.decode())))
    expected = {(r["kind"], r["week"], r["tow"], r["side"]): r for r in baseline}
    grouped = collections.defaultdict(list)
    for row in rows:
        key = (row["kind"], row["week"], row["tow"], row["side"])
        assert key in expected
        original = expected[key]
        assert row["status"] == original["status"]
        for field in ["x_m", "y_m", "z_m"]:
            assert row[field] == original[field]
        grouped[key].append(row)
    assert set(grouped) == set(expected)
    candidates = [r for r in rows if r["status"] == "FLOAT" and int(r["retained_arcs"]) > 0 and int(r["reset_arcs"]) > 0]
    return dict(all_epoch_statuses_and_positions_identical=True,
                paired_FLOAT_epochs=sum(r["status"] == "FLOAT" for r in baseline),
                retained_and_new_system_rows=len(candidates),
                offset_over_1_m=sum(abs(float(r["mean_offset_m"])) > 1 for r in candidates),
                offset_over_10_m=sum(abs(float(r["mean_offset_m"])) > 10 for r in candidates),
                largest_offsets=sorted(candidates,key=lambda r:abs(float(r["mean_offset_m"])),reverse=True)[:10])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root",type=Path,required=True)
    parser.add_argument("--binary",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args();repo=Path(__file__).resolve().parents[1]
    helpers=runpy.run_path(str(repo/"scripts/run_native_ppc_comparison.py"))
    sha,save,retain=(helpers[k] for k in ["sha","save","retain"])
    prior=repo/"docs/comparison/2026-10-21-gps-galileo-float"
    control=json.loads((prior/"predeclared-inputs.json").read_text())
    out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    binary=args.binary.resolve();dataset=args.dataset_root.resolve()
    sources=list((repo/"src").rglob("*.rs"))+[repo/p for p in ["examples/ppc_rtk_gauge.rs","scripts/run_ppc_gauge_probe.py"]]
    inputs={}
    for case in control["cases"]:
        city,run=case.split("-")
        for name in ["rover.obs","base.obs","base.nav"]:
            p=dataset/city/run/name
            expected=next(h for k,h in control["input_sha256"].items() if k.endswith(f"/{city}/{run}/{name}"))
            assert sha(p)==expected;inputs[str(p)]=expected
    manifest=dict(base_commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=repo,text=True).strip(),
        source_sha256={str(p.relative_to(repo)):sha(p) for p in sources},
        input_sha256=inputs,binary_sha256=sha(binary),
        prior_manifest_sha256=sha(prior/"predeclared-inputs.json"),
        cases=control["cases"],profiles=["gps","galileo"],base_ecef_m=control["base_ecef_m"],
        read_only_diagnostic=True,reference_read=False,solver_configuration_changed=False,
        interpretation="Successful updates only. Reset arcs, including global resets, are excluded from the retained-arc mean.")
    save(out/"predeclared-inputs.json",manifest);report={}
    for case in manifest["cases"]:
        city,run=case.split("-");report[case]={}
        for mode in manifest["profiles"]:
            assert all(sha(repo/p)==h for p,h in manifest["source_sha256"].items())
            assert sha(binary)==manifest["binary_sha256"]
            command=[str(binary),str(dataset/city/run),mode,",".join(map(str,manifest["base_ecef_m"][city]))]
            result=subprocess.run(command,capture_output=True,cwd=repo)
            retain(out,case+"."+mode,command,result.stdout,result.stderr,result.returncode)
            assert result.returncode==0,result.stderr.decode()
            with gzip.open(prior/(case+"."+mode+".stdout.csv.gz"),"rt") as stream:baseline=list(csv.DictReader(stream))
            report[case][mode]=summarize(result.stdout,baseline)
            save(out/"report.json",report)
            print(case,mode,report[case][mode]["retained_and_new_system_rows"],flush=True)
    assert all(sha(Path(p))==h for p,h in inputs.items())


if __name__=="__main__":main()
