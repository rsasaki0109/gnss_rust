"""Post-run diagnostics only; no results are fed back to either solver."""
import argparse
import csv
from decimal import Decimal
import gzip
import json
import math
from pathlib import Path
import runpy


def time_key(row):
    return int(row["week"]), float(row["tow"]).hex()


def analyze(archive, dataset):
    repo = archive.parents[2]
    scoring = runpy.run_path(str(repo / "scripts/run_ppc_rtk_float.py"))
    report = json.loads((archive / "report.json").read_text())
    result = dict(post_run_only=True, acceptance_thresholds_defined=False, cases={})
    for case, original in report["cases"].items():
        city, run = case.split("-")
        references = {}
        with (dataset / city / run / "reference.csv").open() as stream:
            for raw in csv.DictReader(stream):
                row = {k.strip(): v.strip() for k, v in raw.items()}
                references[(int(row["GPS Week"]), Decimal(row["GPS TOW (s)"]))] = row
        rows = {}
        for solver in ["native", "rust"]:
            with gzip.open(archive / (case + "." + solver + ".stdout.csv.gz"), "rt") as stream:
                rows[solver] = [r for r in csv.DictReader(stream) if r["kind"] == "pair"]
        native = {time_key(r): r for r in rows["native"] if r["status"] in {"FLOAT", "FIX"}}
        rust = {time_key(r): r for r in rows["rust"] if r["status"] == "FLOAT"}
        common = native.keys() & rust.keys()
        # Status remains intact in the retained CSV. The shared scoring helper
        # selects FLOAT, so use copies solely for native FIX accuracy statistics.
        def score(selected):
            return scoring["score"]([{**r, "status": "FLOAT"} for r in selected], references, city)
        delta = [math.sqrt(sum((float(native[t][k]) - float(rust[t][k])) ** 2
                               for k in ["x_m", "y_m", "z_m"])) for t in sorted(common)]
        slips = {k: sum(int(r[k]) for r in rows["native"])
                 for k in ["code_slips", "doppler_slips", "lli_slips"]}
        result["cases"][case] = dict(
            paired_epochs=original["native"]["paired_epochs"],
            native_RTK_only_epochs=len(native.keys() - rust.keys()),
            rust_FLOAT_only_epochs=len(rust.keys() - native.keys()),
            common_RTK_epochs=len(common),
            neither_RTK_epochs=len(rows["native"]) - len(native.keys() | rust.keys()),
            native_RTK_coverage_among_pairs=len(native) / len(rows["native"]),
            native_status_conditioned_reference_scores={
                status: score(r for r in rows["native"] if r["status"] == status)
                for status in ["FLOAT", "FIX", "SPP"]
            },
            common_epochs_reference_scores={
                "native": score(native[t] for t in sorted(common)),
                "rust": score(rust[t] for t in sorted(common)),
            },
            common_epochs_solver_position_difference_m=dict(
                epochs=len(delta), rms=math.sqrt(sum(d*d for d in delta) / len(delta)), maximum=max(delta)
            ) if delta else None,
            native_slip_telemetry_counts=slips,
        )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root", type=Path, required=True)
    args = parser.parse_args()
    archive = Path(__file__).resolve().parent
    print(json.dumps(analyze(archive, args.dataset_root), indent=2))
