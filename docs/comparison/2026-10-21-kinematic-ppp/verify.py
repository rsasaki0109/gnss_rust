"""Verify native PPP prediction, moving synthetic traces and exact static controls."""
import gzip
import hashlib
import json
from pathlib import Path
import re
import runpy


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return gzip.decompress(p.read_bytes())
def main():
    archive=Path(__file__).resolve().parent;repo=archive.parents[2]
    m=json.loads((archive/"predeclared-inputs.json").read_text());report=json.loads((archive/"report.json").read_text())
    assert not any(m[k] for k in ["default_position_model_changed","whole_native_PPP_executed","real_PPC_PPP_executed","truth_used_by_solver","existing_thresholds_changed","practical_parity_accepted"])
    assert not any(report[k] for k in ["whole_native_PPP_executed","real_PPC_PPP_executed","practical_parity_accepted"])
    assert m["native_predict_component_calls"]==24 and m["original_full_trace_clean_target_m"]==.002
    assert len(m["moving_profiles"])==len(report["moving"])==11 and len(m["static_profiles"])==len(report["static_controls"])==10
    for path,h in m["source_sha256"].items():assert sha(repo/path)==h,path
    helper=runpy.run_path(str(repo/"scripts/run_ppp_kinematic.py"));truth=helper["rows"]((repo/"tests/fixtures/synthetic_ppp_kinematic_truth.csv").read_bytes())
    def capture(stem):
        e=json.loads((archive/(stem+".execution.json")).read_text());data=read(archive/(stem+".stdout.csv.gz"));stderr=read(archive/(stem+".stderr.log.gz"))
        assert e["returncode"]==0 and hashlib.sha256(data).hexdigest()==e["csv_sha256"] and hashlib.sha256(stderr).hexdigest()==e["stderr_sha256"]
        return data,stderr
    for name,c in m["moving_profiles"].items():
        data,stderr=capture(name);assert b"truth_used_by_solver=false" in stderr and b"whole_native_PPP_executed=false" in stderr
        assert helper["moving_metrics"](data,truth,c["moving"])==report["moving"][name]
    assert not report["moving"]["kinematic_clean"]["full_trace_2mm_target_met"]
    assert [x["epoch"] for x in report["moving"]["kinematic_clean"]["two_mm_exceedances"]]==[1,2,3]
    assert report["moving"]["kinematic_unrounded_phase"]["full_trace_2mm_target_met"] and report["moving"]["kinematic_unrounded_both"]["full_trace_2mm_target_met"] and not report["moving"]["kinematic_unrounded_code"]["full_trace_2mm_target_met"]
    for name in m["static_profiles"]:
        old,_=capture("static-old-"+name);new,_=capture("static-new-"+name);assert old==new
        assert report["static_controls"][name]==dict(epochs=64,all_CSV_bytes_identical=True,csv_sha256=hashlib.sha256(new).hexdigest())
        rows=helper["rows"](new);assert len(rows)==64 and all(r["status"]=="PPP_FLOAT" for r in rows)
    old=json.loads((archive/"old-static-build.json").read_text());assert sha(archive/"old-static-build.json")==m["old_static_build_sha256"] and old["base_commit"]==m["base_commit"] and old["builder_sha256"]==sha(repo/old.get("builder_source","scripts/build_ppp_static_control.py"))
    for path,h in old["source_sha256"].items():
        actual=hashlib.sha256(read(archive/"old-static-ppp.rs.gz")).hexdigest() if path=="src/ppp.rs" else sha(repo/path)
        assert actual==h,path
    for c in old["checks"]:assert c["returncode"]==0 and hashlib.sha256(read(archive/("old-static-"+c["log"]))).hexdigest()==c["uncompressed_sha256"]
    replay=json.loads((archive/"replay-verification.json").read_text());assert replay["builder_sha256"]==sha(archive/"replay_static_control.py") and set(replay["cases"])==set(m["static_profiles"])
    for name,c in replay["cases"].items():assert c["returncode"]==0 and c["same_as_archived_static_control"] and c["stdout_sha256"]==report["static_controls"][name]["csv_sha256"]
    rb=json.loads((archive/"replay-static-build.json").read_text());assert rb["source_sha256"]==old["source_sha256"] and rb["original_build_sha256"]==sha(archive/"old-static-build.json") and rb["builder_sha256"]==replay["builder_sha256"]
    for c in rb["checks"]:assert c["returncode"]==0 and hashlib.sha256(read(archive/("replay-static-"+c["log"]))).hexdigest()==c["uncompressed_sha256"]
    native=json.loads((archive/"native-generator.json").read_text());assert native["returncode"]==0 and native["actual_native_PPP_predictState_calls_per_repeat"]==24 and native["unrounded_observable_rows"]==576 and not native["whole_native_PPP_executed"] and not native["upstream_source_patched"]
    assert native["build_log_sha256"]==sha(archive/"native-generator-build.log.gz")
    for path,h in native["source_sha256"].items():assert sha(repo/path)==h,path
    for path,h in native["earlier_native_source_proof"].items():assert sha(repo/path)==h,path
    baseline=json.loads((repo/"docs/comparison/2026-10-21-gps-galileo-float/native-comparison/predeclared-inputs.json").read_text())
    assert all(baseline["native_archive_sha256"][path]==h for path,h in native["native_library_sha256"].items())
    expected={f.name:sha(f) for f in (repo/"tests/fixtures").glob("*ppp_kinematic*")};assert len(expected)==5
    assert len(native["repeat_runs"])==2 and all(r["returncode"]==0 and r["fixture_sha256"]==expected for r in native["repeat_runs"])
    prior=json.loads((archive/"pre-unrounded-fixture-sha256.json").read_text());assert all(sha(repo/f)==h for f,h in prior.items())
    q=json.loads((archive/"qualification.json").read_text());assert q["predeclared_inputs_sha256"]==sha(archive/"predeclared-inputs.json") and {c["name"] for c in q["checks"]}=={"rust","examples","fmt","clippy","release"}
    for c in q["checks"]:
        data=read(archive/c["log"]);assert c["returncode"]==0 and hashlib.sha256(data).hexdigest()==c["uncompressed_sha256"]
        if c["name"]=="rust":assert sum(map(int,re.findall(rb"test result: ok\. (\d+) passed;",data)))==c["passed_tests"]==398
    probe=json.loads((archive/"probes/probes.json").read_text());assert not probe["QR_retained_in_solver"]
    assert probe["QR_function_sha256"]==sha(archive/"probes/unused-whitened-qr.rs")
    for item in probe["files"].values():assert hashlib.sha256(read(archive/"probes"/item["archive"])).hexdigest()==item["stdout_sha256"]
    sealed=json.loads((archive/"evidence-sha256.json").read_text());assert set(sealed)=={str(p.relative_to(archive)) for p in archive.rglob("*") if p.is_file() and p.name!="evidence-sha256.json"}
    for path,h in sealed.items():assert sha(archive/path)==h,path
    print("Verified 24 actual native PPP predictions, 11 moving profiles, 10 byte-identical static controls and 398 Rust tests; real PPC PPP untested")


if __name__=="__main__":main()
