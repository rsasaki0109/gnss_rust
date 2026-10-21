"""Run whole pinned native RTK on decoded-equal GPS/Galileo E1 PPC inputs."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
from itertools import zip_longest
import json
from pathlib import Path
import runpy
import subprocess
import tempfile


def save(path,value):path.write_text(json.dumps(value,indent=2)+"\n")


def variant(original):
    start=original.index("static ObservationData gps1c(")
    end=original.index("static const char* status_name",start)
    helpers=r'''
static bool eligible(const std::pair<SatelliteId,std::string>& key) {
    return (key.first.system==GNSSSystem::GPS && key.second=="1C") ||
        (key.first.system==GNSSSystem::Galileo && (key.second=="1C" || key.second=="1X"));
}
static ObservationData selected_l1(const ObservationData& raw) {
    ObservationData selected(raw.time);
    selected.receiver_position=raw.receiver_position;selected.receiver_clock_bias=raw.receiver_clock_bias;
    std::map<SatelliteId,std::pair<std::string,Observation>> choices;
    auto rank=[](const std::string& code,const Observation& o) {
        bool complete=o.has_pseudorange && std::isfinite(o.pseudorange) && o.pseudorange>0 &&
            o.has_carrier_phase && std::isfinite(o.carrier_phase);
        return std::pair<bool,bool>(!complete,code!="1C");
    };
    for(const auto&[key,o]:raw.rinex_tracking_observations) {
        if(!eligible(key))continue;
        const auto old=choices.find(key.first);
        if(old==choices.end() || rank(key.second,o)<rank(old->second.first,old->second.second))
            choices[key.first]={key.second,o};
    }
    for(const auto&[sat,choice]:choices) {
        (void)sat;selected.observations.push_back(choice.second);
        selected.addRinexTrackingObservation(choice.first,choice.second);
    }
    selected.setRinexFrequencySlot(GNSSSystem::GPS,0,"1C");
    return selected;
}
static bool continuity(const ObservationData& raw) {
    for(const auto&[key,value]:raw.rinex_tracking_observations)
        if(eligible(key) && (value.lli&3))return true;
    return false;
}
static void ledger(std::ostream& out,const char* side,const ObservationData& epoch,bool reset) {
    for(const auto&[key,o]:epoch.rinex_tracking_observations) {
        if(!o.has_pseudorange && !o.has_carrier_phase && !o.has_doppler && o.snr==0.0 && o.lli==0)continue;
        out<<side<<','<<epoch.time.week<<','<<epoch.time.tow<<','<<o.satellite.toString()<<','<<key.second<<','
           <<o.has_pseudorange<<','<<o.pseudorange<<','<<o.has_carrier_phase<<','<<o.carrier_phase<<','
           <<o.has_doppler<<','<<o.doppler<<','<<o.snr<<','<<static_cast<int>(o.lli)<<','<<reset<<'\n';
    }
}
'''
    source=original[:start]+helpers+original[end:]
    replacements={
        "if (argc != 7 && argc != 8)":"if (argc != 8)",
        "[SMOKE_MAX_PAIRS]":"previous|spp",
        "const std::size_t maximum = argc == 8 ? std::stoull(argv[7]) : 0;":
            'const std::string seed_policy=argv[7]; if(seed_policy!="previous" && seed_policy!="spp") throw std::runtime_error("seed policy"); const std::size_t maximum=0;',
        "config.prefer_rover_position_seed = true;":'config.prefer_rover_position_seed = seed_policy=="previous";',
        "side,week,tow,satellite,has_code":"side,week,tow,satellite,tracking_code,has_code",
        "gps1c(rover)":"selected_l1(rover)","gps1c(base)":"selected_l1(base)",
        "GPS_1C_only=true":"GPS_1C_Galileo_E1_only=true",
        '<< " previous_own_RTK_seed=true AR_default_continuous=true performance_accepted=false\\n";':
            '<< " seed_policy=" << seed_policy << " AR_default_continuous=true performance_accepted=false\\n";',
    }
    for old,new in replacements.items():
        assert source.count(old)==1,(old,source.count(old));source=source.replace(old,new)
    return source


def ledger_key(row,helper):
    key=helper(row)
    return key[:4]+[row["tracking_code"]]+key[4:]


def compare(native,rust,helper):
    digest=hashlib.sha256();count=0
    with native.open() as a,rust.open() as b:
        for count,(left,right) in enumerate(zip_longest(csv.DictReader(a),csv.DictReader(b)),1):
            assert left is not None and right is not None,(count,"length")
            key=ledger_key(left,helper)
            assert key==ledger_key(right,helper),(count,left,right)
            digest.update((json.dumps(key,separators=(",",":"))+"\n").encode())
    assert count>0
    return dict(rows=count,normalized_sha256=digest.hexdigest(),all_decoded_values_tracking_and_resets_equal=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root",type=Path,required=True)
    parser.add_argument("--native-build",type=Path,required=True)
    parser.add_argument("--rust-ledger-binary",type=Path,required=True)
    parser.add_argument("--rust-archive",type=Path,required=True)
    args=parser.parse_args();repo=Path(__file__).resolve().parents[1]
    parent=args.rust_archive.resolve();output=parent/"native-comparison";output.mkdir(exist_ok=False)
    primary=json.loads((parent/"predeclared-inputs.json").read_text())
    native_base=repo/"docs/comparison/2026-10-21-native-rtk-common-inputs"
    provenance=json.loads((native_base/"native-build.json").read_text())
    control=json.loads((native_base/"predeclared-inputs.json").read_text())
    runner=runpy.run_path(str(repo/"scripts/run_native_ppc_comparison.py"));sha=runner["sha"]
    scoring=runpy.run_path(str(repo/"scripts/run_ppc_rtk_float.py"))
    for name,h in primary["source_sha256"].items():assert sha(repo/name)==h,name
    for p,h in control["native_archive_sha256"].items():assert sha(Path(p))==h,p
    source=output/"native_ppc_mixed.cpp";source.write_text(variant((repo/"scripts/native_ppc_rtk.cpp").read_text()))
    build=args.native_build.resolve();binary=build/"native_ppc_mixed";rust=args.rust_ledger_binary.resolve();dataset=args.dataset_root.resolve()
    command=provenance["commands"][-1][:];command[command.index("scripts/native_ppc_rtk.cpp")]=str(source);command[-1]=str(binary)
    result=subprocess.run(command,capture_output=True,cwd=repo)
    (output/"build.stderr.log.gz").write_bytes(gzip.compress(result.stderr,mtime=0))
    save(output/"build.json",dict(command=command,returncode=result.returncode,stderr_sha256=hashlib.sha256(result.stderr).hexdigest()))
    assert result.returncode==0,result.stderr.decode()
    names=["scripts/run_native_galileo.py","scripts/native_ppc_rtk.cpp","examples/ppc_mixed_input_ledger.rs","scripts/run_native_ppc_comparison.py","scripts/run_ppc_rtk_float.py"]
    manifest=dict(source_sha256={n:sha(repo/n) for n in names},adapter_sha256=sha(source),
        binary_sha256={str(p):sha(p) for p in [binary,rust]},native_archive_sha256=control["native_archive_sha256"],
        parent_manifest_sha256=sha(parent/"predeclared-inputs.json"),parent_report_sha256=sha(parent/"report.json"),
        cases=primary["cases"],input_sha256=primary["input_sha256"],seed_policies=["previous","spp"],
        tracking="complete GPS 1C/Galileo 1C or 1X, prefer 1C per receiver",native_AR_QC_defaults_preserved=True,
        SPP="native default may use Galileo; Rust SPP seed remains GPS-only",truth_used_by_solver=False,
        timing_memory_acceptance=False,practical_PPC_parity_accepted=False)
    save(output/"predeclared-inputs.json",manifest);report=dict(practical_PPC_parity_accepted=False,cases={})
    with tempfile.TemporaryDirectory(prefix="native-mixed-ppc-") as temporary:
        temporary=Path(temporary)
        for case in primary["cases"]:
            city,run=case.split("-");directory=dataset/city/run;report["cases"][case]={}
            references={}
            with (directory/"reference.csv").open() as stream:
                for raw in csv.DictReader(stream):
                    r={k.strip():v.strip() for k,v in raw.items()};references[(int(r["GPS Week"]),Decimal(r["GPS TOW (s)"]))]=r
            rust_inputs=temporary/(case+".rust.csv")
            with rust_inputs.open("wb") as stream:decoded=subprocess.run([str(rust),str(directory)],stdout=stream,stderr=subprocess.PIPE,cwd=repo)
            assert decoded.returncode==0,decoded.stderr.decode();first_input_hash=None
            for policy in manifest["seed_policies"]:
                native_csv=temporary/(case+"."+policy+".csv");native_inputs=temporary/(case+"."+policy+".inputs.csv")
                command=[str(binary),str(directory),*map(str,primary["base_ecef_m"][city]),str(native_csv),str(native_inputs),policy]
                result=subprocess.run(command,capture_output=True,cwd=repo);data=native_csv.read_bytes() if native_csv.exists() else b""
                runner["retain"](output,case+"."+policy,command,data,result.stderr,result.returncode,
                                 csv_capture="explicit native output file",process_stdout=result.stdout)
                assert result.returncode==0,result.stderr.decode()
                try:comparison=compare(native_inputs,rust_inputs,runner["ledger_key"])
                except AssertionError:
                    for name,path in [("native",native_inputs),("rust",rust_inputs)]:
                        (output/(case+"."+policy+".FAILED."+name+".ledger.csv.gz")).write_bytes(gzip.compress(path.read_bytes(),mtime=0))
                    raise
                comparison.update(native_ledger_sha256=sha(native_inputs),rust_ledger_sha256=sha(rust_inputs))
                assert first_input_hash is None or first_input_hash==sha(native_inputs)
                first_input_hash=sha(native_inputs)
                if policy=="previous":
                    for name,path in [("native",native_inputs),("rust",rust_inputs)]:
                        (output/(case+"."+name+".input-ledger.csv.gz")).write_bytes(gzip.compress(path.read_bytes(),mtime=0))
                rows,rtk,metrics=runner["native_summary"](data)
                assert metrics["valid_outputs_with_mismatched_time"]==0 and b"input_drained=1" in result.stderr
                with gzip.open(parent/(case+".galileo.stdout.csv.gz"),"rt") as stream:rust_rows=list(csv.DictReader(stream))
                keys=lambda rows:[(r["kind"],int(r["week"]),float(r["tow"]).hex(),r["side"]) for r in rows]
                assert keys(rows)==keys(rust_rows)
                metrics["reference_scores"]=scoring["score"](rtk,references,city)
                metrics["FIX_reference_scores"]=scoring["score"]([{**r,"status":"FLOAT"} for r in rows if r["status"]=="FIX"],references,city)
                report["cases"][case][policy]=dict(native=metrics,input_comparison=comparison)
                save(output/"report.json",report);print(case,policy,metrics["status_counts"],"input rows",comparison["rows"],flush=True)
    assert all(sha(repo/p)==h for p,h in manifest["source_sha256"].items())
    assert all(sha(Path(p))==h for p,h in manifest["input_sha256"].items())


if __name__=="__main__":main()
