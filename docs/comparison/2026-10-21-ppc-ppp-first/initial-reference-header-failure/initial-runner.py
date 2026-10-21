"""First full 5 Hz PPC moving PPP diagnostic; all statuses and both origins retained."""
import argparse
import collections
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import runpy
import subprocess
import time


def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def rows(data):return list(csv.DictReader(io.StringIO(data.decode())))
def inspect(data):
    rr=rows(data);assert rr
    statuses=collections.Counter(r['status'] for r in rr);assert set(statuses)<={'PPP_FLOAT','INVALID','EVENT'}
    previous=None;longest=streak=0;last_success=None;max_gap=0.0;uncommitted=0
    for r in rr:
        t=Decimal(r['week'])*604800+Decimal(r['tow']);assert previous is None or t>previous;previous=t
        if r['status']=='PPP_FLOAT':
            assert all(math.isfinite(float(r[k])) for k in ['x_m','y_m','z_m','vx_mps','vy_mps','vz_mps','nis','code_rms_m','phase_rms_m','clock_s','zenith_m','position_covariance_trace_m2','velocity_covariance_trace_m2_s2'])
            assert int(r['satellites'])>=5 and int(r['ambiguities'])>=4 and not r['error']
            assert (r['committed_week'],r['committed_tow'])==(r['week'],r['tow'])
            streak+=1;longest=max(longest,streak)
            if last_success is not None:max_gap=max(max_gap,float(t-last_success))
            last_success=t
        else:
            assert all(not r[k] for k in ['x_m','y_m','z_m','vx_mps','vy_mps','vz_mps']);streak=0
            uncommitted+=not r['committed_week']
    return rr,dict(input_epochs=len(rr),status_counts=dict(statuses),successful_fraction=statuses['PPP_FLOAT']/len(rr),longest_successful_epoch_streak=longest,maximum_success_to_success_gap_s=max_gap if statuses['PPP_FLOAT']>1 else None,epochs_before_any_commit=uncommitted,terminal_error_reasons=dict(collections.Counter(r['error'] for r in rr if r['status']=='INVALID')),initial_success_time=next(([r['week'],r['tow']] for r in rr if r['status']=='PPP_FLOAT'),None))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--dataset-root',type=Path,required=True);p.add_argument('--product-root',type=Path,required=True);p.add_argument('--antex',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
    repo=Path(__file__).resolve().parents[1];out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True);assert not (out/'predeclared-inputs.json').exists()
    binary=repo/'target/release/examples/ppc_ppp';helper=runpy.run_path(str(repo/'scripts/run_native_ppc_comparison.py'));save,retain=[helper[k] for k in ['save','retain']]
    scorer=runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))['score']
    profiles={pair+'-'+phase:dict(pair=pair,phase=phase,code_sigma_per_frequency_m=1.0,phase_sigma_per_frequency_m=.01,elevation_weighting=False,model='ConstantVelocity default',windup='NominalYawApproximateSun',satellite_IF_PCO=True,receiver_calibration=False,satellite_PCV=False,solid_earth_tide=False,thresholds='unchanged PppFilterSettings default') for pair in ['1C2W','1C2L'] for phase in ['calibrated','raw-float']}
    cases={}
    for city in ['nagoya','tokyo']:
        for run in ['run1','run2','run3']:
            stamp='20242050000' if city=='tokyo' else '20242020000' if run=='run2' else '20242160000';root=a.dataset_root/city/run
            cases[city+'-'+run]=dict(city=city,rover=str((root/'rover.obs').resolve()),reference=str((root/'reference.csv').resolve()),orbit=str((a.product_root/f'COD0OPSRAP_{stamp}_01D_05M_ORB.SP3').resolve()),clock=str((a.product_root/f'COD0OPSRAP_{stamp}_01D_30S_CLK.CLK').resolve()),bias=str((a.product_root/f'COD0OPSRAP_{stamp}_01D_01D_OSB.BIA').resolve()))
    source=list((repo/'src').rglob('*.rs'))+[repo/f for f in ['examples/ppc_ppp.rs','scripts/run_ppc_ppp.py','scripts/run_ppc_rtk_float.py','scripts/run_native_ppc_comparison.py','Cargo.toml','Cargo.lock','rust-toolchain.toml']]
    inputs=sorted({Path(c[k]) for c in cases.values() for k in ['rover','reference','orbit','clock','bias']}|{a.antex.resolve()})
    manifest=dict(base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),source_sha256={str(p.relative_to(repo)):sha(p) for p in source},binary_sha256=sha(binary),input_sha256={str(p):sha(p) for p in inputs},cases=cases,profiles=profiles,full_raw_rover_epochs=True,base_observations_used=False,SPP_or_Doppler_fused=False,truth_used_by_solver=False,reference_used_only_after_solver=True,whole_native_PPP_executed=False,runtime_memory_acceptance=False,practical_parity_accepted=False)
    save(out/'predeclared-inputs.json',manifest);report=dict(cases={},common_success_comparisons={},whole_native_PPP_executed=False,practical_parity_accepted=False)
    for name,c in cases.items():
        profiles_rows={}
        with Path(c['reference']).open(newline='') as stream:
            refs=[{k.strip():v for k,v in r.items()} for r in csv.DictReader(stream)]
        references={(int(r['GPSWeek']),Decimal(r['GPSTOW (s)'])):r for r in refs}
        for profile,settings in profiles.items():
            assert all(sha(repo/f)==h for f,h in manifest['source_sha256'].items()) and sha(binary)==manifest['binary_sha256']
            command=[str(binary)]+[c[k] for k in ['rover','orbit','clock','bias']]+[str(a.antex.resolve()),name,settings['pair'],settings['phase']]
            stem=name+'-'+profile;time_file=out/(stem+'.resource.log');started=time.perf_counter();r=subprocess.run(['/usr/bin/time','-v','-o',str(time_file)]+command,capture_output=True,cwd=repo);elapsed=time.perf_counter()-started
            retain(out,stem,command,r.stdout,r.stderr,r.returncode);assert r.returncode==0,r.stderr.decode()
            rr,metrics=inspect(r.stdout);profiles_rows[profile]=rr
            # Read-only view for the shared scorer; retained solver status remains PPP_FLOAT.
            score_rows=[dict(row,status='FLOAT') if row['status']=='PPP_FLOAT' else row for row in rr]
            metrics['reference_errors']=scorer(score_rows,references,c['city']);metrics['wall_seconds_including_load_and_CSV']=elapsed;metrics['resource_log']=time_file.name;metrics['runtime_memory_acceptance']=False
            report['cases'].setdefault(name,{})[profile]=metrics;save(out/'report.json',report);print(name,profile,metrics['status_counts'],metrics['terminal_error_reasons'],flush=True)
        for pair in ['1C2W','1C2L']:
            aa,bb=[{(r['week'],r['tow']):r for r in profiles_rows[pair+'-'+phase] if r['status']=='PPP_FLOAT'} for phase in ['calibrated','raw-float']];common=aa.keys()&bb.keys();diff=[math.sqrt(sum((float(aa[t][k])-float(bb[t][k]))**2 for k in ['x_m','y_m','z_m'])) for t in common]
            report['common_success_comparisons'][name+'-'+pair]=dict(common_epochs=len(common),calibrated_only=len(aa.keys()-bb.keys()),raw_only=len(bb.keys()-aa.keys()),position_difference_3D_RMS_m=math.sqrt(sum(v*v for v in diff)/len(diff)) if diff else None,position_difference_maximum_m=max(diff) if diff else None)
        save(out/'report.json',report)
    assert all(sha(Path(f))==h for f,h in manifest['input_sha256'].items())
    print('Completed all 24 profiles; no whole-native PPP or practical parity acceptance',flush=True)


if __name__=='__main__':main()
