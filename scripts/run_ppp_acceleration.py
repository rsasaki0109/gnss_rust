"""Freeze native PVA prediction, 5 Hz car physics and full-rate PPC CV/PVA controls."""
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
import sys


def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def rows(data):return list(csv.DictReader(io.StringIO(data.decode())))


def acceleration_fields(rr,enabled):
    keys=['ax_mps2','ay_mps2','az_mps2','acceleration_qxx_m2_s4','acceleration_qyy_m2_s4','acceleration_qzz_m2_s4']
    for r in rr:
        if not enabled:assert all(k not in r for k in keys)
        elif r['status']=='PPP_FLOAT':assert all(math.isfinite(float(r[k])) for k in keys) and all(float(r[k])>0 for k in keys[3:])
        else:assert all(not r[k] for k in keys)


def car_metrics(data,truth,noisy,pva):
    rr=rows(data);assert len(rr)==len(truth)==64;acceleration_fields(rr,pva);errors=[];velocities=[];accelerations=[]
    assert {r['status'] for r in rr}<={'PPP_FLOAT','INVALID','EVENT'}
    for i,(r,t) in enumerate(zip(rr,truth)):
        assert (int(r['week']),Decimal(r['tow']).quantize(Decimal('.000001')))==(int(t['week']),Decimal(t['tow']).quantize(Decimal('.000001')))
        if r['status']!='PPP_FLOAT':assert all(not r[k] for k in ['x_m','y_m','z_m','vx_mps','vy_mps','vz_mps']);continue
        assert not r['error']
        assert all(math.isfinite(float(r[k])) for k in ['x_m','y_m','z_m','vx_mps','vy_mps','vz_mps','clock_s','zenith_m','position_covariance_trace_m2','velocity_covariance_trace_m2_s2','code_rms_m','phase_rms_m','nis'])
        errors.append((i,math.sqrt(sum((float(r[k+'_m'])-float(t[k]))**2 for k in ['x','y','z']))))
        velocities.append((i,math.sqrt(sum((float(r[k+'_mps'])-float(t[k]))**2 for k in ['vx','vy','vz']))))
        if pva:accelerations.append((i,math.sqrt(sum((float(r[k+'_mps2'])-float(t[k]))**2 for k in ['ax','ay','az']))))
    def summary(v):return dict(epochs=len(v),RMS=math.sqrt(sum(x*x for _,x in v)/len(v)),maximum=max(x for _,x in v)) if v else None
    bound=1.0 if noisy else .01;velocity_bound=1.0
    position_failures=[dict(epoch=i,error_m=x) for i,x in errors if x>=bound];velocity_failures=[dict(epoch=i,error_mps=x) for i,x in velocities if i>=8 and x>=velocity_bound]
    return dict(input_epochs=64,status_counts=dict(collections.Counter(r['status'] for r in rr)),position_3D_m=summary(errors),velocity_3D_mps=summary(velocities),velocity_after_eight_epochs_mps=summary([(i,x) for i,x in velocities if i>=8]),acceleration_3D_mps2=summary(accelerations),position_bound_m=bound,velocity_after_eight_bound_mps=velocity_bound,position_bound_failures=position_failures,velocity_bound_failures=velocity_failures,predeclared_targets_met=len(errors)==64 and not position_failures and not velocity_failures,terminal_error_reasons=dict(collections.Counter(r['error'] for r in rr if r['status']=='INVALID')))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--native-fixtures-dir',type=Path,required=True);a=p.parse_args();repo=Path(__file__).resolve().parents[1];out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True);assert not (out/'predeclared-inputs.json').exists()
    prior=repo/'docs/comparison/2026-10-21-ppp-rejected-arcs';old=json.loads((prior/'predeclared-inputs.json').read_text());cases=old['cases'];antex=next(f for f in old['input_sha256'] if f.endswith('.atx'))
    native=a.native_fixtures_dir.resolve();ng=json.loads((native/'native-generator.json').read_text());assert ng['returncode']==0 and len(ng['repeat_runs'])==2 and ng['repeat_runs'][0]['fixture_sha256']==ng['repeat_runs'][1]['fixture_sha256']
    assert all(sha(repo/f)==h for f,h in ng['source_sha256'].items());assert all(sha(repo/'tests/fixtures'/f)==h for f,h in ng['repeat_runs'][0]['fixture_sha256'].items())
    binaries={n:repo/'target/release/examples'/n for n in ['ppc_ppp','ppp_kinematic']}
    sources=list((repo/'src').rglob('*.rs'))+[repo/f for f in ['examples/ppc_ppp.rs','examples/ppp_kinematic.rs','tests/ppp_kinematic.rs','scripts/run_ppp_acceleration.py','scripts/build_ppp_acceleration_fixtures.py','scripts/upstream_ppp_acceleration.cpp','scripts/upstream_precise_transmit.cpp','scripts/upstream_precise.cpp','scripts/upstream_navigation.cpp','scripts/run_ppc_ppp.py','scripts/run_ppc_ppp_initial_prior.py','scripts/run_ppc_rtk_float.py','scripts/run_native_ppc_comparison.py','scripts/run_ppp_kinematic.py','scripts/measure_ppp_process.py','docs/ppp-acceleration.md','Cargo.toml','Cargo.lock','rust-toolchain.toml']]+list((repo/'tests/fixtures').glob('*ppp_acceleration*'))
    profiles={head+'-'+clock+'-'+motion:dict(head=head,clock=clock,motion=motion,arguments=([] if head=='default' else ['native-head'])+([] if clock=='random-walk' else ['-epoch-clock'])+([] if motion=='CV' else ['-pva'])) for head in ['default','native-head'] for clock in ['random-walk','independent'] for motion in ['CV','PVA']}
    manifest=dict(base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),source_sha256={str(f.relative_to(repo)):sha(f) for f in sources},binary_sha256={n:sha(f) for n,f in binaries.items()},input_sha256=old['input_sha256'],cases=cases,pairs=old['pairs'],profiles=profiles,prior_archive=str(prior.relative_to(repo)),prior_evidence_sha256=sha(prior/'evidence-sha256.json'),native_generator_sha256=sha(native/'native-generator.json'),full_raw_rover_epochs=True,thresholds_changed=False,default_position_model_changed=False,ambiguity_prior_changed=False,SPP_or_Doppler_fused=False,truth_used_by_solver=False,whole_native_PPP_executed=False,practical_parity_accepted=False,retention_enabled=False,native_head_velocity_variance_by_motion=dict(CV=25,PVA=100),synthetic_position_bound_m=dict(clean=.01,noisy=1.0),synthetic_velocity_after_eight_bound_mps=1.0)
    helper=runpy.run_path(str(repo/'scripts/run_native_ppc_comparison.py'));save,retain=[helper[k] for k in ['save','retain']];baseline=runpy.run_path(str(repo/'scripts/run_ppc_ppp.py'));initial=runpy.run_path(str(repo/'scripts/run_ppc_ppp_initial_prior.py'));scorer=runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))['score'];moving=runpy.run_path(str(repo/'scripts/run_ppp_kinematic.py'))
    assert all(sha(f)==h for f,h in manifest['input_sha256'].items());save(out/'predeclared-inputs.json',manifest)
    (out/'native-generator.json').write_bytes((native/'native-generator.json').read_bytes());(out/'native-build.log.gz').write_bytes((native/'build.log.gz').read_bytes())
    for i in range(2):
        for stream in ['stdout','stderr']:(out/f'native-run{i}.{stream}.log.gz').write_bytes((native/f'run{i}.{stream}.log.gz').read_bytes())
    report=dict(synthetic={},old_synthetic_controls={},cases={},comparisons={},practical_parity_accepted=False)
    def execute(stem,name,arguments,measure=False):
        assert all(sha(repo/f)==h for f,h in manifest['source_sha256'].items()) and sha(binaries[name])==manifest['binary_sha256'][name]
        command=[str(binaries[name])]+arguments;wrapper=[sys.executable,str(repo/'scripts/measure_ppp_process.py'),str(out/(stem+'.resource.json'))] if measure else []
        r=subprocess.run(wrapper+command,capture_output=True,cwd=repo);retain(out,stem,command,r.stdout,r.stderr,r.returncode);assert r.returncode==0,r.stderr.decode();return r.stdout
    truth=rows((repo/'tests/fixtures/synthetic_ppp_acceleration_truth.csv').read_bytes())
    for kind in ['car-clean','car-noisy','car-unrounded']:
        for motion in ['kinematic','pva']:
            stem='synthetic-'+kind+'-'+motion;data=execute(stem,'ppp_kinematic',[motion,kind]);metrics=car_metrics(data,truth,kind=='car-noisy',motion=='pva');report['synthetic'][stem]=metrics;save(out/'report.json',report)
            if motion=='pva':assert metrics['predeclared_targets_met'],stem
            print(stem,metrics['status_counts'],metrics['position_3D_m'],flush=True)
    for kind in ['clean','noisy']:
        stem='old-synthetic-'+kind;data=execute(stem,'ppp_kinematic',['kinematic',kind]);previous=repo/'docs/comparison/2026-10-21-ppp-epoch-clock'/('synthetic-'+kind+'-random-walk.stdout.csv.gz');assert data==gzip.decompress(previous.read_bytes());report['old_synthetic_controls'][stem]=dict(all_CSV_bytes_identical=True,previous=str(previous.relative_to(repo)),epochs=64);save(out/'report.json',report)
    for case,c in cases.items():
        with Path(c['reference']).open(newline='') as f:refs=[{k.strip():v for k,v in r.items()} for r in csv.DictReader(f)]
        references={(int(r['GPS Week']),Decimal(r['GPS TOW (s)'])):r for r in refs};assert len(references)==len(refs)
        for pair in manifest['pairs']:
            outputs={}
            for profile,s in profiles.items():
                stem=case+'-'+pair+'-'+profile;data=execute(stem,'ppc_ppp',[c[k] for k in ['rover','orbit','clock','bias']]+[antex,case,pair,'calibrated']+s['arguments'],True);rr,metrics=baseline['inspect'](data);acceleration_fields(rr,s['motion']=='PVA');outputs[profile]=rr;metrics['reference_errors']=scorer(initial['score_view'](rr),references,c['city']);metrics['resource_log']=stem+'.resource.json';metrics['runtime_memory_acceptance']=False
                if s['motion']=='CV':
                    previous=prior/(case+'-'+pair+'-'+s['head']+'-'+s['clock']+'-interrupt.stdout.csv.gz');assert data==gzip.decompress(previous.read_bytes()),stem;metrics['all_prior_CSV_bytes_identical']=True
                else:report['comparisons'][stem]=initial['compare'](rr,outputs[s['head']+'-'+s['clock']+'-CV'],references,c['city'],scorer)
                report['cases'].setdefault(case,{})[pair+'-'+profile]=metrics;save(out/'report.json',report);print(case,pair,profile,metrics['status_counts'],metrics['reference_errors']['hypotheses'].get('reference_is_antenna'),flush=True)
    assert all(sha(f)==h for f,h in manifest['input_sha256'].items());print('Completed 96 PPC PVA/CV profiles, 48 exact controls and 8 synthetic runs; practical parity unaccepted',flush=True)


if __name__=='__main__':main()
