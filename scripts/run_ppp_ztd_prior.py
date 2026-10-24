"""Freeze initial-ZTD uncertainty only: 36 PPC profiles and nine car traces."""
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
import subprocess
import sys


def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def prior_diagnostic(err,sigma,synthetic=False):
    pattern=rb'ZTD_initial_prior_sigma_m=(\S+) variance_m2=(\S+) initial_mean_m=(\S+) process_noise_m2_per_s=(\S+) ZTD_still_estimated=true other_settings_unchanged=true'
    match=re.search(pattern,err)
    if sigma is None:assert match is None;return None
    assert match;values=list(map(float,match.groups()));assert values[0]==sigma and values[1]==sigma*sigma and values[3]==1e-8 and all(math.isfinite(v) for v in values)
    if synthetic:assert values[2]==2.3
    return dict(sigma_m=values[0],variance_m2=values[1],mean_m=values[2],process_noise_m2_per_s=values[3],ZTD_still_estimated=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--native-fixtures-dir',type=Path,required=True);a=p.parse_args();repo=Path(__file__).resolve().parents[1];out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    prior=repo/'docs/comparison/2026-10-21-ppp-standardized-screening';old=json.loads((prior/'predeclared-inputs.json').read_text());cases=old['cases'];antex=next(f for f in old['input_sha256'] if f.endswith('.atx'));native=a.native_fixtures_dir.resolve();ng=json.loads((native/'native-generator.json').read_text());assert ng['returncode']==0 and len(ng['repeat_runs'])==2 and ng['repeat_runs'][0]['fixture_sha256']==ng['repeat_runs'][1]['fixture_sha256'];assert all(sha(repo/f)==h for f,h in ng['source_sha256'].items());assert all(sha(repo/'tests/fixtures'/f)==h for f,h in ng['repeat_runs'][0]['fixture_sha256'].items())
    binaries={n:repo/'target/release/examples'/n for n in ['ppc_ppp','ppp_kinematic']};profiles={name:dict(sigma_m=sigma,variance_m2=25.0 if sigma is None else sigma*sigma,arguments=[] if sigma is None else ['-ztd-sigma',str(sigma)]) for name,sigma in [('default',None),('sigma-0.6',0.6),('sigma-0.12',0.12)]}
    sources=list((repo/'src').rglob('*.rs'))+[repo/f for f in ['examples/ppc_ppp.rs','examples/ppp_kinematic.rs','examples/support/ppp_seed.rs','examples/support/ppp_ztd_prior.rs','scripts/run_ppp_ztd_prior.py','scripts/build_ppp_ztd_prior_fixtures.py','scripts/upstream_ppp_ztd_prior.cpp','scripts/run_ppp_acceleration.py','scripts/run_ppp_motion_seeds.py','scripts/run_ppc_ppp.py','scripts/run_ppc_ppp_initial_prior.py','scripts/run_ppc_rtk_float.py','scripts/run_native_ppc_comparison.py','scripts/measure_ppp_process.py','docs/ppp-ztd-prior.md','Cargo.toml','Cargo.lock','rust-toolchain.toml','tests/fixtures/upstream_ppp_ztd_prior.csv']]
    tracked=subprocess.check_output(['git','ls-tree','-r','--name-only','HEAD','src','tests','Cargo.toml','Cargo.lock','rust-toolchain.toml'],cwd=repo,text=True).splitlines();library={}
    for f in tracked:
        if f.endswith('.md'):continue
        assert subprocess.check_output(['git','show','HEAD:'+f],cwd=repo)==(repo/f).read_bytes(),f;library[f]=sha(repo/f)
    manifest=dict(base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),source_sha256={str(f.relative_to(repo)):sha(f) for f in sources},binary_sha256={n:sha(f) for n,f in binaries.items()},library_and_old_test_input_sha256=library,library_processing_unchanged=True,input_sha256=old['input_sha256'],cases=cases,pairs=old['pairs'],profiles=profiles,prior_archive=str(prior.relative_to(repo)),prior_evidence_sha256=sha(prior/'evidence-sha256.json'),native_generator_sha256=sha(native/'native-generator.json'),full_raw_rover_epochs=True,expected_PPC_profiles=36,expected_PPC_rows=351804,expected_exact_controls=12,expected_synthetic_profiles=9,synthetic_default_variance_m2=100.0,synthetic_mean_ZTD_m=2.3,synthetic_truth_ZTD_m=2.4,prior_mean_changed=False,process_or_measurement_noise_changed=False,existing_gates_changed=False,standardized_screening_enabled=False,retention_enabled=False,reference_or_IMU_used_by_solver=False,whole_native_PPP_executed=False,practical_parity_accepted=False,runtime_memory_acceptance=False)
    helper=runpy.run_path(str(repo/'scripts/run_native_ppc_comparison.py'));save,retain=[helper[k] for k in ['save','retain']];baseline=runpy.run_path(str(repo/'scripts/run_ppc_ppp.py'));initial=runpy.run_path(str(repo/'scripts/run_ppc_ppp_initial_prior.py'));seed=runpy.run_path(str(repo/'scripts/run_ppp_motion_seeds.py'));car=runpy.run_path(str(repo/'scripts/run_ppp_acceleration.py'));scorer=runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))['score'];assert all(sha(f)==h for f,h in manifest['input_sha256'].items());save(out/'predeclared-inputs.json',manifest)
    for name in ['native-generator.json','build.log.gz']+[f'run{i}.{stream}.log.gz' for i in range(2) for stream in ['stdout','stderr']]:(out/('native-'+name)).write_bytes((native/name).read_bytes())
    report=dict(synthetic={},invalid_sigma_CLI=[],cases={},comparisons={},practical_parity_accepted=False)
    def execute(stem,name,args,measure=False):
        assert all(sha(repo/f)==h for f,h in manifest['source_sha256'].items()) and sha(binaries[name])==manifest['binary_sha256'][name]
        command=[str(binaries[name])]+args;wrapper=[sys.executable,str(repo/'scripts/measure_ppp_process.py'),str(out/(stem+'.resource.json'))] if measure else [];r=subprocess.run(wrapper+command,capture_output=True,cwd=repo);retain(out,stem,command,r.stdout,r.stderr,r.returncode);return r
    for name in binaries:
        for bad in ['0','-1','NaN','inf','1e309','1e-300']:
            args=['pva','car-clean'] if name=='ppp_kinematic' else ['missing']*6+['1C2W','calibrated'];stem='invalid-'+name+'-'+bad.replace('-','minus');r=execute(stem,name,args+['-ztd-sigma',bad]);assert r.returncode!=0 and b'must be finite and positive' in r.stderr;report['invalid_sigma_CLI'].append(dict(stem=stem,returncode=r.returncode,solver_was_not_run=True))
    truth=car['rows']((repo/'tests/fixtures/synthetic_ppp_acceleration_truth.csv').read_bytes())
    for kind in ['car-clean','car-noisy','car-unrounded']:
        for profile,s in profiles.items():
            stem='synthetic-'+kind+'-'+profile;r=execute(stem,'ppp_kinematic',['pva',kind]+s['arguments']);assert r.returncode==0;metrics=car['car_metrics'](r.stdout,truth,kind=='car-noisy',True);metrics['prior_diagnostic']=prior_diagnostic(r.stderr,s['sigma_m'],True)
            if profile=='default':
                previous=repo/'docs/comparison/2026-10-21-ppp-acceleration'/('synthetic-'+kind+'-pva.stdout.csv.gz');assert r.stdout==gzip.decompress(previous.read_bytes()) and metrics['predeclared_targets_met'];metrics['all_prior_CSV_bytes_identical']=True
            report['synthetic'][stem]=metrics;save(out/'report.json',report);print(stem,metrics['status_counts'],metrics['position_3D_m'],'existing targets',metrics['predeclared_targets_met'],flush=True)
    total=controls=0
    for case,c in cases.items():
        with Path(c['reference']).open(newline='') as f:refs=[{k.strip():v for k,v in r.items()} for r in csv.DictReader(f)]
        references={(int(r['GPS Week']),Decimal(r['GPS TOW (s)'])):r for r in refs};assert len(references)==len(refs);seed_control=None;mean=None
        for pair in manifest['pairs']:
            outputs={}
            for profile,s in profiles.items():
                stem=case+'-'+pair+'-'+profile;args=[c[k] for k in ['rover','orbit','clock','bias']]+[antex,case,pair,'calibrated','native-head','-epoch-clock','-pva','-seed-doppler',c['navigation']]+s['arguments'];r=execute(stem,'ppc_ppp',args,True);assert r.returncode==0,r.stderr.decode();rr,metrics=baseline['inspect'](r.stdout);total+=len(rr);outputs[profile]=rr;assert all(len(r)==72 for r in rr);metrics['prior_diagnostic']=prior_diagnostic(r.stderr,s['sigma_m']);metrics['reference_errors']=scorer(initial['score_view'](rr),references,c['city']);metrics['seeds_and_service']=seed['seed_metrics'](rr,references,c['city'],scorer);metrics['resource_log']=stem+'.resource.json';metrics['runtime_memory_acceptance']=False
                cfg=re.search(rb'initial_head_variance_m2: \[([^\]]+)\]',r.stderr);assert cfg;head=list(map(float,cfg[1].split(b',')));assert head==[100.,100.,100.,1e8,s['variance_m2']];observed_mean=float(re.search(rb'initial_zenith_delay_m: ([^,}]+)',r.stderr)[1]);mean=observed_mean if mean is None else mean;assert observed_mean==mean;metrics['initial_ZTD_mean_m']=mean
                keys=[k for k in rr[0] if k.startswith('seed_')][:23];ss=[[r[k] for k in keys] for r in rr]
                if seed_control is None:seed_control=ss
                else:assert ss==seed_control
                metrics['all_seed_source_fields_identical']=True
                if profile=='default':assert r.stdout==gzip.decompress((prior/(case+'-'+pair+'-off.stdout.csv.gz')).read_bytes());metrics['all_72_prior_CSV_bytes_identical']=True;controls+=1
                else:report['comparisons'][stem]=initial['compare'](rr,outputs['default'],references,c['city'],scorer)
                report['cases'].setdefault(case,{})[pair+'-'+profile]=metrics;save(out/'report.json',report);print(case,pair,profile,metrics['status_counts'],metrics['reference_errors']['hypotheses'].get('reference_is_antenna'),flush=True)
    assert total==351804 and controls==12;assert all(sha(f)==h for f,h in manifest['input_sha256'].items());print('Completed 36 PPC / 351804 rows, 12 exact controls, nine car traces and 12 invalid sigma checks; parity unaccepted',flush=True)


if __name__=='__main__':main()
