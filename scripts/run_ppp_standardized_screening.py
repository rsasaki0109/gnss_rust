"""Freeze 24 full-rate PPC controls/native-sigma-4 PPP screening comparisons."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
from pathlib import Path
import runpy
import subprocess
import sys


def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def legacy_columns(data,fields):
    rr=list(csv.DictReader(io.StringIO(data.decode())));out=io.StringIO();w=csv.DictWriter(out,fieldnames=fields,extrasaction='ignore',lineterminator='\n');w.writeheader()
    for r in rr:
        if 'seed_first_system_clock_s' in r:r['seed_gps_clock_s']=r.pop('seed_first_system_clock_s')
        if ';SPP_first_clock_system=' in r['seed_source']:
            r['seed_source'],clock=r['seed_source'].rsplit(';SPP_first_clock_system=',1);assert clock=='Gps'
        w.writerow(r)
    return out.getvalue().encode()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--native-fixtures-dir',type=Path,required=True);a=p.parse_args();repo=Path(__file__).resolve().parents[1];out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    prior=repo/'docs/comparison/2026-10-21-ppp-motion-seeds';old=json.loads((prior/'predeclared-inputs.json').read_text());cases=old['cases'];antex=next(f for f in old['input_sha256'] if f.endswith('.atx'));binary=repo/'target/release/examples/ppc_ppp';native=a.native_fixtures_dir.resolve();ng=json.loads((native/'native-generator.json').read_text())
    assert ng['returncode']==0 and len(ng['repeat_runs'])==2 and ng['repeat_runs'][0]['fixture_sha256']==ng['repeat_runs'][1]['fixture_sha256'];assert all(sha(repo/f)==h for f,h in ng['source_sha256'].items());assert all(sha(repo/'tests/fixtures'/f)==h for f,h in ng['repeat_runs'][0]['fixture_sha256'].items())
    sources=list((repo/'src').rglob('*.rs'))+[repo/f for f in ['examples/ppc_ppp.rs','examples/support/ppp_seed.rs','tests/ppp_kinematic.rs','tests/ppp_screening.rs','tests/phase_arc.rs','scripts/run_ppp_standardized_screening.py','scripts/build_ppp_screening_fixtures.py','scripts/upstream_ppp_screening.cpp','scripts/run_ppp_motion_seeds.py','scripts/run_ppc_ppp.py','scripts/run_ppc_ppp_initial_prior.py','scripts/run_ppc_rtk_float.py','scripts/run_native_ppc_comparison.py','scripts/measure_ppp_process.py','docs/ppp-standardized-screening.md','Cargo.toml','Cargo.lock','rust-toolchain.toml']]+list((repo/'tests/fixtures').glob('*ppp_screening*'))
    profiles={n:dict(arguments=[] if n=='off' else ['-native-screen'],threshold_sigma=None if n=='off' else 4.0) for n in ['off','screen']};manifest=dict(base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),source_sha256={str(f.relative_to(repo)):sha(f) for f in sources},binary_sha256=sha(binary),input_sha256=old['input_sha256'],cases=cases,pairs=old['pairs'],profiles=profiles,prior_archive=str(prior.relative_to(repo)),prior_evidence_sha256=sha(prior/'evidence-sha256.json'),native_generator_sha256=sha(native/'native-generator.json'),full_raw_rover_epochs=True,expected_profiles=24,expected_rows=234536,expected_controls=12,seed_mode='GPS/Galileo native-variance-QC SPP/Doppler',head='native-head',clock='IndependentEpoch',motion='PVA',phase='CODE calibrated',retention_enabled=False,other_GNSS_gates_or_noise_changed=False,reference_or_IMU_used_by_solver=False,whole_native_PPP_executed=False,practical_parity_accepted=False,runtime_memory_acceptance=False)
    helper=runpy.run_path(str(repo/'scripts/run_native_ppc_comparison.py'));save,retain=[helper[k] for k in ['save','retain']];baseline=runpy.run_path(str(repo/'scripts/run_ppc_ppp.py'));initial=runpy.run_path(str(repo/'scripts/run_ppc_ppp_initial_prior.py'));seed_helper=runpy.run_path(str(repo/'scripts/run_ppp_motion_seeds.py'));scorer=runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))['score'];assert all(sha(f)==h for f,h in manifest['input_sha256'].items());save(out/'predeclared-inputs.json',manifest)
    for name in ['native-generator.json','build.log.gz']+[f'run{i}.{stream}.log.gz' for i in range(2) for stream in ['stdout','stderr']]:(out/('native-'+name)).write_bytes((native/name).read_bytes())
    report=dict(cases={},comparisons={},practical_parity_accepted=False);total=controls=0
    for case,c in cases.items():
        with Path(c['reference']).open(newline='') as f:refs=[{k.strip():v for k,v in r.items()} for r in csv.DictReader(f)]
        references={(int(r['GPS Week']),Decimal(r['GPS TOW (s)'])):r for r in refs};assert len(references)==len(refs);seed_control=None
        for pair in manifest['pairs']:
            outputs={}
            for profile,s in profiles.items():
                assert all(sha(repo/f)==h for f,h in manifest['source_sha256'].items()) and sha(binary)==manifest['binary_sha256']
                stem=case+'-'+pair+'-'+profile;command=[str(binary)]+[c[k] for k in ['rover','orbit','clock','bias']]+[antex,case,pair,'calibrated','native-head','-epoch-clock','-pva','-seed-doppler',c['navigation']]+s['arguments'];wrapper=[sys.executable,str(repo/'scripts/measure_ppp_process.py'),str(out/(stem+'.resource.json'))];r=subprocess.run(wrapper+command,capture_output=True,cwd=repo);retain(out,stem,command,r.stdout,r.stderr,r.returncode);assert r.returncode==0,r.stderr.decode()
                rr,metrics=baseline['inspect'](r.stdout);total+=len(rr);outputs[profile]=rr;metrics['reference_errors']=scorer(initial['score_view'](rr),references,c['city']);metrics['seeds_and_service']=seed_helper['seed_metrics'](rr,references,c['city'],scorer);metrics['resource_log']=stem+'.resource.json';metrics['runtime_memory_acceptance']=False
                keys=[k for k in rr[0] if k.startswith('seed_')][:23];seeds=[[r[k] for k in keys] for r in rr]
                if seed_control is None:seed_control=seeds
                else:assert seeds==seed_control
                metrics['all_seed_source_fields_identical']=True
                if profile=='off':
                    data=gzip.decompress((prior/(case+'-'+pair+'-doppler.stdout.csv.gz')).read_bytes());fields=next(csv.reader(io.StringIO(data.decode())));assert len(fields)==72 and legacy_columns(r.stdout,fields)==data,stem;metrics['all_72_prior_columns_identical_after_audited_label_normalization']=True;controls+=1
                else:
                    assert all(len(r)==75 for r in rr);metrics['screen_attempts']=sum(int(r['screen_attempts']) for r in rr);metrics['epochs_with_screen_attempts']=sum(int(r['screen_attempts'])>0 for r in rr);metrics['removed_satellites']=sum(len(r['screened_satellites'].split(';')) for r in rr if r['screened_satellites']);metrics['successful_epochs_after_exclusion']=sum(r['status']=='PPP_FLOAT' and bool(r['screened_satellites']) for r in rr);report['comparisons'][stem]=initial['compare'](rr,outputs['off'],references,c['city'],scorer)
                report['cases'].setdefault(case,{})[pair+'-'+profile]=metrics;save(out/'report.json',report);print(case,pair,profile,metrics['status_counts'],metrics.get('screen_attempts'),metrics['reference_errors']['hypotheses'].get('reference_is_antenna'),flush=True)
    assert total==234536 and controls==12;assert all(sha(f)==h for f,h in manifest['input_sha256'].items());print('Completed 24 PPC profiles / 234536 rows and 12 exact legacy controls; parity unaccepted',flush=True)


if __name__=='__main__':main()
