"""Freeze 36 full-rate PPC shadow/position/Doppler-seeded PVA comparisons."""
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


def projection(data,fields):
    rr=list(csv.DictReader(io.StringIO(data.decode())));out=io.StringIO();w=csv.DictWriter(out,fieldnames=fields,extrasaction='ignore',lineterminator='\n');w.writeheader();w.writerows(rr);return out.getvalue().encode()


def seed_metrics(rr,references,city,scorer):
    def view(which):
        result=[]
        for r in rr:
            status=r['service_status'] if which=='combined' else r['seed_status'] if which=='SPP_all' else 'SPP' if r['service_status']=='SPP' else 'INVALID'
            source='service' if which in ['combined','SPP_fallback_only'] else 'seed'
            result.append(dict(week=r['week'],tow=str(Decimal(r['tow']).quantize(Decimal('.000001'))),status='FLOAT' if status in ['SPP','PPP_FLOAT'] else 'INVALID',**{k+'_m':r[source+'_'+k+'_m'] for k in ['x','y','z']}))
        return result
    for r in rr:
        assert r['seed_status'] in ['SPP','INVALID'] and r['service_status'] in ['SPP','PPP_FLOAT','INVALID']
        if r['seed_status']=='SPP':assert all(math.isfinite(float(r[k])) for k in ['seed_x_m','seed_y_m','seed_z_m','seed_gps_clock_s','seed_position_covariance_trace_m2'])
        if r['seed_committed']=='true':assert r['status']=='PPP_FLOAT'
        if r['status']=='PPP_FLOAT':assert r['service_status']=='PPP_FLOAT' and all(r['service_'+k]==r[k] for k in ['x_m','y_m','z_m','vx_mps','vy_mps','vz_mps'])
        elif r['seed_status']=='SPP':assert r['service_status']=='SPP' and all(r['service_'+k]==r['seed_'+k] for k in ['x_m','y_m','z_m','vx_mps','vy_mps','vz_mps'])
        else:assert r['service_status']=='INVALID'
    return dict(seed_status_counts=dict(collections.Counter(r['seed_status'] for r in rr)),seed_velocity_available=sum(bool(r['seed_vx_mps']) for r in rr),seed_errors=dict(collections.Counter(r['seed_error'] for r in rr if r['seed_error'])),proposed_initializations=sum(r['seed_initialized']=='true' for r in rr),committed_initializations=sum(r['seed_initialized']=='true' and r['seed_committed']=='true' for r in rr),proposed_restarts=sum(r['seed_restarted']=='true' for r in rr),committed_restarts=sum(r['seed_restarted']=='true' and r['seed_committed']=='true' for r in rr),velocity_update_actions=dict(collections.Counter(r['seed_velocity_update'].split(' {')[0] for r in rr)),service_status_counts=dict(collections.Counter(r['service_status'] for r in rr)),reference_errors_by_output={kind:scorer(view(kind),references,city) for kind in ['SPP_all','SPP_fallback_only','combined']},SPP_relabelled_PPP=False)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--native-fixtures-dir',type=Path,required=True);a=p.parse_args();repo=Path(__file__).resolve().parents[1];out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    prior=repo/'docs/comparison/2026-10-21-ppp-acceleration';old=json.loads((prior/'predeclared-inputs.json').read_text());cases=old['cases'];antex=next(f for f in old['input_sha256'] if f.endswith('.atx'));binary=repo/'target/release/examples/ppc_ppp'
    for c in cases.values():c['navigation']=str(Path(c['rover']).with_name('base.nav'))
    native=a.native_fixtures_dir.resolve();ng=json.loads((native/'native-generator.json').read_text());assert ng['returncode']==0 and len(ng['repeat_runs'])==2 and ng['repeat_runs'][0]['fixture_sha256']==ng['repeat_runs'][1]['fixture_sha256']
    assert all(sha(repo/f)==h for f,h in ng['source_sha256'].items());assert all(sha(repo/'tests/fixtures'/f)==h for f,h in ng['repeat_runs'][0]['fixture_sha256'].items())
    sources=list((repo/'src').rglob('*.rs'))+[repo/f for f in ['examples/ppc_ppp.rs','examples/support/ppp_seed.rs','examples/ppc_ppp_seed_probe.rs','tests/ppp_kinematic.rs','tests/ppp_motion.rs','scripts/run_ppp_motion_seeds.py','scripts/build_ppp_motion_fixtures.py','scripts/upstream_ppp_motion.cpp','scripts/run_ppc_ppp.py','scripts/run_ppc_ppp_initial_prior.py','scripts/run_ppc_rtk_float.py','scripts/run_native_ppc_comparison.py','scripts/measure_ppp_process.py','docs/ppp-motion-seeds.md','Cargo.toml','Cargo.lock','rust-toolchain.toml']]+list((repo/'tests/fixtures').glob('*ppp_motion*'))
    inputs=dict(old['input_sha256']);inputs.update({c['navigation']:sha(c['navigation']) for c in cases.values()})
    profiles={n:dict(arguments=['native-head','-epoch-clock','-pva','-seed-'+n],SPP_position_used=n!='shadow',Doppler_used=n=='doppler') for n in ['shadow','position','doppler']}
    manifest=dict(base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),source_sha256={str(f.relative_to(repo)):sha(f) for f in sources},binary_sha256=sha(binary),input_sha256=inputs,cases=cases,pairs=old['pairs'],profiles=profiles,prior_archive=str(prior.relative_to(repo)),prior_evidence_sha256=sha(prior/'evidence-sha256.json'),native_generator_sha256=sha(native/'native-generator.json'),full_raw_rover_epochs=True,expected_profiles=36,expected_output_rows=351804,expected_exact_shadow_controls=12,clock='IndependentEpoch',head='native-head',motion='PVA',phase='CODE calibrated',retention_enabled=False,GNSS_thresholds_or_noise_changed=False,truth_used_by_solver=False,whole_native_PPP_executed=False,practical_parity_accepted=False,runtime_memory_acceptance=False,SPP_fallback_is_PPP=False,seed_coordinate_compatibility='broadcast WGS84 approximately compatible with IGS20: explicit caller assertion, no transformation',seed_receiver_reference='InstantaneousMarker;receiver-calibration-unknown',seed_reuse='same-epoch data-dependent heuristic, not independent position factor')
    helper=runpy.run_path(str(repo/'scripts/run_native_ppc_comparison.py'));save,retain=[helper[k] for k in ['save','retain']];baseline=runpy.run_path(str(repo/'scripts/run_ppc_ppp.py'));initial=runpy.run_path(str(repo/'scripts/run_ppc_ppp_initial_prior.py'));scorer=runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))['score']
    assert all(sha(f)==h for f,h in inputs.items());save(out/'predeclared-inputs.json',manifest)
    for name in ['native-generator.json','build.log.gz']+[f'run{i}.{stream}.log.gz' for i in range(2) for stream in ['stdout','stderr']]:(out/('native-'+name)).write_bytes((native/name).read_bytes())
    report=dict(cases={},comparisons={},practical_parity_accepted=False);total=controls=0
    for case,c in cases.items():
        with Path(c['reference']).open(newline='') as f:refs=[{k.strip():v for k,v in r.items()} for r in csv.DictReader(f)]
        references={(int(r['GPS Week']),Decimal(r['GPS TOW (s)'])):r for r in refs};assert len(references)==len(refs);seed_control=None
        for pair in manifest['pairs']:
            outputs={}
            for profile,s in profiles.items():
                assert all(sha(repo/f)==h for f,h in manifest['source_sha256'].items()) and sha(binary)==manifest['binary_sha256']
                stem=case+'-'+pair+'-'+profile;command=[str(binary)]+[c[k] for k in ['rover','orbit','clock','bias']]+[antex,case,pair,'calibrated']+s['arguments']+[c['navigation']];wrapper=[sys.executable,str(repo/'scripts/measure_ppp_process.py'),str(out/(stem+'.resource.json'))]
                r=subprocess.run(wrapper+command,capture_output=True,cwd=repo);retain(out,stem,command,r.stdout,r.stderr,r.returncode);assert r.returncode==0,r.stderr.decode()
                rr,metrics=baseline['inspect'](r.stdout);total+=len(rr);outputs[profile]=rr;metrics['reference_errors']=scorer(initial['score_view'](rr),references,c['city']);metrics['seeds_and_service']=seed_metrics(rr,references,c['city'],scorer);metrics['resource_log']=stem+'.resource.json';metrics['runtime_memory_acceptance']=False
                seed_fields=[k for k in rr[0] if k.startswith('seed_')][:23];seeds=[[r[k] for k in seed_fields] for r in rr]
                if seed_control is None:seed_control=seeds
                else:assert seeds==seed_control,(stem,'seed source depends on PPP mode or IF pair')
                metrics['all_independent_seed_fields_identical_across_profiles']=True
                if profile=='shadow':
                    previous=prior/(case+'-'+pair+'-native-head-independent-PVA.stdout.csv.gz');data=gzip.decompress(previous.read_bytes());fields=next(csv.reader(io.StringIO(data.decode())));assert len(fields)==35 and projection(r.stdout,fields)==data,stem;metrics['all_35_prior_PPP_columns_byte_identical']=True;controls+=1
                    assert all(not r['seed_committed'] for r in rr)
                else:report['comparisons'][stem]=initial['compare'](rr,outputs['shadow'],references,c['city'],scorer)
                report['cases'].setdefault(case,{})[pair+'-'+profile]=metrics;save(out/'report.json',report);print(case,pair,profile,metrics['status_counts'],metrics['seeds_and_service']['service_status_counts'],metrics['reference_errors']['hypotheses'].get('reference_is_antenna'),flush=True)
    assert total==351804 and controls==12;assert all(sha(f)==h for f,h in inputs.items());print('Completed 36 PPC profiles, 351804 rows and 12 exact 35-column controls; practical parity unaccepted',flush=True)


if __name__=='__main__':main()
