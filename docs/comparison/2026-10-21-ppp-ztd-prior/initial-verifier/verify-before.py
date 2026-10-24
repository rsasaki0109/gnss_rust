"""Verify initial-ZTD-only parameter study, unchanged library and prior controls."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import runpy


def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verify-local-inputs',action='store_true');a=p.parse_args();root=Path(__file__).resolve().parent;repo=root.parents[2];m=json.loads((root/'predeclared-inputs.json').read_text());report=json.loads((root/'report.json').read_text());helper=runpy.run_path(str(repo/'scripts/run_ppp_ztd_prior.py'));baseline=runpy.run_path(str(repo/'scripts/run_ppc_ppp.py'));initial=runpy.run_path(str(repo/'scripts/run_ppc_ppp_initial_prior.py'));seed=runpy.run_path(str(repo/'scripts/run_ppp_motion_seeds.py'));car=runpy.run_path(str(repo/'scripts/run_ppp_acceleration.py'));scorer=runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))['score']
    assert list(m['profiles'])==['default','sigma-0.6','sigma-0.12'] and len(m['cases'])==6 and len(m['pairs'])==2
    assert m['expected_PPC_profiles']==36 and m['expected_PPC_rows']==351804 and m['expected_exact_controls']==12 and m['expected_synthetic_profiles']==9 and m['full_raw_rover_epochs'] and m['library_processing_unchanged']
    for k in ['prior_mean_changed','process_or_measurement_noise_changed','existing_gates_changed','standardized_screening_enabled','retention_enabled','reference_or_IMU_used_by_solver','whole_native_PPP_executed','practical_parity_accepted','runtime_memory_acceptance']:assert not m[k]
    for f,h in dict(m['source_sha256'],**m['library_and_old_test_input_sha256']).items():assert sha(repo/f)==h,f
    prior=repo/m['prior_archive'];assert sha(prior/'evidence-sha256.json')==m['prior_evidence_sha256']
    if a.verify_local_inputs:
        for f,h in m['input_sha256'].items():assert sha(f)==h,f
    def capture(stem,ok=True):
        e=json.loads((root/(stem+'.execution.json')).read_text());data=gzip.decompress((root/(stem+'.stdout.csv.gz')).read_bytes());err=gzip.decompress((root/(stem+'.stderr.log.gz')).read_bytes());assert (e['returncode']==0)==ok and hashlib.sha256(data).hexdigest()==e['csv_sha256'] and hashlib.sha256(err).hexdigest()==e['stderr_sha256'];return data,err
    assert len(report['invalid_sigma_CLI'])==12
    for item in report['invalid_sigma_CLI']:
        data,err=capture(item['stem'],False);assert not data and b'must be finite and positive' in err and item['solver_was_not_run']
    truth=car['rows']((repo/'tests/fixtures/synthetic_ppp_acceleration_truth.csv').read_bytes());assert len(report['synthetic'])==9
    for kind in ['car-clean','car-noisy','car-unrounded']:
        for profile,s in m['profiles'].items():
            stem='synthetic-'+kind+'-'+profile;data,err=capture(stem);metrics=car['car_metrics'](data,truth,kind=='car-noisy',True);actual=report['synthetic'][stem];assert all(actual[k]==v for k,v in metrics.items());assert helper['prior_diagnostic'](err,s['sigma_m'],True)==actual['prior_diagnostic']
            if profile=='default':assert actual['all_prior_CSV_bytes_identical'] and actual['predeclared_targets_met'] and data==gzip.decompress((repo/'docs/comparison/2026-10-21-ppp-acceleration'/('synthetic-'+kind+'-pva.stdout.csv.gz')).read_bytes())
            elif kind!='car-noisy':assert not actual['predeclared_targets_met'] and actual['position_bound_failures']
    total=controls=profiles=0
    for case,c in m['cases'].items():
        references=None
        if a.verify_local_inputs:
            with Path(c['reference']).open(newline='') as f:refs=[{k.strip():v for k,v in r.items()} for r in csv.DictReader(f)]
            references={(int(r['GPS Week']),Decimal(r['GPS TOW (s)'])):r for r in refs};assert len(references)==len(refs)
        seed_control=None;mean=None
        for pair in m['pairs']:
            outputs={}
            for profile,s in m['profiles'].items():
                stem=case+'-'+pair+'-'+profile;data,err=capture(stem);rr,metrics=baseline['inspect'](data);actual=report['cases'][case][pair+'-'+profile];outputs[profile]=rr;assert all(actual[k]==v for k,v in metrics.items()) and not actual['runtime_memory_acceptance'] and all(len(r)==72 for r in rr);assert helper['prior_diagnostic'](err,s['sigma_m'])==actual['prior_diagnostic'];assert b'truth_used_by_solver=false' in err and b'standardized_residual_screen_sigma=' not in err
                cfg=re.search(rb'initial_head_variance_m2: \[([^\]]+)\]',err);assert list(map(float,cfg[1].split(b',')))==[100.,100.,100.,1e8,s['variance_m2']];observed=float(re.search(rb'initial_zenith_delay_m: ([^,}]+)',err)[1]);mean=observed if mean is None else mean;assert observed==mean==actual['initial_ZTD_mean_m']
                keys=[k for k in rr[0] if k.startswith('seed_')][:23];ss=[[r[k] for k in keys] for r in rr]
                if seed_control is None:seed_control=ss
                else:assert seed_control==ss
                assert actual['all_seed_source_fields_identical']
                if profile=='default':assert actual['all_72_prior_CSV_bytes_identical'] and data==gzip.decompress((prior/(case+'-'+pair+'-off.stdout.csv.gz')).read_bytes());controls+=1
                elif references is not None:assert initial['compare'](rr,outputs['default'],references,c['city'],scorer)==report['comparisons'][stem]
                if references is not None:assert scorer(initial['score_view'](rr),references,c['city'])==actual['reference_errors'] and seed['seed_metrics'](rr,references,c['city'],scorer)==actual['seeds_and_service']
                resource=json.loads((root/actual['resource_log']).read_text());assert resource['returncode']==0 and resource['peak_RSS_KiB']>0 and not resource['performance_accepted'];total+=len(rr);profiles+=1
    assert total==351804 and profiles==36 and controls==12 and len(report['comparisons'])==24
    ng=json.loads((root/'native-native-generator.json').read_text());assert sha(root/'native-native-generator.json')==m['native_generator_sha256'] and ng['returncode']==0 and ng['actual_initialTroposphereVariance_calls_per_repeat']==12 and ng['pinned_upstream_commit']=='72f3b7c2c4c088dc575286a6f0edf4e407bc499a' and not ng['upstream_source_patched'] and not ng['whole_native_PPP_executed'] and len(ng['repeat_runs'])==2 and ng['repeat_runs'][0]['fixture_sha256']==ng['repeat_runs'][1]['fixture_sha256']
    for f,h in ng['source_sha256'].items():assert sha(repo/f)==h,f
    for f,h in ng['repeat_runs'][0]['fixture_sha256'].items():assert sha(repo/'tests/fixtures'/f)==h,f
    with (repo/'tests/fixtures/upstream_ppp_ztd_prior.csv').open() as f:table={(int(r['madoca_per_frequency']),int(r['broadcast_model']),float(r['configured_variance_m2'])):float(r['returned_variance_m2']) for r in csv.DictReader(f)}
    assert len(table)==12 and table[0,0,0.36]==m['profiles']['default']['variance_m2'] and table[0,1,0.36]==m['profiles']['sigma-0.6']['variance_m2'] and table[1,0,0.36]==m['profiles']['sigma-0.12']['variance_m2']
    origin=json.loads((root/'prior-origin.json').read_text());assert origin['pinned_upstream_commit']==ng['pinned_upstream_commit'] and origin['source_archive_sha256']==ng['source_archive_sha256'] and origin['only_prior_values_borrowed'] and not origin['native_broadcast_or_MADOCA_whole_solver_executed']
    if a.verify_local_inputs:
        for f,h in dict(ng['native_library_sha256'],**ng['native_private_header_sha256'],**origin['headers_sha256']).items():assert sha(f)==h,f
    q=json.loads((root/'qualification.json').read_text());assert q['predeclared_inputs_sha256']==sha(root/'predeclared-inputs.json') and q['passed_tests']==446 and q['library_processing_unchanged'] and {c['name'] for c in q['checks']}=={'tests','fmt','clippy','release','examples'}
    for c in q['checks']:
        data=gzip.decompress((root/c['log']).read_bytes());assert c['returncode']==0 and hashlib.sha256(data).hexdigest()==c['uncompressed_sha256']
        if c['name']=='tests':assert sum(map(int,re.findall(rb'test result: ok\. (\d+) passed;',data)))==446
    d=json.loads((root/'post-hoc-diagnostics.json').read_text());assert d['post_hoc_diagnostic_only'] and not d['used_to_select_profile'] and not d['SPP_relabelled_PPP'] and not d['practical_parity_accepted'] and len(d['native_GPS_1C2W_comparisons'])==36 and len(d['service_vs_default'])==24
    assert d['source_sha256']==sha(repo/'scripts/analyze_ppp_ztd_prior.py') and d['predeclared_inputs_sha256']==sha(root/'predeclared-inputs.json') and sha(repo/d['frozen_native_archive']/'evidence-sha256.json')==d['frozen_native_evidence_sha256']
    if a.verify_local_inputs:assert runpy.run_path(str(repo/'scripts/analyze_ppp_ztd_prior.py'))['derive'](root,repo)==d
    seal=json.loads((root/'evidence-sha256.json').read_text());assert set(seal)=={str(f.relative_to(root)) for f in root.rglob('*') if f.is_file() and f.name!='evidence-sha256.json'}
    for f,h in seal.items():assert sha(root/f)==h,f
    print('Verified ZTD-only variance study: 36 PPC profiles / 351804 rows, 12 exact controls, nine car traces, 12 invalid CLI checks and 446 tests; parity unaccepted')


if __name__=='__main__':main()
