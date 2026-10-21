"""Verify native standardized-residual traces, fixed PPC controls and all origins."""
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


def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verify-local-inputs',action='store_true');a=p.parse_args();root=Path(__file__).resolve().parent;repo=root.parents[2];m=json.loads((root/'predeclared-inputs.json').read_text());report=json.loads((root/'report.json').read_text());helper=runpy.run_path(str(repo/'scripts/run_ppp_standardized_screening.py'));baseline=runpy.run_path(str(repo/'scripts/run_ppc_ppp.py'));initial=runpy.run_path(str(repo/'scripts/run_ppc_ppp_initial_prior.py'));seed=runpy.run_path(str(repo/'scripts/run_ppp_motion_seeds.py'));scorer=runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))['score']
    assert list(m['profiles'])==['off','screen'] and len(m['cases'])==6 and len(m['pairs'])==2
    assert m['expected_profiles']==24 and m['expected_rows']==234536 and m['expected_controls']==12 and m['full_raw_rover_epochs']
    for k in ['retention_enabled','other_GNSS_gates_or_noise_changed','reference_or_IMU_used_by_solver','whole_native_PPP_executed','practical_parity_accepted','runtime_memory_acceptance']:assert not m[k]
    for f,h in m['source_sha256'].items():assert sha(repo/f)==h,f
    prior=repo/m['prior_archive'];assert sha(prior/'evidence-sha256.json')==m['prior_evidence_sha256']
    legacy=json.loads((root/'legacy-fixtures.json').read_text());assert legacy['base_commit']==m['base_commit'] and legacy['all_numeric_fixture_bytes_unchanged']
    for f,h in legacy['sha256'].items():assert sha(repo/f)==h,f
    if a.verify_local_inputs:
        for f,h in m['input_sha256'].items():assert sha(f)==h,f
    def capture(stem):
        e=json.loads((root/(stem+'.execution.json')).read_text());data=gzip.decompress((root/(stem+'.stdout.csv.gz')).read_bytes());err=gzip.decompress((root/(stem+'.stderr.log.gz')).read_bytes());assert e['returncode']==0 and hashlib.sha256(data).hexdigest()==e['csv_sha256'] and hashlib.sha256(err).hexdigest()==e['stderr_sha256'];assert b'truth_used_by_solver=false' in err;return data,err
    total=controls=profiles=0
    for case,c in m['cases'].items():
        references=None
        if a.verify_local_inputs:
            with Path(c['reference']).open(newline='') as f:refs=[{k.strip():v for k,v in r.items()} for r in csv.DictReader(f)]
            references={(int(r['GPS Week']),Decimal(r['GPS TOW (s)'])):r for r in refs};assert len(references)==len(refs)
        seed_control=None
        for pair in m['pairs']:
            outputs={}
            for profile in m['profiles']:
                stem=case+'-'+pair+'-'+profile;data,err=capture(stem);rr,metrics=baseline['inspect'](data);actual=report['cases'][case][pair+'-'+profile];outputs[profile]=rr;assert all(actual[k]==v for k,v in metrics.items()) and not actual['runtime_memory_acceptance']
                keys=[k for k in rr[0] if k.startswith('seed_')][:23];ss=[[r[k] for k in keys] for r in rr]
                if seed_control is None:seed_control=ss
                else:assert seed_control==ss
                assert actual['all_seed_source_fields_identical']
                if profile=='off':
                    previous=gzip.decompress((prior/(case+'-'+pair+'-doppler.stdout.csv.gz')).read_bytes());fields=next(csv.reader(io.StringIO(previous.decode())));assert len(fields)==72 and helper['legacy_columns'](data,fields)==previous and actual['all_72_prior_columns_identical_after_audited_label_normalization'];assert all(len(r)==72 for r in rr);controls+=1
                else:
                    assert b'standardized_residual_screen_sigma=4' in err and all(len(r)==75 for r in rr)
                    for r in rr:
                        d=r['screen_diagnostics'];n=int(r['screen_attempts']);removed=r['screened_satellites'].split(';') if r['screened_satellites'] else [];assert d.count('PppResidualScreeningDiagnostic {')==n and d.count('threshold_sigma: 4.0')==n and d.count('removed: true')==len(removed) and len(set(removed))==len(removed)
                        stats=list(map(float,re.findall(r'standardized_residual: ([^,]+)',d)));assert len(stats)==n and all(math.isfinite(v) and abs(v)>4 for v in stats)
                        if r['status']=='PPP_FLOAT':assert d.count('removed: false')==0
                    assert sum(int(r['screen_attempts']) for r in rr)==actual['screen_attempts'] and sum(int(r['screen_attempts'])>0 for r in rr)==actual['epochs_with_screen_attempts'] and sum(len(r['screened_satellites'].split(';')) for r in rr if r['screened_satellites'])==actual['removed_satellites'] and sum(r['status']=='PPP_FLOAT' and bool(r['screened_satellites']) for r in rr)==actual['successful_epochs_after_exclusion']
                    if references is not None:assert initial['compare'](rr,outputs['off'],references,c['city'],scorer)==report['comparisons'][stem]
                if references is not None:
                    assert scorer(initial['score_view'](rr),references,c['city'])==actual['reference_errors'];assert seed['seed_metrics'](rr,references,c['city'],scorer)==actual['seeds_and_service']
                r=json.loads((root/actual['resource_log']).read_text());assert r['returncode']==0 and r['peak_RSS_KiB']>0 and not r['performance_accepted'];total+=len(rr);profiles+=1
    assert total==234536 and profiles==24 and controls==12 and len(report['comparisons'])==12
    ng=json.loads((root/'native-native-generator.json').read_text());assert sha(root/'native-native-generator.json')==m['native_generator_sha256'] and ng['returncode']==0 and ng['pinned_upstream_commit']=='72f3b7c2c4c088dc575286a6f0edf4e407bc499a' and not ng['upstream_source_patched'] and not ng['whole_native_PPP_executed']
    assert ng['actual_worstStandardizedResidualRow_calls_per_repeat']==128 and ng['explicit_inverse_helper_calls_per_repeat']==64 and ng['production_LDLT_helper_calls_per_repeat']==64 and len(ng['repeat_runs'])==2 and ng['repeat_runs'][0]['fixture_sha256']==ng['repeat_runs'][1]['fixture_sha256']
    for f,h in ng['source_sha256'].items():assert sha(repo/f)==h,f
    for f,h in ng['repeat_runs'][0]['fixture_sha256'].items():assert sha(repo/'tests/fixtures'/f)==h,f
    if a.verify_local_inputs:
        for f,h in dict(ng['native_library_sha256'],**ng['native_private_header_sha256']).items():assert sha(f)==h,f
    old=json.loads((root/'initial-64-call-generator/native-generator.json').read_text());assert old['actual_worstStandardizedResidualRow_calls_per_repeat']==64 and old['repeat_runs'][0]['fixture_sha256']==old['repeat_runs'][1]['fixture_sha256'];assert old['repeat_runs'][0]['fixture_sha256']['upstream_ppp_screening.csv']==ng['repeat_runs'][0]['fixture_sha256']['upstream_ppp_screening.csv']
    q=json.loads((root/'qualification.json').read_text());assert q['predeclared_inputs_sha256']==sha(root/'predeclared-inputs.json') and q['passed_tests']==446 and {c['name'] for c in q['checks']}=={'tests','fmt','clippy','release','examples'}
    for c in q['checks']:
        data=gzip.decompress((root/c['log']).read_bytes());assert c['returncode']==0 and hashlib.sha256(data).hexdigest()==c['uncompressed_sha256']
        if c['name']=='tests':assert sum(map(int,re.findall(rb'test result: ok\. (\d+) passed;',data)))==446
    fail=root/'initial-inverse-validation';assert b'InvalidCovariance' in gzip.decompress((fail/'component-failure.log.gz').read_bytes()) and b'production Gaussian case 36' in gzip.decompress((fail/'production-LDLT-before-refinement.log.gz').read_bytes()) and b'3 passed' in gzip.decompress((fail/'all-inverse-refined-tests.log.gz').read_bytes())
    oracle=json.loads((fail/'high-precision-case36.json').read_text());assert oracle[0]['precision_decimal_digits']==90 and oracle[0]['exact_binary64_inputs']
    d=json.loads((root/'post-hoc-diagnostics.json').read_text());assert d['post_hoc_diagnostic_only'] and not d['used_to_select_profile'] and not d['SPP_relabelled_PPP'] and not d['practical_parity_accepted'] and len(d['native_GPS_1C2W_comparisons'])==24 and len(d['service_vs_off'])==12
    assert d['source_sha256']==sha(repo/'scripts/analyze_ppp_standardized_screening.py') and d['predeclared_inputs_sha256']==sha(root/'predeclared-inputs.json') and sha(repo/d['frozen_native_archive']/'evidence-sha256.json')==d['frozen_native_evidence_sha256']
    if a.verify_local_inputs:assert runpy.run_path(str(repo/'scripts/analyze_ppp_standardized_screening.py'))['derive'](root,repo)==d
    seal=json.loads((root/'evidence-sha256.json').read_text());assert set(seal)=={str(f.relative_to(root)) for f in root.rglob('*') if f.is_file() and f.name!='evidence-sha256.json'}
    for f,h in seal.items():assert sha(root/f)==h,f
    print('Verified 128 native screening calls per repeat, 446 tests, 24 PPC profiles / 234536 rows and 12 exact controls; parity unaccepted')


if __name__=='__main__':main()
