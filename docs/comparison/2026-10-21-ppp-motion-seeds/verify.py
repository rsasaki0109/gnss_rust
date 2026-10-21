"""Verify native motion components, 36 PPC captures and separate PPP/SPP scores."""
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
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verify-local-inputs',action='store_true');a=p.parse_args();root=Path(__file__).resolve().parent;repo=root.parents[2];m=json.loads((root/'predeclared-inputs.json').read_text());report=json.loads((root/'report.json').read_text());helper=runpy.run_path(str(repo/'scripts/run_ppp_motion_seeds.py'));baseline=runpy.run_path(str(repo/'scripts/run_ppc_ppp.py'));initial=runpy.run_path(str(repo/'scripts/run_ppc_ppp_initial_prior.py'));scorer=runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))['score']
    assert len(m['cases'])==6 and len(m['pairs'])==2 and list(m['profiles'])==['shadow','position','doppler']
    assert m['expected_profiles']==36 and m['expected_output_rows']==351804 and m['expected_exact_shadow_controls']==12 and m['full_raw_rover_epochs']
    for k in ['retention_enabled','GNSS_thresholds_or_noise_changed','truth_used_by_solver','whole_native_PPP_executed','practical_parity_accepted','runtime_memory_acceptance','SPP_fallback_is_PPP']:assert not m[k]
    label=json.loads((root/'clock-label-fix/audit.json').read_text());assert label['metadata_only_change'] and not label['GNSS_numeric_processing_changed'] and label['library_sources_unchanged']
    for f,h in m['source_sha256'].items():
        if f in label['changed_sources']:
            v=label['changed_sources'][f];assert sha(root/v['captured_path'])==h==v['before_sha256'] and sha(repo/f)==v['after_sha256'],f
        else:assert sha(repo/f)==h,f
    prior=repo/m['prior_archive'];assert sha(prior/'evidence-sha256.json')==m['prior_evidence_sha256']
    legacy=json.loads((root/'legacy-fixtures.json').read_text());assert legacy['base_commit']==m['base_commit'] and legacy['all_current_bytes_identical_to_base']
    for f,h in legacy['sha256'].items():assert sha(repo/f)==h,f
    assert sha(repo/legacy['documentation_change']['path'])==legacy['documentation_change']['current_sha256'] and not legacy['documentation_change']['numeric_fixture_bytes_changed']
    if a.verify_local_inputs:
        for f,h in m['input_sha256'].items():assert sha(f)==h,f
    def capture(stem):
        e=json.loads((root/(stem+'.execution.json')).read_text());data=gzip.decompress((root/(stem+'.stdout.csv.gz')).read_bytes());err=gzip.decompress((root/(stem+'.stderr.log.gz')).read_bytes());assert e['returncode']==0 and hashlib.sha256(data).hexdigest()==e['csv_sha256'] and hashlib.sha256(err).hexdigest()==e['stderr_sha256'];assert b'truth_used_by_solver=false' in err and b'seed_clock_imported=false' in err;return data
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
                stem=case+'-'+pair+'-'+profile;data=capture(stem);rr,metrics=baseline['inspect'](data);outputs[profile]=rr;actual=report['cases'][case][pair+'-'+profile];assert all(actual[k]==v for k,v in metrics.items()) and not actual['runtime_memory_acceptance'];assert all(len(r)==72 and None not in r for r in rr)
                keys=[k for k in rr[0] if k.startswith('seed_')][:23];seeds=[[r[k] for k in keys] for r in rr]
                if seed_control is None:seed_control=seeds
                else:assert seeds==seed_control
                assert actual['all_independent_seed_fields_identical_across_profiles']
                if profile=='shadow':
                    previous=gzip.decompress((prior/(case+'-'+pair+'-native-head-independent-PVA.stdout.csv.gz')).read_bytes());fields=next(csv.reader(io.StringIO(previous.decode())));assert len(fields)==35 and helper['projection'](data,fields)==previous and actual['all_35_prior_PPP_columns_byte_identical'];assert all(not r['seed_committed'] for r in rr);controls+=1
                else:
                    for r in rr:
                        if r['status']=='PPP_FLOAT':assert r['seed_status']=='SPP' and r['seed_committed']=='true'
                        if r['seed_committed']=='true':assert r['status']=='PPP_FLOAT'
                        if profile=='position' and r['seed_velocity_update']:assert r['seed_velocity_update']=='Disabled'
                    if references is not None:assert initial['compare'](rr,outputs['shadow'],references,c['city'],scorer)==report['comparisons'][stem]
                if references is not None:
                    assert scorer(initial['score_view'](rr),references,c['city'])==actual['reference_errors'];assert helper['seed_metrics'](rr,references,c['city'],scorer)==actual['seeds_and_service']
                r=json.loads((root/actual['resource_log']).read_text());assert r['returncode']==0 and r['peak_RSS_KiB']>0 and not r['performance_accepted'];total+=len(rr);profiles+=1
    assert total==351804 and profiles==36 and controls==12 and len(report['comparisons'])==24
    ng=json.loads((root/'native-native-generator.json').read_text());assert sha(root/'native-native-generator.json')==m['native_generator_sha256'] and ng['returncode']==0 and ng['pinned_upstream_commit']=='72f3b7c2c4c088dc575286a6f0edf4e407bc499a' and not ng['upstream_source_patched'] and not ng['whole_native_PPP_executed']
    for key in ['actual_predictState_dt_zero_calls_per_repeat','actual_initializeFilter_calls_per_repeat','direct_velocity_helper_calls_per_repeat','direct_position_disagreement_calls_per_repeat']:assert ng[key]==32
    assert not gzip.decompress((root/'native-build.log.gz').read_bytes()) and len(ng['repeat_runs'])==2 and ng['repeat_runs'][0]['fixture_sha256']==ng['repeat_runs'][1]['fixture_sha256']
    for f,h in ng['source_sha256'].items():assert sha(repo/f)==h,f
    for f,h in ng['repeat_runs'][0]['fixture_sha256'].items():assert sha(repo/'tests/fixtures'/f)==h,f
    assert sha(repo/ng['prior_proof'])==ng['prior_proof_sha256']
    if a.verify_local_inputs:
        for f,h in dict(ng['native_library_sha256'],**ng['native_private_header_sha256']).items():assert sha(f)==h,f
    q=json.loads((root/'qualification.json').read_text());assert q['predeclared_inputs_sha256']==sha(root/'predeclared-inputs.json') and q['passed_tests']==437 and {c['name'] for c in q['checks']}=={'tests','fmt','clippy','release','examples'}
    for c in q['checks']:
        data=gzip.decompress((root/c['log']).read_bytes());assert c['returncode']==0 and hashlib.sha256(data).hexdigest()==c['uncompressed_sha256']
        if c['name']=='tests':assert sum(map(int,re.findall(rb'test result: ok\. (\d+) passed;',data)))==437
    assert not q['library_sources_changed_after_full_tests'] and {c['name'] for c in q['metadata_fix_checks']}=={'fmt','clippy','examples'}
    for c in q['metadata_fix_checks']:
        data=gzip.decompress((root/c['log']).read_bytes());assert c['returncode']==0 and hashlib.sha256(data).hexdigest()==c['uncompressed_sha256']
    assert label['source_sha256']==sha(repo/'scripts/audit_ppp_seed_clock_labels.py') and set(label['cases'])==set(m['cases'])
    audited=0
    for case,metrics in label['cases'].items():
        stem=case+'-1C2W-shadow';old=list(csv.DictReader(io.StringIO(gzip.decompress((root/(stem+'.stdout.csv.gz')).read_bytes()).decode())))
        e=json.loads((root/'clock-label-fix'/(case+'.execution.json')).read_text());raw=gzip.decompress((root/'clock-label-fix'/(case+'.stdout.csv.gz')).read_bytes());err=gzip.decompress((root/'clock-label-fix'/(case+'.stderr.log.gz')).read_bytes());assert e['returncode']==0 and hashlib.sha256(raw).hexdigest()==e['csv_sha256'] and hashlib.sha256(err).hexdigest()==e['stderr_sha256']
        new=list(csv.DictReader(io.StringIO(raw.decode())));assert len(old)==len(new)==metrics['input_epochs'] and metrics['all_72_fields_equal_after_metadata_normalization'];clocks={}
        for aa,bb in zip(old,new):
            bb['seed_gps_clock_s']=bb.pop('seed_first_system_clock_s')
            if bb['seed_status']=='SPP':
                source,clock=bb['seed_source'].rsplit(';SPP_first_clock_system=',1);assert clock in ['Gps','Galileo'];clocks[clock]=clocks.get(clock,0)+1;bb['seed_source']=source
            assert aa==bb
        assert clocks==metrics['first_SPP_clock_system_counts'];audited+=len(new)
    assert audited==58634
    failure=root/'initial-native-harness';fm=json.loads((failure/'failed-build.json').read_text());assert fm['returncode']!=0 and sha(failure/'source.cpp')==fm['source_sha256'] and b'misleading-indentation' in gzip.decompress((failure/'failed-build.log.gz').read_bytes())
    probe=root/'read-only-worst-point-probe';pm=json.loads((probe/'execution.json').read_text());pe=json.loads((probe/'reference-evaluation.json').read_text());assert pm['returncode']==0 and not pm['PPP_filter_updated'] and not pm['truth_used_by_solver'] and sha(repo/'examples/ppc_ppp_seed_probe.rs')==pm['source_sha256'] and pe['all_four_profiles_retained'] and len(pe['scores'])==4
    repeat=json.loads((probe/'repeat-current-binary.json').read_text());assert repeat['repeat_stdout_and_stderr_byte_identical'] and repeat['returncode']==0 and not repeat['PPP_filter_updated'] and repeat['source_sha256']==pm['source_sha256']
    d=json.loads((root/'post-hoc-diagnostics.json').read_text());assert d['post_hoc_diagnostic_only'] and not d['used_to_select_profile'] and not d['SPP_relabelled_PPP'] and not d['practical_parity_accepted'] and len(d['native_GPS_1C2W_comparisons'])==36 and len(d['service_vs_shadow'])==24
    assert d['source_sha256']==sha(repo/'scripts/analyze_ppp_motion_seeds.py') and d['predeclared_inputs_sha256']==sha(root/'predeclared-inputs.json') and sha(repo/d['frozen_native_archive']/'evidence-sha256.json')==d['frozen_native_evidence_sha256']
    if a.verify_local_inputs:
        assert runpy.run_path(str(repo/'scripts/analyze_ppp_motion_seeds.py'))['derive'](root,repo)==d
        assert sha(m['cases']['nagoya-run3']['reference'])==pe['reference_sha256']
        with Path(m['cases']['nagoya-run3']['reference']).open() as f:refs=[{k.strip():v for k,v in r.items()} for r in csv.DictReader(f)]
        references={(int(r['GPS Week']),Decimal(r['GPS TOW (s)'])):r for r in refs}
        for r in pe['raw_solver_rows']:
            view=dict(r,tow=str(Decimal(r['tow']).quantize(Decimal('.000001'))),status='FLOAT' if r['status']=='SPP' else 'INVALID');assert scorer([view],references,'nagoya')==pe['scores'][r['profile']]['reference_errors']
    seal=json.loads((root/'evidence-sha256.json').read_text());assert set(seal)=={str(f.relative_to(root)) for f in root.rglob('*') if f.is_file() and f.name!='evidence-sha256.json'}
    for f,h in seal.items():assert sha(root/f)==h,f
    print('Verified native motion components, 437 tests, 36 PPC profiles / 351804 rows, 12 exact controls and separate PPP/SPP evidence; parity unaccepted')


if __name__=='__main__':main()
