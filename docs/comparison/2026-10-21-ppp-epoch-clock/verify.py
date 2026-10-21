"""Verify every retained epoch-clock profile without changing solver results."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import json
from pathlib import Path
import re
import runpy


def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verify-local-inputs',action='store_true');a=p.parse_args();root=Path(__file__).resolve().parent;repo=root.parents[2]
    m=json.loads((root/'predeclared-inputs.json').read_text());report=json.loads((root/'report.json').read_text());prior=repo/m['prior_archive'];assert sha(prior/'evidence-sha256.json')==m['prior_evidence_sha256']
    for f,h in m['source_sha256'].items():assert sha(repo/f)==h,f
    if a.verify_local_inputs:
        for f,h in m['input_sha256'].items():assert sha(f)==h,f
    for k in ['thresholds_changed','ambiguity_prior_changed','SPP_or_Doppler_fused','truth_used_by_solver','default_clock_model_changed','whole_native_PPP_executed','practical_parity_accepted']:assert not m[k]
    assert m['full_raw_rover_epochs'] and len(m['cases'])==len(report['cases'])==6 and len(m['profiles'])==4 and len(report['clock_comparisons'])==24 and len(report['synthetic'])==8
    helper=runpy.run_path(str(repo/'scripts/run_ppp_epoch_clock.py'));baseline=runpy.run_path(str(repo/'scripts/run_ppc_ppp.py'));initial=runpy.run_path(str(repo/'scripts/run_ppc_ppp_initial_prior.py'));moving=runpy.run_path(str(repo/'scripts/run_ppp_kinematic.py'));scorer=runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))['score'];rows=baseline['rows'];truth=rows((repo/'tests/fixtures/synthetic_ppp_kinematic_truth.csv').read_bytes())
    def capture(stem):
        e=json.loads((root/(stem+'.execution.json')).read_text());data=gzip.decompress((root/(stem+'.stdout.csv.gz')).read_bytes());err=gzip.decompress((root/(stem+'.stderr.log.gz')).read_bytes());assert e['returncode']==0 and hashlib.sha256(data).hexdigest()==e['csv_sha256'] and hashlib.sha256(err).hexdigest()==e['stderr_sha256'];return data,err
    for variant in ['clean','noisy']:
        independent={}
        for jumped in [False,True]:
            kind=('jumped-' if jumped else '')+variant
            for clock in ['random-walk','independent']:
                stem='synthetic-'+kind+'-'+clock;data,err=capture(stem);assert moving['moving_metrics'](data,truth,True)==report['synthetic'][stem]
                if clock=='independent':
                    assert b'clock_model=IndependentEpoch' in err;independent[jumped]=rows(data);assert report['synthetic'][stem]['status_counts']=={'PPP_FLOAT':64}
                elif jumped:assert report['synthetic'][stem]['status_counts']=={'PPP_FLOAT':22,'INVALID':42}
        assert helper['stress_compare'](independent[False],independent[True])==report['stress_comparisons'][variant]
    totals={pair+'-'+profile:0 for pair in m['pairs'] for profile in m['profiles']}
    for case,c in m['cases'].items():
        references=None
        if a.verify_local_inputs:
            with Path(c['reference']).open(newline='') as f:refs=[{k.strip():v for k,v in r.items()} for r in csv.DictReader(f)]
            references={(int(r['GPS Week']),Decimal(r['GPS TOW (s)'])):r for r in refs}
        for pair in m['pairs']:
            outputs={}
            for profile,settings in m['profiles'].items():
                stem=case+'-'+pair+'-'+profile;data,err=capture(stem);rr,metrics=baseline['inspect'](data);actual=report['cases'][case][pair+'-'+profile];assert all(actual[k]==v for k,v in metrics.items());outputs[profile]=rr;totals[pair+'-'+profile]+=len(rr)
                if settings['clock_model']=='random-walk':assert actual['all_prior_CSV_bytes_identical'] and data==gzip.decompress((prior/(case+'-'+pair+'-'+settings['head_prior']+'.stdout.csv.gz')).read_bytes())
                else:
                    assert b'clock_model=IndependentEpoch' in err and b'SPP_fused=false' in err and b'thresholds_changed=false' in err
                    comparison=report['clock_comparisons'][stem];control=outputs[settings['head_prior']+'-random-walk'];old_keys={(r['week'],r['tow']) for r in control if r['status']=='PPP_FLOAT'};new_keys={(r['week'],r['tow']) for r in rr if r['status']=='PPP_FLOAT'}
                    assert (len(old_keys&new_keys),len(new_keys-old_keys),len(old_keys-new_keys))==(comparison['common_epochs'],comparison['newly_accepted_epochs'],comparison['lost_epochs'])
                    if references is not None:assert initial['compare'](rr,control,references,c['city'],scorer)==comparison
                score=actual['reference_errors'];assert score['valid_positions_without_reference']==0 and all(v['epochs']==metrics['status_counts'].get('PPP_FLOAT',0) for v in score['hypotheses'].values())
                if references is not None:assert scorer(initial['score_view'](rr),references,c['city'])==score
                resource=json.loads((root/actual['resource_log']).read_text());assert resource['returncode']==0 and resource['peak_RSS_KiB']>0 and not resource['performance_accepted'] and not actual['runtime_memory_acceptance']
    assert set(totals.values())=={58634}
    native=json.loads((root/'native-generator.json').read_text());assert native['build_returncode']==0 and native['actual_native_scalar_reset_calls_per_repeat']==16 and len(native['repeat_runs'])==2 and not native['whole_native_PPP_executed'] and not native['upstream_source_patched'] and native['source_sha256']==sha(repo/'scripts/upstream_ppp_epoch_clock.cpp')
    for r in native['repeat_runs']:assert r['returncode']==0 and r['fixture_sha256']==sha(repo/'tests/fixtures/upstream_ppp_epoch_clock.csv')
    q=json.loads((root/'qualification.json').read_text());assert q['predeclared_inputs_sha256']==sha(root/'predeclared-inputs.json') and {c['name'] for c in q['checks']}=={'rust','examples','fmt','clippy','release'}
    for c in q['checks']:
        data=gzip.decompress((root/c['log']).read_bytes());assert c['returncode']==0 and hashlib.sha256(data).hexdigest()==c['uncompressed_sha256']
        if c['name']=='rust':assert sum(map(int,re.findall(rb'test result: ok\. (\d+) passed;',data)))==c['passed_tests']==413
    worst=json.loads((root/'worst-points.json').read_text());assert worst['post_hoc_diagnostic_only'] and not worst['used_to_select_solver_profile'] and not worst['thresholds_changed'] and len(worst['records'])==60
    assert worst['source_sha256']==sha(repo/'scripts/analyze_ppp_epoch_clock.py') and worst['predeclared_inputs_sha256']==sha(root/'predeclared-inputs.json')
    if a.verify_local_inputs:assert runpy.run_path(str(repo/'scripts/analyze_ppp_epoch_clock.py'))['derive'](root,repo)==worst
    seal=json.loads((root/'evidence-sha256.json').read_text());assert set(seal)=={str(f.relative_to(root)) for f in root.rglob('*') if f.is_file() and f.name!='evidence-sha256.json'}
    for f,h in seal.items():assert sha(root/f)==h,f
    print('Verified 8 synthetic clock profiles, 48 full-rate PPC profiles, 24 exact controls, 16 native reset cases and 413 tests; practical parity unaccepted')


if __name__=='__main__':main()
