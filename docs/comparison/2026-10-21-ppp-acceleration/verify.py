"""Verify native PVA component provenance, synthetic bounds and full PPC controls."""
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
    m=json.loads((root/'predeclared-inputs.json').read_text());report=json.loads((root/'report.json').read_text());helper=runpy.run_path(str(repo/'scripts/run_ppp_acceleration.py'));baseline=runpy.run_path(str(repo/'scripts/run_ppc_ppp.py'));initial=runpy.run_path(str(repo/'scripts/run_ppc_ppp_initial_prior.py'));scorer=runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))['score']
    assert len(m['profiles'])==8 and len(m['cases'])==6 and m['full_raw_rover_epochs'] and len(m['pairs'])==2
    for k in ['thresholds_changed','default_position_model_changed','ambiguity_prior_changed','SPP_or_Doppler_fused','truth_used_by_solver','whole_native_PPP_executed','practical_parity_accepted','retention_enabled']:assert not m[k]
    for f,h in m['source_sha256'].items():assert sha(repo/f)==h,f
    prior=repo/m['prior_archive'];assert sha(prior/'evidence-sha256.json')==m['prior_evidence_sha256']
    legacy=json.loads((root/'legacy-ppp-inputs.json').read_text());assert legacy['collected_after_runs'] and legacy['base_commit']==m['base_commit'] and legacy['all_current_bytes_identical_to_base']
    for f,h in legacy['sha256'].items():assert sha(repo/f)==h,f
    if a.verify_local_inputs:
        for f,h in m['input_sha256'].items():assert sha(f)==h,f
    def capture(stem):
        e=json.loads((root/(stem+'.execution.json')).read_text());data=gzip.decompress((root/(stem+'.stdout.csv.gz')).read_bytes());err=gzip.decompress((root/(stem+'.stderr.log.gz')).read_bytes());assert e['returncode']==0 and hashlib.sha256(data).hexdigest()==e['csv_sha256'] and hashlib.sha256(err).hexdigest()==e['stderr_sha256'];return data,err
    truth=helper['rows']((repo/'tests/fixtures/synthetic_ppp_acceleration_truth.csv').read_bytes());assert len(report['synthetic'])==6 and len(report['old_synthetic_controls'])==2
    for kind in ['car-clean','car-noisy','car-unrounded']:
        for motion in ['kinematic','pva']:
            stem='synthetic-'+kind+'-'+motion;data,_=capture(stem);metrics=helper['car_metrics'](data,truth,kind=='car-noisy',motion=='pva');assert metrics==report['synthetic'][stem]
            if motion=='pva':assert metrics['predeclared_targets_met']
    for stem,c in report['old_synthetic_controls'].items():
        data,_=capture(stem);assert c['all_CSV_bytes_identical'] and data==gzip.decompress((repo/c['previous']).read_bytes())
    total=controls=profiles=0
    for case,c in m['cases'].items():
        references=None
        if a.verify_local_inputs:
            with Path(c['reference']).open(newline='') as f:refs=[{k.strip():v for k,v in r.items()} for r in csv.DictReader(f)]
            references={(int(r['GPS Week']),Decimal(r['GPS TOW (s)'])):r for r in refs};assert len(references)==len(refs)
        for pair in m['pairs']:
            outputs={}
            for profile,s in m['profiles'].items():
                stem=case+'-'+pair+'-'+profile;data,err=capture(stem);assert b'truth_used_by_solver=false' in err;rr,metrics=baseline['inspect'](data);helper['acceleration_fields'](rr,s['motion']=='PVA');actual=report['cases'][case][pair+'-'+profile];assert all(actual[k]==v for k,v in metrics.items());outputs[profile]=rr
                if s['motion']=='CV':
                    assert actual['all_prior_CSV_bytes_identical'] and data==gzip.decompress((prior/(case+'-'+pair+'-'+s['head']+'-'+s['clock']+'-interrupt.stdout.csv.gz')).read_bytes());controls+=1
                elif references is not None:assert initial['compare'](rr,outputs[s['head']+'-'+s['clock']+'-CV'],references,c['city'],scorer)==report['comparisons'][stem]
                if references is not None:assert scorer(initial['score_view'](rr),references,c['city'])==actual['reference_errors']
                resource=json.loads((root/actual['resource_log']).read_text());assert resource['returncode']==0 and resource['peak_RSS_KiB']>0 and not resource['performance_accepted'] and not actual['runtime_memory_acceptance'];total+=len(rr);profiles+=1
    assert total==938144 and profiles==96 and controls==48 and len(report['comparisons'])==48
    ng=json.loads((root/'native-generator.json').read_text());assert sha(root/'native-generator.json')==m['native_generator_sha256'] and ng['returncode']==0 and ng['pinned_upstream_commit']=='72f3b7c2c4c088dc575286a6f0edf4e407bc499a' and not ng['upstream_source_patched'] and not ng['whole_native_PPP_executed']
    assert ng['actual_PVA_prediction_calls_per_repeat']==48 and ng['actual_PVA_initialization_calls_per_repeat']==54 and ng['explicit_PPC_header_initializations_per_repeat']==6
    assert not gzip.decompress((root/'native-build.log.gz').read_bytes()) and len(ng['repeat_runs'])==2 and ng['repeat_runs'][0]['fixture_sha256']==ng['repeat_runs'][1]['fixture_sha256']
    for f,h in ng['source_sha256'].items():assert sha(repo/f)==h,f
    for f,h in ng['repeat_runs'][0]['fixture_sha256'].items():assert sha(repo/'tests/fixtures'/f)==h,f
    assert sha(repo/ng['prior_proof'])==ng['prior_proof_sha256']
    if a.verify_local_inputs:
        for f,h in ng['input_sha256'].items():assert sha(f)==h,f
        for f,h in ng['native_library_sha256'].items():assert sha(f)==h,f
    q=json.loads((root/'qualification.json').read_text());assert q['predeclared_inputs_sha256']==sha(root/'predeclared-inputs.json') and q['passed_tests']==427 and {c['name'] for c in q['checks']}=={'tests','fmt','clippy','release','examples'}
    for c in q['checks']:
        data=gzip.decompress((root/c['log']).read_bytes());assert c['returncode']==0 and hashlib.sha256(data).hexdigest()==c['uncompressed_sha256']
        if c['name']=='tests':assert sum(map(int,re.findall(rb'test result: ok\. (\d+) passed;',data)))==427
    partial=root/'initial-scope-count';pm=json.loads((partial/'predeclared-inputs.json').read_text());interrupted=json.loads((partial/'interruption.json').read_text());assert interrupted['intentional_interruption'] and not interrupted['numerical_solver_changed'] and interrupted['exit_code']==130
    assert all(pm[k]==m[k] for k in ['binary_sha256','cases','pairs','profiles','input_sha256'])
    for f,h in pm['source_sha256'].items():
        if f not in ['scripts/run_ppp_acceleration.py','docs/ppp-acceleration.md']:assert sha(repo/f)==h,f
    assert sha(partial/'runner-before-scope-count.py')==pm['source_sha256']['scripts/run_ppp_acceleration.py'] and sha(partial/'contract-before-scope-count.md')==pm['source_sha256']['docs/ppp-acceleration.md']
    for f in partial.glob('*.stdout.csv.gz'):assert gzip.decompress(f.read_bytes())==gzip.decompress((root/f.name).read_bytes()),f.name
    failure=root/'initial-harness';fm=json.loads((failure/'failed-build.json').read_text());assert fm['returncode']==1 and sha(failure/'source.cpp')==fm['source_sha256'] and b'navigation_internal.hpp' in gzip.decompress((failure/'failed-build.log.gz').read_bytes())
    d=json.loads((root/'post-hoc-diagnostics.json').read_text());assert d['post_hoc_diagnostic_only'] and not d['used_to_select_solver_profile'] and not d['native_SPP_relabelled_PPP'] and not d['practical_parity_accepted'] and len(d['native_GPS_1C2W_PPP_only_comparisons'])==48
    assert d['source_sha256']==sha(repo/'scripts/analyze_ppp_acceleration.py') and d['predeclared_inputs_sha256']==sha(root/'predeclared-inputs.json') and sha(repo/d['frozen_native_archive']/'evidence-sha256.json')==d['frozen_native_evidence_sha256']
    if a.verify_local_inputs:assert runpy.run_path(str(repo/'scripts/analyze_ppp_acceleration.py'))['derive'](root,repo)==d
    seal=json.loads((root/'evidence-sha256.json').read_text());assert set(seal)=={str(f.relative_to(root)) for f in root.rglob('*') if f.is_file() and f.name!='evidence-sha256.json'}
    for f,h in seal.items():assert sha(root/f)==h,f
    print('Verified PVA component/native fixtures, 8 synthetic runs, 96 PPC profiles / 938144 rows, 48 exact controls and 427 tests; parity unaccepted')


if __name__=='__main__':main()
