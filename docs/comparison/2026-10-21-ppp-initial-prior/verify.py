"""Read-only recomputation of all 48 initial-prior PPC profiles and controls."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import runpy


def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verify-local-inputs',action='store_true');a=p.parse_args()
    root=Path(__file__).resolve().parent;repo=root.parents[2];m=json.loads((root/'predeclared-inputs.json').read_text());report=json.loads((root/'report.json').read_text())
    prior=repo/m['prior_archive'];assert sha(prior/'evidence-sha256.json')==m['prior_evidence_sha256']
    for f,h in m['source_sha256'].items():assert sha(repo/f)==h,f
    if a.verify_local_inputs:
        for f,h in m['input_sha256'].items():assert sha(f)==h,f
    for k in ['thresholds_changed','process_noise_changed','ambiguity_prior_changed','truth_used_by_solver','whole_native_PPP_executed','practical_parity_accepted']:assert not m[k]
    assert m['full_raw_rover_epochs'] and len(m['profiles'])==4 and len(m['pairs'])==2 and len(m['cases'])==len(report['cases'])==6 and len(report['comparisons'])==36
    helper=runpy.run_path(str(repo/'scripts/run_ppc_ppp_initial_prior.py'));inspect=runpy.run_path(str(repo/'scripts/run_ppc_ppp.py'))['inspect'];scorer=runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))['score']
    totals={pair+'-'+profile:0 for pair in m['pairs'] for profile in m['profiles']}
    for case,c in m['cases'].items():
        references=None
        if a.verify_local_inputs:
            with Path(c['reference']).open(newline='') as f:refs=[{k.strip():v for k,v in r.items()} for r in csv.DictReader(f)]
            references={(int(r['GPS Week']),Decimal(r['GPS TOW (s)'])):r for r in refs}
        for pair in m['pairs']:
            outputs={}
            for profile,settings in m['profiles'].items():
                stem=case+'-'+pair+'-'+profile;e=json.loads((root/(stem+'.execution.json')).read_text());data=gzip.decompress((root/(stem+'.stdout.csv.gz')).read_bytes());err=gzip.decompress((root/(stem+'.stderr.log.gz')).read_bytes());assert e['returncode']==0 and hashlib.sha256(data).hexdigest()==e['csv_sha256'] and hashlib.sha256(err).hexdigest()==e['stderr_sha256']
                rr,metrics=inspect(data);actual=report['cases'][case][pair+'-'+profile];assert all(actual[k]==v for k,v in metrics.items());outputs[profile]=rr;totals[pair+'-'+profile]+=len(rr)
                if profile=='default':assert actual['all_prior_CSV_bytes_identical'] and data==gzip.decompress((prior/(case+'-'+pair+'-calibrated.stdout.csv.gz')).read_bytes())
                else:assert ('initial_head_profile='+profile).encode() in err and b'thresholds_changed=false' in err and b'ambiguity_prior_changed=false' in err
                resource=json.loads((root/actual['resource_log']).read_text());assert resource['returncode']==0 and resource['peak_RSS_KiB']>0 and not resource['performance_accepted'] and not actual['runtime_memory_acceptance']
                score=actual['reference_errors'];assert score['valid_positions_without_reference']==0 and all(r['epochs']==metrics['status_counts'].get('PPP_FLOAT',0) for r in score['hypotheses'].values())
                if references is not None:assert scorer(helper['score_view'](rr),references,c['city'])==score
            aa={(r['week'],r['tow']):r for r in outputs['default'] if r['status']=='PPP_FLOAT'}
            for profile in m['profiles']:
                if profile=='default':continue
                actual=report['comparisons'][case+'-'+pair+'-'+profile]
                bb={(r['week'],r['tow']):r for r in outputs[profile] if r['status']=='PPP_FLOAT'};common=aa.keys()&bb.keys();added=bb.keys()-aa.keys();lost=aa.keys()-bb.keys()
                assert (len(common),len(added),len(lost))==(actual['common_epochs'],actual['newly_accepted_epochs'],actual['lost_epochs'])
                diff=[math.sqrt(sum((float(aa[t][k])-float(bb[t][k]))**2 for k in ['x_m','y_m','z_m'])) for t in sorted(common)]
                assert actual['common_position_difference_RMS_m']==(math.sqrt(sum(v*v for v in diff)/len(diff)) if diff else None) and actual['common_position_difference_maximum_m']==(max(diff) if diff else None)
                if references is not None:assert helper['compare'](outputs[profile],outputs['default'],references,c['city'],scorer)==actual
    assert set(totals.values())=={58634}
    native=json.loads((root/'native-fixtures.json').read_text());assert len(native['builds'])==2 and len(native['repeat_runs'])==4 and not native['whole_native_PPP_executed'] and not native['upstream_source_patched']
    for b in native['builds']:assert b['returncode']==0 and sha(root/b['log'])==b['log_sha256'] and sha(repo/b['source'])==b['source_sha256']
    for r in native['repeat_runs']:assert r['returncode']==0 and r['actual_native_component_calls']==(312 if r['kind']=='climatology' else 6) and sha(repo/'tests/fixtures'/r['fixture'])==r['fixture_sha256']
    assert native['upstream_commit']=='72f3b7c2c4c088dc575286a6f0edf4e407bc499a'
    q=json.loads((root/'qualification.json').read_text());assert q['predeclared_inputs_sha256']==sha(root/'predeclared-inputs.json') and {c['name'] for c in q['checks']}=={'rust','examples','fmt','clippy','release'}
    for c in q['checks']:
        data=gzip.decompress((root/c['log']).read_bytes());assert c['returncode']==0 and hashlib.sha256(data).hexdigest()==c['uncompressed_sha256']
        if c['name']=='rust':assert sum(map(int,re.findall(rb'test result: ok\. (\d+) passed;',data)))==c['passed_tests']==410
    sealed=json.loads((root/'evidence-sha256.json').read_text());assert set(sealed)=={str(f.relative_to(root)) for f in root.rglob('*') if f.is_file() and f.name!='evidence-sha256.json'}
    for f,h in sealed.items():assert sha(root/f)==h,f
    print('Verified 48 full 5 Hz initial-prior profiles, 469072 output rows, 12 exact controls, 312+6 native component cases and 410 tests; practical parity unaccepted')


if __name__=='__main__':main()
