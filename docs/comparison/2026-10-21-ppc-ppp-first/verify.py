"""Recompute the first full-rate PPC PPP diagnostic from every retained CSV."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import runpy
import math
from decimal import Decimal
import csv


def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verify-local-inputs',action='store_true');a=p.parse_args()
    root=Path(__file__).resolve().parent;repo=root.parents[2];m=json.loads((root/'predeclared-inputs.json').read_text());report=json.loads((root/'report.json').read_text());helper=runpy.run_path(str(repo/'scripts/run_ppc_ppp.py'))
    for k in ['base_observations_used','SPP_or_Doppler_fused','truth_used_by_solver','whole_native_PPP_executed','runtime_memory_acceptance','practical_parity_accepted']:assert not m[k]
    assert m['full_raw_rover_epochs'] and m['reference_used_only_after_solver'] and not report['whole_native_PPP_executed'] and not report['practical_parity_accepted']
    assert len(m['cases'])==len(report['cases'])==6 and len(m['profiles'])==4 and len(report['common_success_comparisons'])==12
    for f,h in m['source_sha256'].items():assert sha(repo/f)==h,f
    if a.verify_local_inputs:
        for f,h in m['input_sha256'].items():assert sha(Path(f))==h,f
    all_rows={};totals={profile:0 for profile in m['profiles']}
    for case,c in m['cases'].items():
        assert set(report['cases'][case])==set(m['profiles'])
        for profile in m['profiles']:
            stem=case+'-'+profile;e=json.loads((root/(stem+'.execution.json')).read_text());data=gzip.decompress((root/(stem+'.stdout.csv.gz')).read_bytes());err=gzip.decompress((root/(stem+'.stderr.log.gz')).read_bytes())
            assert e['returncode']==0 and hashlib.sha256(data).hexdigest()==e['csv_sha256'] and hashlib.sha256(err).hexdigest()==e['stderr_sha256']
            assert b'truth_used_by_solver=false' in err and b'receiver_calibration_known=false' in err and b'satellite_PCV_applied=false' in err and b'solid_earth_tide_applied=false' in err
            rr,metrics=helper['inspect'](data);actual=report['cases'][case][profile];assert all(actual[k]==v for k,v in metrics.items());totals[profile]+=len(rr);all_rows[case,profile]=rr
            rounding=max(abs(Decimal(r['tow']).quantize(Decimal('.000001'))-Decimal(r['tow'])) for r in rr);assert str(rounding)==actual['reference_time_rounding_maximum_s'] and rounding<=Decimal('.000001')
            resource=json.loads((root/actual['resource_log']).read_text());assert resource['returncode']==0 and resource['peak_RSS_KiB']>0 and resource['wall_seconds']>0 and not resource['performance_accepted'] and not actual['runtime_memory_acceptance']
            success=metrics['status_counts'].get('PPP_FLOAT',0);scores=actual['reference_errors'];assert scores['valid_positions_without_reference']==0
            expected={'reference_is_antenna','reference_is_IMU_nagoya'} if c['city']=='nagoya' else {'reference_is_antenna','reference_is_IMU_tokyo_top_README','reference_is_IMU_tokyo_city_README'}
            assert set(scores['hypotheses'])==expected if success else scores['hypotheses']=={}
            assert all(v['epochs']==success for v in scores['hypotheses'].values())
            if a.verify_local_inputs:
                with Path(c['reference']).open(newline='') as f:refs=[{k.strip():v for k,v in r.items()} for r in csv.DictReader(f)]
                refs={(int(r['GPS Week']),Decimal(r['GPS TOW (s)'])):r for r in refs};view=[dict(r,tow=str(Decimal(r['tow']).quantize(Decimal('.000001'))),status='FLOAT') if r['status']=='PPP_FLOAT' else r for r in rr]
                scorer=runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))['score'];assert scorer(view,refs,c['city'])==scores
        for pair in ['1C2W','1C2L']:
            aa,bb=[{(r['week'],r['tow']):r for r in all_rows[case,pair+'-'+phase] if r['status']=='PPP_FLOAT'} for phase in ['calibrated','raw-float']];common=aa.keys()&bb.keys();diff=[math.sqrt(sum((float(aa[t][k])-float(bb[t][k]))**2 for k in ['x_m','y_m','z_m'])) for t in common]
            expected=dict(common_epochs=len(common),calibrated_only=len(aa.keys()-bb.keys()),raw_only=len(bb.keys()-aa.keys()),position_difference_3D_RMS_m=math.sqrt(sum(v*v for v in diff)/len(diff)) if diff else None,position_difference_maximum_m=max(diff) if diff else None)
            actual=report['common_success_comparisons'][case+'-'+pair]
            for key,value in expected.items():
                if isinstance(value,float):assert math.isclose(actual[key],value,rel_tol=1e-14,abs_tol=1e-20)
                else:assert actual[key]==value
    assert set(totals.values())=={58634}
    q=json.loads((root/'qualification.json').read_text());assert q['predeclared_inputs_sha256']==sha(root/'predeclared-inputs.json') and {c['name'] for c in q['checks']}=={'examples','fmt','clippy','release'}
    for c in q['checks']:
        data=gzip.decompress((root/c['log']).read_bytes());assert c['returncode']==0 and hashlib.sha256(data).hexdigest()==c['uncompressed_sha256']
    inherited=q['unchanged_library_tests'];prior=repo/inherited['archive'];assert sha(prior/'qualification.json')==inherited['qualification_sha256'] and inherited['passed_tests']==406
    old=json.loads((prior/'predeclared-inputs.json').read_text());assert all(m['source_sha256'][s]==h for s,h in old['source_sha256'].items() if s.startswith('src/'))
    native=json.loads((root/'native-initial.json').read_text());assert native['returncode']==native['build_returncode']==0 and native['actual_initializeFilter_calls']==6 and not native['whole_native_PPP_executed'] and not native['upstream_source_patched']
    for f,h in native['source_sha256'].items():assert sha(repo/f)==h,f
    for label in ['stdout','stderr']:
        data=gzip.decompress((root/('native-initial.'+label+('.csv.gz' if label=='stdout' else '.log.gz'))).read_bytes());assert hashlib.sha256(data).hexdigest()==native[label+'_sha256']
    initial=helper['rows'](gzip.decompress((root/'native-initial.stdout.csv.gz').read_bytes()));assert len(initial)==6
    for r in initial:assert float(r['clock_m'])==0 and float(r['clock_q_m2'])==1e8 and float(r['zenith_q_m2'])==float(r['velocity_qxx_m2_s2'])==25 and all(float(r[k])==100 for k in ['pos_qxx_m2','pos_qyy_m2','pos_qzz_m2'])
    seal=json.loads((root/'evidence-sha256.json').read_text());assert set(seal)=={str(f.relative_to(root)) for f in root.rglob('*') if f.is_file() and f.name!='evidence-sha256.json'}
    for f,h in seal.items():assert sha(root/f)==h,f
    print('Verified 24 full 5 Hz PPC PPP diagnostics, 234536 output rows, both origin hypotheses and unchanged 406-test library; practical parity unaccepted')


if __name__=='__main__':main()
