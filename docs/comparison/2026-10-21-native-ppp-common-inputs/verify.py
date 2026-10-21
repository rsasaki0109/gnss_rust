"""Recompute whole native/Rust PPC results and shared decoded-input readbacks."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import json
from pathlib import Path
import runpy


def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verify-local-inputs',action='store_true');a=p.parse_args();root=Path(__file__).resolve().parent;repo=root.parents[2]
    m=json.loads((root/'predeclared-inputs.json').read_text());report=json.loads((root/'report.json').read_text());helper=runpy.run_path(str(repo/'scripts/run_native_ppc_ppp_comparison.py'));baseline=runpy.run_path(str(repo/'scripts/run_ppc_ppp.py'));scorer=runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))['score'];compare=runpy.run_path(str(repo/'scripts/run_ppc_ppp_initial_prior.py'))['compare']
    assert m['whole_native_PPP_processEpoch_executed'] and m['shared_decoded_input'] and m['full_raw_rover_epochs']
    for k in ['native_RINEX_observation_pipeline_executed','native_solver_source_patched','truth_used_by_solver','practical_parity_accepted','timing_scopes_equal']:assert not m[k]
    assert len(m['profiles'])==4 and len(m['cases'])==len(report['cases'])==6 and len(report['input_readbacks'])==len(report['bias_readbacks'])==12 and len(report['comparisons'])==24
    for f,h in m['source_sha256'].items():assert sha(repo/f)==h,f
    prior=repo/m['prior_archive'];assert sha(prior/'evidence-sha256.json')==m['prior_evidence_sha256']
    if a.verify_local_inputs:
        for f,h in m['input_sha256'].items():assert sha(f)==h,f
    def capture(stem):
        e=json.loads((root/(stem+'.execution.json')).read_text());data=gzip.decompress((root/(stem+'.stdout.csv.gz')).read_bytes());err=gzip.decompress((root/(stem+'.stderr.log.gz')).read_bytes());assert e['returncode']==0 and hashlib.sha256(data).hexdigest()==e['csv_sha256'] and hashlib.sha256(err).hexdigest()==e['stderr_sha256']
        if 'console_stdout_sha256' in e:assert hashlib.sha256(gzip.decompress((root/(stem+'.console-stdout.log.gz')).read_bytes())).hexdigest()==e['console_stdout_sha256']
        return data,err
    def score(rr,refs,city,allowed):return scorer([dict(r,status='FLOAT') if r['status'] in allowed else dict(r,status='UNSCORED') for r in helper['normalize'](rr)],refs,city)
    totals={n:0 for n in m['profiles']};native_calls=0;observation_groups=0
    for case,c in m['cases'].items():
        source,_=capture(case+'-input-export');bias,_=capture(case+'-bias-export');epochs=[r for r in helper['rows'](source) if r['kind']=='epoch'];assert all(r['flag']=='0' and not r['receiver_offset_s'] for r in epochs)
        observation_groups+=sum(r['kind']=='observation' for r in helper['rows'](source))
        refs=None
        if a.verify_local_inputs:
            with Path(c['reference']).open(newline='') as f:rr=[{k.strip():v for k,v in r.items()} for r in csv.DictReader(f)]
            refs={(int(r['GPS Week']),Decimal(r['GPS TOW (s)'])):r for r in rr}
        outputs={}
        for profile,s in m['profiles'].items():
            stem=case+'-'+profile;data,err=capture(stem);actual=report['cases'][case][profile]
            if s['backend']=='rust':
                rr,metrics=baseline['inspect'](data);metrics['filter_processing']=helper['timing'](rr,'filter_ms');metrics['epoch_stage_processing']=helper['timing'](rr,'epoch_stage_ms');assert b'truth_used_by_solver=false' in err
                if profile=='rust_default':assert actual['all_previous_nontiming_CSV_bytes_identical'] and helper['project_rust'](data)==gzip.decompress((repo/'docs/comparison/2026-10-21-ppc-ppp-first'/(case+'-1C2W-raw-float.stdout.csv.gz')).read_bytes())
            else:
                rr,metrics=helper['inspect_native'](data);native_calls+=len(rr)
                assert b'whole_native_PPP_processEpoch=true' in err and b'native_RINEX_observation_pipeline_executed=false' in err and b'upstream_source_patched=false' in err and b'native_SPP_precise_loaded=0 native_SPP_DCB_loaded=0' in err
                readback=gzip.decompress((root/(stem+'.input-readback.csv.gz')).read_bytes());assert helper['readback_audit'](source,readback)==report['input_readbacks'][stem]
                bias_back=gzip.decompress((root/(stem+'.bias-readback.csv.gz')).read_bytes());assert helper['bias_audit'](bias,bias_back)==report['bias_readbacks'][stem]
                if refs is not None:
                    assert score(rr,refs,c['city'],{'SPP'})==actual['SPP_reference_errors'] and score(rr,refs,c['city'],{'SPP','PPP_FLOAT'})==actual['all_valid_reference_errors']
            assert len(rr)==len(epochs) and all(actual[k]==v for k,v in metrics.items());totals[profile]+=len(rr);outputs[profile]=helper['normalize'](rr)
            assert actual['PPP_reference_errors']['valid_positions_without_reference']==0 and not actual['runtime_memory_acceptance']
            if refs is not None:assert score(rr,refs,c['city'],{'PPP_FLOAT'})==actual['PPP_reference_errors']
            resource=json.loads((root/actual['resource_log']).read_text());assert resource['returncode']==0 and resource['peak_RSS_KiB']>0 and not resource['performance_accepted']
        if refs is not None:
            for control,current in [('native_builtin','rust_default'),('native_builtin','rust_native_head'),('native_header_zero','rust_native_head'),('native_builtin','native_header_zero')]:
                assert compare(outputs[current],outputs[control],refs,c['city'],scorer)==report['comparisons'][case+'__'+control+'__'+current]
    assert set(totals.values())=={58634} and native_calls==117268 and observation_groups==894040
    build=json.loads((root/'native-build.json').read_text());assert build['returncode']==0 and sha(root/'native-build.json')==m['native_build_sha256'] and not build['upstream_source_patched'] and sha(root/'native-build.log.gz')==build['log_sha256']
    assert build['pinned_upstream_commit']=='72f3b7c2c4c088dc575286a6f0edf4e407bc499a' and sha(repo/build['source'])==build['source_sha256'] and sha(repo/build['builder_source'])==build['builder_sha256']
    assert sha(repo/build['earlier_native_build_proof'])==build['earlier_native_build_proof_sha256']
    if a.verify_local_inputs:
        assert all(sha(f)==h for f,h in build['native_library_sha256'].items())
        assert sha(Path(m['work_dir']).parent/'native-ppp-benchmark-20261021/native_ppc_ppp')==build['binary_sha256']==m['binary_sha256']['native']
    q=json.loads((root/'qualification.json').read_text());assert q['predeclared_inputs_sha256']==sha(root/'predeclared-inputs.json') and {c['name'] for c in q['checks']}=={'examples','fmt','clippy','release'}
    for c in q['checks']:
        data=gzip.decompress((root/c['log']).read_bytes());assert c['returncode']==0 and hashlib.sha256(data).hexdigest()==c['uncompressed_sha256']
    inherited=q['unchanged_library_tests'];old=repo/inherited['archive'];assert inherited['passed_tests']==413 and sha(old/'qualification.json')==inherited['qualification_sha256'];old_m=json.loads((old/'predeclared-inputs.json').read_text());assert all(m['source_sha256'][f]==h for f,h in old_m['source_sha256'].items() if f.startswith('src/'))
    pilot=root/'pilot';pm=json.loads((pilot/'predeclared-smoke.json').read_text());audit=json.loads((pilot/'audit.json').read_text())
    assert sha(pilot/'source-before-timing.cpp')==pm['source_sha256']['scripts/native_ppc_ppp.cpp'] and sha(repo/'examples/ppc_ppp_inputs.rs')==pm['source_sha256']['examples/ppc_ppp_inputs.rs']
    source=gzip.decompress((root/'nagoya-run1-input-export.stdout.csv.gz').read_bytes());assert hashlib.sha256(source).hexdigest()==pm['input_sha256'][pm['command'][1]]
    assert (pilot/'native-builtin.returncode').read_text().strip()=='0'
    assert helper['readback_audit'](source,gzip.decompress((pilot/'native-builtin-readback.csv.gz').read_bytes()),128)==audit['input_readback']
    assert helper['bias_audit'](gzip.decompress((pilot/'bias-export.csv.gz').read_bytes()),gzip.decompress((pilot/'native-builtin-bias.csv.gz').read_bytes()))==audit['bias_readback']
    pr=helper['rows'](gzip.decompress((pilot/'native-builtin.csv.gz').read_bytes()));fr=helper['rows'](gzip.decompress((root/'nagoya-run1-native_builtin.stdout.csv.gz').read_bytes()))[:128];keys=list(pr[0]);assert len(pr)==len(fr)==128 and len(keys)==30
    assert all({k:x[k] for k in keys}==y for x,y in zip(fr,pr)) and all(r['status']=='PPP_FLOAT' for r in pr) and float(pr[-1]['filter_zenith_m'])<0
    failure=gzip.decompress((root/'initial-build-failure/failed.log.gz').read_bytes());assert b'misleading-indentation' in failure and b'uninitialized' in failure
    assert not gzip.decompress((root/'initial-build-failure/corrected-pilot.log.gz').read_bytes())
    sealed=json.loads((root/'evidence-sha256.json').read_text());assert set(sealed)=={str(f.relative_to(root)) for f in root.rglob('*') if f.is_file() and f.name!='evidence-sha256.json'}
    for f,h in sealed.items():assert sha(root/f)==h,f
    print('Verified 117268 actual native PPP processEpoch calls, 24 profiles, 12 exact decoded/bias readbacks and unchanged 413-test library; GPS-only diagnostic, parity unaccepted')


if __name__=='__main__':main()
