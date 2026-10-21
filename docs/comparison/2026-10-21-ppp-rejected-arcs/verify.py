"""Recheck full-rate rejected-update retention evidence; optionally rescore local truth."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import json
from pathlib import Path
import re
import runpy


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-local-inputs',action='store_true')
    args=parser.parse_args();root=Path(__file__).resolve().parent;repo=root.parents[2]
    m=json.loads((root/'predeclared-inputs.json').read_text());report=json.loads((root/'report.json').read_text())
    inspect=runpy.run_path(str(repo/'scripts/run_ppc_ppp.py'))['inspect']
    metrics=runpy.run_path(str(repo/'scripts/run_ppp_rejected_arcs.py'))['retention_metrics']
    initial=runpy.run_path(str(repo/'scripts/run_ppc_ppp_initial_prior.py'))
    scorer=runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))['score']
    assert m['full_raw_rover_epochs'] and len(m['profiles'])==8 and len(m['cases'])==6
    for k in ['thresholds_changed','default_policy_changed','ambiguity_prior_changed','SPP_or_Doppler_fused','truth_used_by_solver','whole_native_PPP_rerun','practical_parity_accepted']:
        assert not m[k]
    for f,h in m['source_sha256'].items():assert sha(repo/f)==h,f
    prior=repo/m['prior_archive'];assert sha(prior/'evidence-sha256.json')==m['prior_evidence_sha256']
    assert sha(repo/m['native_archive']/'evidence-sha256.json')==m['native_evidence_sha256']
    if args.verify_local_inputs:
        for f,h in m['input_sha256'].items():assert sha(f)==h,f
    total=controls=profiles=0
    for case,c in m['cases'].items():
        references=None
        if args.verify_local_inputs:
            with Path(c['reference']).open(newline='') as f:refs=[{k.strip():v for k,v in r.items()} for r in csv.DictReader(f)]
            references={(int(r['GPS Week']),Decimal(r['GPS TOW (s)'])):r for r in refs};assert len(references)==len(refs)
        for pair in m['pairs']:
            outputs={}
            for profile,s in m['profiles'].items():
                stem=case+'-'+pair+'-'+profile;execution=json.loads((root/(stem+'.execution.json')).read_text())
                data=gzip.decompress((root/(stem+'.stdout.csv.gz')).read_bytes());err=gzip.decompress((root/(stem+'.stderr.log.gz')).read_bytes())
                assert execution['returncode']==0 and hashlib.sha256(data).hexdigest()==execution['csv_sha256'] and hashlib.sha256(err).hexdigest()==execution['stderr_sha256']
                assert b'truth_used_by_solver=false' in err
                rr,actual=inspect(data);outputs[profile]=rr;summary=report['cases'][case][pair+'-'+profile]
                assert all(summary[k]==v for k,v in actual.items())
                assert metrics(rr,s['retain'])==summary['retention'] and not summary['runtime_memory_acceptance']
                if s['retain']:
                    assert b'rejected_update_policy=RetainContinuousArcs rejected_measurements_accepted=false thresholds_changed=false' in err
                    if references is not None:
                        assert initial['compare'](rr,outputs[s['head']+'-'+s['clock']+'-interrupt'],references,c['city'],scorer)==report['comparisons'][stem]
                else:
                    assert summary['all_prior_CSV_bytes_identical']
                    previous=prior/(case+'-'+pair+'-'+s['head']+'-'+s['clock']+'.stdout.csv.gz')
                    assert data==gzip.decompress(previous.read_bytes());controls+=1
                if references is not None:assert scorer(initial['score_view'](rr),references,c['city'])==summary['reference_errors']
                resource=json.loads((root/summary['resource_log']).read_text());assert resource['returncode']==0 and resource['peak_RSS_KiB']>0 and not resource['performance_accepted']
                total+=len(rr);profiles+=1
    assert total==938144 and profiles==96 and controls==48 and len(report['comparisons'])==48 and not report['practical_parity_accepted']
    qualification=json.loads((root/'qualification.json').read_text());assert qualification['predeclared_inputs_sha256']==sha(root/'predeclared-inputs.json')
    assert {c['name'] for c in qualification['checks']}=={'tests','fmt','clippy','release','examples'}
    for c in qualification['checks']:
        data=gzip.decompress((root/c['log']).read_bytes());assert c['returncode']==0 and hashlib.sha256(data).hexdigest()==c['uncompressed_sha256']
        if c['name']=='tests':assert sum(map(int,re.findall(rb'test result: ok\. (\d+) passed;',data)))==420
    outliers=json.loads((root/'worst-points.json').read_text())
    assert outliers['post_hoc_diagnostic_only'] and not outliers['used_to_select_solver_profile'] and not outliers['thresholds_changed']
    assert outliers['source_sha256']==sha(repo/'scripts/analyze_ppp_rejected_arcs.py') and outliers['predeclared_inputs_sha256']==sha(root/'predeclared-inputs.json')
    if args.verify_local_inputs:
        assert runpy.run_path(str(repo/'scripts/analyze_ppp_rejected_arcs.py'))['derive'](root,repo)==outliers
    seal=json.loads((root/'evidence-sha256.json').read_text());assert set(seal)=={str(f.relative_to(root)) for f in root.rglob('*') if f.is_file() and f.name!='evidence-sha256.json'}
    for f,h in seal.items():assert sha(root/f)==h,f
    print('Verified 96 profiles / 938144 rows, 48 byte-identical controls and 420 tests; practical parity unaccepted')


if __name__=='__main__':main()
