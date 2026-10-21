"""Read-only verification of raw FLOAT gauges, controls and product admission."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import runpy


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--verify-local-inputs',action='store_true');args=parser.parse_args()
    root=Path(__file__).resolve().parent;repo=root.parents[2]
    m=json.loads((root/'predeclared-inputs.json').read_text());report=json.loads((root/'report.json').read_text())
    for key in ['default_mode_changed','phase_OSB_zero_fallback','thresholds_changed','whole_native_PPP_executed','real_PPC_PPP_executed','practical_parity_accepted']:assert not m[key]
    assert not report['real_PPC_PPP_executed'] and not report['practical_parity_accepted']
    for name,h in m['source_sha256'].items():assert sha(repo/name)==h,name
    prior=repo/m['prior_archive'];assert sha(prior/'evidence-sha256.json')==m['prior_evidence_sha256']
    helper=runpy.run_path(str(repo/'scripts/run_ppp_raw_phase.py'));rows=helper['rows']
    def capture(name):
        e=json.loads((root/(name+'.execution.json')).read_text());data=gzip.decompress((root/(name+'.stdout.csv.gz')).read_bytes());err=gzip.decompress((root/(name+'.stderr.log.gz')).read_bytes())
        assert e['returncode']==0 and hashlib.sha256(data).hexdigest()==e['csv_sha256'] and hashlib.sha256(err).hexdigest()==e['stderr_sha256'];return data
    expected={r['satellite']:float(r['unmodelled_phase_IF_bias_m']) for r in rows((repo/'tests/fixtures/upstream_ppp_raw_phase.csv').read_bytes())}
    assert len(report['gauge_profiles'])==4 and len(report['default_controls'])==21
    for variant in ['clean','noisy']:
        calibrated=rows(capture('calibrated-'+variant))
        for mode in ['raw-phase','raw-dsb']:
            assert helper['compare'](calibrated,rows(capture(mode+'-'+variant)),expected)==report['gauge_profiles'][mode+'-'+variant]
    for name,c in report['default_controls'].items():
        data=capture('control-'+name);assert data==gzip.decompress((prior/c['prior_capture']).read_bytes())
        assert c['all_CSV_bytes_identical'] and hashlib.sha256(data).hexdigest()==c['csv_sha256']
    admission=rows(capture('code-products'));assert admission==report['code_product_admission'] and len(admission)==9
    assert all(r['frame']=='IGS20' and r['GPS_IF_C1W_C2W_declared']==r['bias_APC_IGS20']=='true' and int(r['GPS_sample_states_with_clock'])==32 and int(r['GPS_sample_states_unavailable'])==int(r['gauge_pairs_missing'])==0 and float(r['max_abs_GPS_IF_code_gauge_m'])<1e-4 for r in admission)
    for stamp in ['20242020000','20242050000','20242160000']:
        selected=rows(capture('antex-'+stamp));assert len(selected)==32 and all(r['status']=='AVAILABLE' for r in selected)
    assert len(report['retained_SVN_conflicts'])==3 and all(r['family']=='CAS' and r['satellite']=='G01' and r['ANTEX_SVN']=='G049' and r['bias_SVNs']==['G063'] for r in report['retained_SVN_conflicts'])
    selection=json.loads((root/'antex-selection.json').read_text());assert len(selection['selected_blocks'])==116 and len(selection['excluded_blocks'])==810 and selection['all_selected_block_bytes_unchanged'] and not selection['full_ANTEX_format_support_claimed']
    assert hashlib.sha256(gzip.decompress((root/'gps-satellite-igs20.atx.gz').read_bytes())).hexdigest()==selection['output_sha256']
    for name,item in m['acquisition'].items():assert sha(root/(name+'-acquisition.json'))==item['acquisition_sha256']
    code=json.loads((root/'code-acquisition.json').read_text());assert code['acquired_count']==27 and all(i['acquired'] for i in code['products'])
    if args.verify_local_inputs:
        for i in code['products']:
            for name,key in [('filename','compressed_sha256'),('expanded_filename','expanded_sha256')]:assert sha(Path(m['acquisition']['code']['root'])/i[name])==i[key]
        assert sha(Path(selection['source']))==selection['source_sha256']
    native=json.loads((root/'native-generator.json').read_text());assert native['returncode']==0 and not native['whole_native_PPP_executed'] and native['generator_sha256']==sha(repo/'scripts/upstream_ppp_raw_phase.cpp')
    assert len(native['repeat_runs'])==2
    for r in native['repeat_runs']:
        assert r['returncode']==0
        for name,h in r['fixture_sha256'].items():assert sha(repo/'tests/fixtures'/name)==h
    q=json.loads((root/'qualification.json').read_text());assert q['predeclared_inputs_sha256']==sha(root/'predeclared-inputs.json') and {c['name'] for c in q['checks']}=={'rust','examples','fmt','clippy','release'}
    for c in q['checks']:
        data=gzip.decompress((root/c['log']).read_bytes());assert c['returncode']==0 and hashlib.sha256(data).hexdigest()==c['uncompressed_sha256']
        if c['name']=='rust':assert sum(map(int,re.findall(rb'test result: ok\. (\d+) passed;',data)))==c['passed_tests']==406
    seal=json.loads((root/'evidence-sha256.json').read_text());assert set(seal)=={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p.name!='evidence-sha256.json'}
    for name,h in seal.items():assert sha(root/name)==h,name
    print('Verified 4 raw FLOAT gauge comparisons, 21 byte-identical controls, 9 CODE days, 96 ANTEX admissions and 406 tests; real PPC PPP still pending')


if __name__=='__main__':main()
