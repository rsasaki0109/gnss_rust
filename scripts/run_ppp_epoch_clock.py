"""Freeze independent epoch-clock stress tests and full-rate PPC clock comparisons."""
import argparse
from decimal import Decimal
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import runpy
import subprocess
import sys
import tempfile


def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def stress_compare(control,jumped):
    assert len(control)==len(jumped)==64
    position=clock=0.0
    for i,(a,b) in enumerate(zip(control,jumped)):
        assert a['status']==b['status']=='PPP_FLOAT'
        offset=0.0 if i<8 else [30000.0,-20000.0,50000.0,0.0][(i-8)%4]
        position=max(position,math.sqrt(sum((float(a[k])-float(b[k]))**2 for k in ['x_m','y_m','z_m'])))
        clock=max(clock,abs((float(b['clock_s'])-float(a['clock_s']))*299792458-offset))
    assert position<.001 and clock<.001
    return dict(epochs=64,maximum_position_difference_m=position,maximum_clock_offset_error_m=clock,bound_m=.001)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args();repo=Path(__file__).resolve().parents[1];out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True);assert not (out/'predeclared-inputs.json').exists()
    prior=repo/'docs/comparison/2026-10-21-ppp-initial-prior';old=json.loads((prior/'predeclared-inputs.json').read_text());cases=old['cases'];antex=next(f for f in old['input_sha256'] if f.endswith('.atx'))
    baseline=runpy.run_path(str(repo/'scripts/run_ppc_ppp.py'));inspect=baseline['inspect'];rows=baseline['rows'];scorer=runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))['score'];initial=runpy.run_path(str(repo/'scripts/run_ppc_ppp_initial_prior.py'));save,retain=[runpy.run_path(str(repo/'scripts/run_native_ppc_comparison.py'))[k] for k in ['save','retain']]
    moving=runpy.run_path(str(repo/'scripts/run_ppp_kinematic.py'));truth=moving['rows']((repo/'tests/fixtures/synthetic_ppp_kinematic_truth.csv').read_bytes())
    binaries={name:repo/'target/release/examples'/name for name in ['ppc_ppp','ppp_kinematic']}
    source=list((repo/'src').rglob('*.rs'))+[repo/f for f in ['examples/ppc_ppp.rs','examples/ppp_kinematic.rs','tests/ppp_kinematic.rs','tests/fixtures/upstream_ppp_epoch_clock.csv','scripts/upstream_ppp_epoch_clock.cpp','scripts/run_ppp_epoch_clock.py','scripts/run_ppc_ppp_initial_prior.py','scripts/run_ppc_ppp.py','scripts/run_ppp_kinematic.py','scripts/run_ppc_rtk_float.py','scripts/run_native_ppc_comparison.py','scripts/measure_ppp_process.py','docs/ppp-epoch-clock.md','Cargo.toml','Cargo.lock','rust-toolchain.toml']]+list((repo/'tests/fixtures').glob('*ppp_kinematic*'))
    profiles={head+'-'+clock:dict(head_prior=head,clock_model=clock,arguments=([] if head=='default' else ['native-head'])+([] if clock=='random-walk' else ['-epoch-clock']),clock_variance_m2=1e8) for head in ['default','native-head'] for clock in ['random-walk','independent']}
    manifest=dict(base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),source_sha256={str(f.relative_to(repo)):sha(f) for f in source},binary_sha256={n:sha(f) for n,f in binaries.items()},input_sha256=old['input_sha256'],cases=cases,pairs=old['pairs'],profiles=profiles,prior_archive=str(prior.relative_to(repo)),prior_evidence_sha256=sha(prior/'evidence-sha256.json'),full_raw_rover_epochs=True,thresholds_changed=False,ambiguity_prior_changed=False,SPP_or_Doppler_fused=False,truth_used_by_solver=False,default_clock_model_changed=False,whole_native_PPP_executed=False,practical_parity_accepted=False)
    assert all(sha(f)==h for f,h in manifest['input_sha256'].items());save(out/'predeclared-inputs.json',manifest)
    report=dict(synthetic={},stress_comparisons={},cases={},clock_comparisons={},whole_native_PPP_executed=False,practical_parity_accepted=False)
    def execute(stem,name,args,measure=False):
        assert all(sha(repo/f)==h for f,h in manifest['source_sha256'].items()) and sha(binaries[name])==manifest['binary_sha256'][name]
        command=[str(binaries[name])]+args;wrapper=[sys.executable,str(repo/'scripts/measure_ppp_process.py'),str(out/(stem+'.resource.json'))] if measure else []
        r=subprocess.run(wrapper+command,capture_output=True,cwd=repo);retain(out,stem,command,r.stdout,r.stderr,r.returncode);assert r.returncode==0,r.stderr.decode();return r.stdout
    for variant in ['clean','noisy']:
        independent={}
        for jumped in [False,True]:
            kind=('jumped-' if jumped else '')+variant
            for clock in ['random-walk','independent']:
                stem='synthetic-'+kind+'-'+clock;data=execute(stem,'ppp_kinematic',['kinematic',kind]+([] if clock=='random-walk' else ['-epoch-clock']));rr=rows(data);report['synthetic'][stem]=moving['moving_metrics'](data,truth,True)
                if clock=='independent':independent[jumped]=rr
                elif jumped:assert report['synthetic'][stem]['status_counts'].get('INVALID',0)>10
        report['stress_comparisons'][variant]=stress_compare(independent[False],independent[True]);save(out/'report.json',report)
    for case,c in cases.items():
        with Path(c['reference']).open(newline='') as f:refs=[{k.strip():v for k,v in r.items()} for r in csv.DictReader(f)]
        references={(int(r['GPS Week']),Decimal(r['GPS TOW (s)'])):r for r in refs};assert len(references)==len(refs)
        for pair in manifest['pairs']:
            outputs={}
            for profile,settings in profiles.items():
                stem=case+'-'+pair+'-'+profile;args=[c[k] for k in ['rover','orbit','clock','bias']]+[antex,case,pair,'calibrated']+settings['arguments'];data=execute(stem,'ppc_ppp',args,True);rr,metrics=inspect(data);outputs[profile]=rr
                metrics['reference_errors']=scorer(initial['score_view'](rr),references,c['city']);metrics['resource_log']=stem+'.resource.json';metrics['runtime_memory_acceptance']=False
                if settings['clock_model']=='random-walk':
                    capture=case+'-'+pair+'-'+settings['head_prior']+'.stdout.csv.gz';assert data==gzip.decompress((prior/capture).read_bytes()),stem;metrics['all_prior_CSV_bytes_identical']=True
                else:report['clock_comparisons'][stem]=initial['compare'](rr,outputs[settings['head_prior']+'-random-walk'],references,c['city'],scorer)
                report['cases'].setdefault(case,{})[pair+'-'+profile]=metrics;save(out/'report.json',report);print(case,pair,profile,metrics['status_counts'],metrics['reference_errors']['hypotheses'].get('reference_is_antenna'),flush=True)
    native_prior=json.loads((repo/'docs/comparison/2026-10-21-ppc-ppp-first/native-initial.json').read_text());native=dict(upstream_commit=native_prior['upstream_commit'],source_archive_sha256=native_prior['source_archive_sha256'],native_library_sha256=native_prior['native_library_sha256'],actual_native_scalar_reset_calls_per_repeat=16,whole_native_PPP_executed=False,upstream_source_patched=False,repeat_runs=[])
    with tempfile.TemporaryDirectory(prefix='ppp-epoch-clock-',dir='/workspace/.cloud-setup') as temp:
        temp=Path(temp);binary=temp/'native';source=repo/'scripts/upstream_ppp_epoch_clock.cpp';command=['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror','-I/workspace/.cloud-setup/reference/gnssplusplus-library/include','-I/workspace/.cloud-setup/eigen-3.4.0',str(source),'-Wl,--start-group']+list(native['native_library_sha256'])+['-Wl,--end-group','-o',str(binary)]
        assert all(sha(f)==h for f,h in native['native_library_sha256'].items());r=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(out/'native-build.log.gz').write_bytes(gzip.compress(r.stdout,mtime=0));assert r.returncode==0,r.stdout.decode();native.update(command=command,build_returncode=0,source_sha256=sha(source),binary_sha256=sha(binary),compiler=subprocess.check_output(['g++','--version'],text=True).splitlines()[0])
        for _ in range(2):
            r=subprocess.run([str(binary),str(temp)],capture_output=True);assert r.returncode==0 and (temp/'upstream_ppp_epoch_clock.csv').read_bytes()==(repo/'tests/fixtures/upstream_ppp_epoch_clock.csv').read_bytes();native['repeat_runs'].append(dict(returncode=0,fixture_sha256=sha(temp/'upstream_ppp_epoch_clock.csv')))
    save(out/'native-generator.json',native);assert all(sha(f)==h for f,h in manifest['input_sha256'].items());print('Completed 8 synthetic clock profiles, 48 full-rate PPC profiles and 24 exact controls',flush=True)


if __name__=='__main__':main()
