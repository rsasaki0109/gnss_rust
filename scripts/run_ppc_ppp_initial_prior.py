"""Freeze the 2x2 initial-head prior experiment on both PPC GPS IF pairs."""
import argparse
from decimal import Decimal
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
def score_view(rows):return [dict(r,tow=str(Decimal(r['tow']).quantize(Decimal('.000001'))),status='FLOAT') if r['status']=='PPP_FLOAT' else r for r in rows]
def compare(current,control,references,city,scorer):
    aa,bb=[{(r['week'],r['tow']):r for r in rr if r['status']=='PPP_FLOAT'} for rr in [control,current]];common=aa.keys()&bb.keys();added=bb.keys()-aa.keys();lost=aa.keys()-bb.keys()
    deltas=[math.sqrt(sum((float(aa[t][k])-float(bb[t][k]))**2 for k in ['x_m','y_m','z_m'])) for t in sorted(common)]
    return dict(common_epochs=len(common),newly_accepted_epochs=len(added),lost_epochs=len(lost),common_position_difference_RMS_m=math.sqrt(sum(v*v for v in deltas)/len(deltas)) if deltas else None,common_position_difference_maximum_m=max(deltas) if deltas else None,common_control_errors=scorer(score_view([aa[t] for t in sorted(common)]),references,city),common_current_errors=scorer(score_view([bb[t] for t in sorted(common)]),references,city),newly_accepted_errors=scorer(score_view([bb[t] for t in sorted(added)]),references,city),lost_control_errors=scorer(score_view([aa[t] for t in sorted(lost)]),references,city))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
    repo=Path(__file__).resolve().parents[1];out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True);assert not (out/'predeclared-inputs.json').exists()
    prior=repo/'docs/comparison/2026-10-21-ppc-ppp-first';old=json.loads((prior/'predeclared-inputs.json').read_text());cases=old['cases'];binary=repo/'target/release/examples/ppc_ppp';antex=next(f for f in old['input_sha256'] if f.endswith('.atx'))
    native_helper=runpy.run_path(str(repo/'scripts/run_native_ppc_comparison.py'));save,retain=[native_helper[k] for k in ['save','retain']]
    baseline=runpy.run_path(str(repo/'scripts/run_ppc_ppp.py'));inspect=baseline['inspect'];scorer=runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))['score']
    sources=list((repo/'src').rglob('*.rs'))+[repo/f for f in ['Cargo.toml','Cargo.lock','rust-toolchain.toml','examples/ppc_ppp.rs','scripts/run_ppc_ppp_initial_prior.py','scripts/run_ppc_ppp.py','scripts/measure_ppp_process.py','scripts/run_ppc_rtk_float.py','scripts/run_native_ppc_comparison.py','scripts/native_ppc_ppp_initial.cpp','scripts/upstream_zenith_climatology.cpp','tests/climatology.rs','tests/ppp_kinematic.rs','tests/fixtures/upstream_zenith_climatology.csv','tests/fixtures/upstream_ppp_initial.csv','docs/ppp-initial-prior.md']]
    profiles={n:dict(native_head_variances=n in ['native-variances','native-head'],climatological_initial_ZTD=n in ['climatology','native-head'],arguments=[] if n=='default' else [n]) for n in ['default','native-variances','climatology','native-head']}
    manifest=dict(base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),source_sha256={str(f.relative_to(repo)):sha(f) for f in sources},binary_sha256=sha(binary),input_sha256=old['input_sha256'],prior_archive=str(prior.relative_to(repo)),prior_evidence_sha256=sha(prior/'evidence-sha256.json'),profiles=profiles,cases=cases,pairs=['1C2W','1C2L'],full_raw_rover_epochs=True,phase='calibrated CODE OSB',thresholds_changed=False,process_noise_changed=False,ambiguity_prior_changed=False,truth_used_by_solver=False,whole_native_PPP_executed=False,practical_parity_accepted=False)
    assert all(sha(f)==h for f,h in manifest['input_sha256'].items());save(out/'predeclared-inputs.json',manifest)
    report=dict(cases={},comparisons={},whole_native_PPP_executed=False,practical_parity_accepted=False)
    for case,c in cases.items():
        import csv
        with Path(c['reference']).open(newline='') as f:refs=[{k.strip():v for k,v in r.items()} for r in csv.DictReader(f)]
        references={(int(r['GPS Week']),Decimal(r['GPS TOW (s)'])):r for r in refs};assert len(references)==len(refs)
        for pair in manifest['pairs']:
            outputs={}
            for profile,settings in profiles.items():
                assert all(sha(repo/f)==h for f,h in manifest['source_sha256'].items()) and sha(binary)==manifest['binary_sha256']
                stem=case+'-'+pair+'-'+profile;resources=out/(stem+'.resource.json');command=[str(binary)]+[c[k] for k in ['rover','orbit','clock','bias']]+[antex,case,pair,'calibrated']+settings['arguments']
                r=subprocess.run([sys.executable,str(repo/'scripts/measure_ppp_process.py'),str(resources)]+command,capture_output=True,cwd=repo);retain(out,stem,command,r.stdout,r.stderr,r.returncode);assert r.returncode==0,r.stderr.decode()
                rr,metrics=inspect(r.stdout);outputs[profile]=rr
                metrics['reference_errors']=scorer(score_view(rr),references,c['city']);metrics['resource_log']=resources.name;metrics['runtime_memory_acceptance']=False
                if profile=='default':
                    previous=gzip.decompress((prior/(case+'-'+pair+'-calibrated.stdout.csv.gz')).read_bytes());assert previous==r.stdout,(case,pair,'changed default CSV');metrics['all_prior_CSV_bytes_identical']=True
                else:report['comparisons'][stem]=compare(rr,outputs['default'],references,c['city'],scorer)
                report['cases'].setdefault(case,{})[pair+'-'+profile]=metrics;save(out/'report.json',report)
                print(case,pair,profile,metrics['status_counts'],metrics['reference_errors']['hypotheses'].get('reference_is_antenna'),flush=True)
    # Independent component fixtures, both recreated twice using actual native code.
    native_prior=json.loads((prior/'native-initial.json').read_text());upstream=Path('/workspace/.cloud-setup/reference/gnssplusplus-library');native=dict(upstream_commit=native_prior['upstream_commit'],source_archive_sha256=native_prior['source_archive_sha256'],whole_native_PPP_executed=False,upstream_source_patched=False,builds=[],repeat_runs=[])
    with tempfile.TemporaryDirectory(prefix='ppp-initial-fixtures-',dir='/workspace/.cloud-setup') as temp:
        temp=Path(temp)
        for label,source in [('climatology','scripts/upstream_zenith_climatology.cpp'),('initial','scripts/native_ppc_ppp_initial.cpp')]:
            binary=temp/label;cmd=['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror','-I'+str(upstream/'include'),'-I/workspace/.cloud-setup/eigen-3.4.0',str(repo/source)]
            if label=='initial':cmd+=['-Wl,--start-group']+list(native_prior['native_library_sha256'])+['-Wl,--end-group'];assert all(sha(f)==h for f,h in native_prior['native_library_sha256'].items())
            cmd+=['-o',str(binary)];r=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);log='native-'+label+'-build.log.gz';(out/log).write_bytes(gzip.compress(r.stdout,mtime=0));assert r.returncode==0,r.stdout.decode();native['builds'].append(dict(command=cmd,returncode=0,log=log,log_sha256=sha(out/log),binary_sha256=sha(binary),source=source,source_sha256=sha(repo/source)))
            for repeat in range(2):
                command=[str(binary),str(temp) if label=='climatology' else '/workspace/.cloud-setup/datasets/PPC-Dataset'];r=subprocess.run(command,capture_output=True);assert r.returncode==0
                data=(temp/'upstream_zenith_climatology.csv').read_bytes() if label=='climatology' else r.stdout;name='upstream_zenith_climatology.csv' if label=='climatology' else 'upstream_ppp_initial.csv';assert data==(repo/'tests/fixtures'/name).read_bytes()
                native['repeat_runs'].append(dict(kind=label,command=command,returncode=0,fixture=name,fixture_sha256=sha(repo/'tests/fixtures'/name),actual_native_component_calls=312 if label=='climatology' else 6))
    save(out/'native-fixtures.json',native);assert all(sha(f)==h for f,h in manifest['input_sha256'].items())
    print('Completed 48 full-rate initialization profiles, 12 byte-identical controls and repeated native fixtures',flush=True)


if __name__=='__main__':main()
