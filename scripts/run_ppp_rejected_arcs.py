"""Freeze rejected-update phase retention on both full-rate PPC GPS IF pairs."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import json
from pathlib import Path
import runpy
import subprocess
import sys


def sha(p):
    with Path(p).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def retention_metrics(rr, enabled):
    if not enabled:
        assert all('retained_arcs' not in r for r in rr)
        return dict(enabled=False)
    counts=[]
    for r in rr:
        count=int(r['retained_arcs']);satellites=r['retained_satellites'].split(';') if count else []
        assert len(satellites)==len(set(satellites))==count
        if count:
            assert r['status']=='INVALID' and r['committed_week']
        counts.append(count)
    return dict(enabled=True,rejected_epochs_with_retained_arcs=sum(c>0 for c in counts),total_retained_arc_epochs=sum(counts),maximum_retained_arcs=max(counts),rejected_epochs_relabelled_as_solutions=False)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
    repo=Path(__file__).resolve().parents[1];out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True);assert not (out/'predeclared-inputs.json').exists()
    prior=repo/'docs/comparison/2026-10-21-ppp-epoch-clock';old=json.loads((prior/'predeclared-inputs.json').read_text());cases=old['cases'];antex=next(f for f in old['input_sha256'] if f.endswith('.atx'));binary=repo/'target/release/examples/ppc_ppp'
    helpers=runpy.run_path(str(repo/'scripts/run_native_ppc_comparison.py'));save,retain=[helpers[k] for k in ['save','retain']]
    inspect=runpy.run_path(str(repo/'scripts/run_ppc_ppp.py'))['inspect'];initial=runpy.run_path(str(repo/'scripts/run_ppc_ppp_initial_prior.py'));scorer=runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))['score']
    sources=list((repo/'src').rglob('*.rs'))+[repo/f for f in ['examples/ppc_ppp.rs','tests/phase_arc.rs','tests/ppp.rs','tests/ppp_kinematic.rs','scripts/run_ppp_rejected_arcs.py','scripts/run_ppc_ppp.py','scripts/run_ppc_ppp_initial_prior.py','scripts/run_ppc_rtk_float.py','scripts/run_native_ppc_comparison.py','scripts/measure_ppp_process.py','docs/ppp-rejected-arcs.md','Cargo.toml','Cargo.lock','rust-toolchain.toml']]
    profiles={head+'-'+clock+'-'+policy:dict(head=head,clock=clock,retain=policy=='retain',arguments=([] if head=='default' else ['native-head'])+([] if clock=='random-walk' else ['-epoch-clock'])+([] if policy=='interrupt' else ['-retain-arcs'])) for head in ['default','native-head'] for clock in ['random-walk','independent'] for policy in ['interrupt','retain']}
    native=repo/'docs/comparison/2026-10-21-native-ppp-common-inputs'
    manifest=dict(base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),source_sha256={str(f.relative_to(repo)):sha(f) for f in sources},binary_sha256=sha(binary),input_sha256=old['input_sha256'],cases=cases,pairs=old['pairs'],profiles=profiles,prior_archive=str(prior.relative_to(repo)),prior_evidence_sha256=sha(prior/'evidence-sha256.json'),native_archive=str(native.relative_to(repo)),native_evidence_sha256=sha(native/'evidence-sha256.json'),full_raw_rover_epochs=True,thresholds_changed=False,default_policy_changed=False,ambiguity_prior_changed=False,SPP_or_Doppler_fused=False,truth_used_by_solver=False,whole_native_PPP_rerun=False,practical_parity_accepted=False,execution_policy='sequential; runtime/RSS descriptive')
    assert all(sha(f)==h for f,h in manifest['input_sha256'].items());save(out/'predeclared-inputs.json',manifest);report=dict(cases={},comparisons={},practical_parity_accepted=False)
    for case,c in cases.items():
        with Path(c['reference']).open(newline='') as f:refs=[{k.strip():v for k,v in r.items()} for r in csv.DictReader(f)]
        references={(int(r['GPS Week']),Decimal(r['GPS TOW (s)'])):r for r in refs};assert len(references)==len(refs)
        for pair in manifest['pairs']:
            outputs={}
            for profile,s in profiles.items():
                assert all(sha(repo/f)==h for f,h in manifest['source_sha256'].items()) and sha(binary)==manifest['binary_sha256']
                stem=case+'-'+pair+'-'+profile;resource=out/(stem+'.resource.json');command=[str(binary)]+[c[k] for k in ['rover','orbit','clock','bias']]+[antex,case,pair,'calibrated']+s['arguments']
                r=subprocess.run([sys.executable,str(repo/'scripts/measure_ppp_process.py'),str(resource)]+command,capture_output=True,cwd=repo);retain(out,stem,command,r.stdout,r.stderr,r.returncode);assert r.returncode==0,r.stderr.decode()
                rr,metrics=inspect(r.stdout);outputs[profile]=rr;metrics['retention']=retention_metrics(rr,s['retain']);metrics['reference_errors']=scorer(initial['score_view'](rr),references,c['city']);metrics['resource_log']=resource.name;metrics['runtime_memory_acceptance']=False
                if not s['retain']:
                    before=prior/(case+'-'+pair+'-'+s['head']+'-'+s['clock']+'.stdout.csv.gz');assert r.stdout==gzip.decompress(before.read_bytes()),stem;metrics['all_prior_CSV_bytes_identical']=True
                else:
                    control=outputs[s['head']+'-'+s['clock']+'-interrupt'];report['comparisons'][stem]=initial['compare'](rr,control,references,c['city'],scorer)
                report['cases'].setdefault(case,{})[pair+'-'+profile]=metrics;save(out/'report.json',report)
                print(case,pair,profile,metrics['status_counts'],metrics['retention'],metrics['reference_errors']['hypotheses'].get('reference_is_antenna'),flush=True)
    assert all(sha(f)==h for f,h in manifest['input_sha256'].items());print('Completed 96 full-rate profiles, 48 exact controls; practical parity unaccepted',flush=True)


if __name__=='__main__':main()
