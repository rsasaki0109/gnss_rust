"""Whole pinned native PPP versus Rust on common decoded GPS 1C/2W PPC inputs."""
import argparse
import collections
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import runpy
import subprocess
import sys


def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def rows(data):return list(csv.DictReader(io.StringIO(data.decode())))
def normalize(rr):return [dict(r,week=str(int(r['week'])),tow=str(Decimal(r['tow']).quantize(Decimal('.000001')))) for r in rr]
def timing(rr,field):
    values=[float(r[field]) for r in rr];assert all(math.isfinite(x) and x>=0 for x in values)
    ordered=sorted(values)
    def quantile(p):
        position=(len(ordered)-1)*p;i=int(position);return ordered[i]+(ordered[min(i+1,len(ordered)-1)]-ordered[i])*(position-i)
    return dict(calls=len(values),total_ms=sum(values),p50_ms=quantile(.5),p95_ms=quantile(.95),maximum_ms=max(values),timing_acceptance=False)
def inspect_native(data):
    rr=rows(data);assert rr;counts=collections.Counter(r['status'] for r in rr);assert set(counts)<={'PPP_FLOAT','SPP','INVALID'}
    seen=set();longest=streak=0
    for r in rr:
        key=(int(r['week']),Decimal(r['tow']).quantize(Decimal('.000001')));assert key not in seen;seen.add(key)
        if r['status']!='INVALID':
            assert r['native_is_valid']=='1' and all(math.isfinite(float(r[k])) for k in ['x_m','y_m','z_m']) and not r['error']
            assert (int(r['solution_week']),Decimal(r['solution_tow']).quantize(Decimal('.000001')))==key
        else:assert all(not r[k] for k in ['x_m','y_m','z_m'])
        if r['status']=='PPP_FLOAT':streak+=1;longest=max(longest,streak)
        else:streak=0
    negative=lambda r:bool(r['filter_zenith_m']) and math.isfinite(float(r['filter_zenith_m'])) and float(r['filter_zenith_m'])<0
    return rr,dict(input_epochs=len(rr),status_counts=dict(counts),PPP_fraction=counts['PPP_FLOAT']/len(rr),all_valid_fraction=(counts['PPP_FLOAT']+counts['SPP'])/len(rr),longest_PPP_epoch_streak=longest,negative_filter_ZTD_epochs=sum(negative(r) for r in rr),PPP_negative_ZTD_epochs=sum(negative(r) for r in rr if r['status']=='PPP_FLOAT'),invalid_reasons=dict(collections.Counter(r['error'] for r in rr if r['status']=='INVALID')),reported_processing=timing(rr,'native_reported_processing_ms'),wrapper_processing=timing(rr,'wrapper_process_ms'))
def project_rust(data):
    rr=rows(data);keys=list(rr[0])[:29];out=io.StringIO();writer=csv.DictWriter(out,fieldnames=keys,lineterminator='\n');writer.writeheader();writer.writerows({k:r[k] for k in keys} for r in rr);return out.getvalue().encode()
def readback_audit(source,readback,maximum_epochs=None):
    a,b=rows(source),rows(readback)
    if maximum_epochs is not None:
        stop=0;selected=[]
        for r in a:
            if r['kind']=='epoch':
                if stop==maximum_epochs:break
                stop+=1
            selected.append(r)
        a=selected
    assert len(a)==len(b)
    numeric={'tow','code_m','phase_cycles','doppler_hz','snr_db_hz','receiver_offset_s','header_x_m','header_y_m','header_z_m'}
    for i,(x,y) in enumerate(zip(a,b)):
        assert x.keys()==y.keys()
        for k in x:
            assert bool(x[k])==bool(y[k]),(i,k,'presence mismatch')
            if k in numeric and x[k]:assert float(x[k]).hex()==float(y[k]).hex(),(i,k,x[k],y[k])
            else:assert x[k]==y[k],(i,k,x[k],y[k])
    return dict(rows=len(a),epochs=sum(r['kind']=='epoch' for r in a),observation_groups=sum(r['kind']=='observation' for r in a),all_float_bits_values_flags_and_source_options_identical=True,native_LLI_SSI_presence_retained_as_external_metadata=True)
def bias_audit(source,readback):
    aa,bb=rows(source),rows(readback);assert len(aa)==len(bb)==64
    for a,b in zip(aa,bb):
        assert (a['satellite'],a['observation'],a['unit'])==(b['satellite'],b['observation'],b['unit'])
        assert float(a['value']).hex()==float(b['value']).hex() and float(a['sigma']).hex()==float(b['sigma']).hex()
        assert abs(float(a['metres'])-float(b['metres']))<1e-10
    return dict(rows=64,numeric_values_sigmas_units_identical=True,maximum_metres_difference=max(abs(float(a['metres'])-float(b['metres'])) for a,b in zip(aa,bb)),bound_m=1e-10)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--work-dir',type=Path,required=True);p.add_argument('--native-build-dir',type=Path,required=True);a=p.parse_args();repo=Path(__file__).resolve().parents[1];out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True);work=a.work_dir.resolve();work.mkdir(parents=True,exist_ok=False);assert not (out/'predeclared-inputs.json').exists()
    earlier=repo/'docs/comparison/2026-10-21-ppp-epoch-clock';old=json.loads((earlier/'predeclared-inputs.json').read_text());cases=old['cases'];antex=next(f for f in old['input_sha256'] if f.endswith('.atx'));build_root=a.native_build_dir.resolve();build=json.loads((build_root/'build.json').read_text());native=build_root/'native_ppc_ppp';assert build['returncode']==0 and sha(native)==build['binary_sha256']
    binaries={'native':native,'rust':repo/'target/release/examples/ppc_ppp','export':repo/'target/release/examples/ppc_ppp_inputs'}
    sources=list((repo/'src').rglob('*.rs'))+[repo/f for f in ['examples/ppc_ppp.rs','examples/ppc_ppp_inputs.rs','scripts/native_ppc_ppp.cpp','scripts/build_native_ppc_ppp.py','scripts/run_native_ppc_ppp_comparison.py','scripts/run_ppc_ppp.py','scripts/run_ppc_ppp_initial_prior.py','scripts/run_ppc_rtk_float.py','scripts/run_native_ppc_comparison.py','scripts/measure_ppp_process.py','docs/ppc-native-ppp.md','Cargo.toml','Cargo.lock','rust-toolchain.toml']]
    inputs=dict(old['input_sha256']);inputs.update({str(Path(c['rover']).parent/'base.nav'):sha(Path(c['rover']).parent/'base.nav') for c in cases.values()})
    profiles=dict(rust_default=dict(backend='rust',arguments=[]),rust_native_head=dict(backend='rust',arguments=['native-head']),native_builtin=dict(backend='native',seed='builtin'),native_header_zero=dict(backend='native',seed='header-zero'))
    manifest=dict(base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),source_sha256={str(f.relative_to(repo)):sha(f) for f in sources},binary_sha256={n:sha(f) for n,f in binaries.items()},input_sha256=inputs,cases=cases,profiles=profiles,pair='1C2W',phase='raw carrier, no phase OSBs',native_build_sha256=sha(build_root/'build.json'),prior_archive=str(earlier.relative_to(repo)),prior_evidence_sha256=sha(earlier/'evidence-sha256.json'),work_dir=str(work),full_raw_rover_epochs=True,whole_native_PPP_processEpoch_executed=True,native_RINEX_observation_pipeline_executed=False,shared_decoded_input=True,native_solver_source_patched=False,truth_used_by_solver=False,practical_parity_accepted=False,native_GNSS_environment_removed=sorted(k for k in os.environ if k.startswith('GNSS_')),timing_scopes_equal=False)
    assert all(sha(f)==h for f,h in inputs.items())
    helpers=runpy.run_path(str(repo/'scripts/run_native_ppc_comparison.py'));save,retain=[helpers[k] for k in ['save','retain']];scorer=runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))['score'];baseline=runpy.run_path(str(repo/'scripts/run_ppc_ppp.py'));comparison=runpy.run_path(str(repo/'scripts/run_ppc_ppp_initial_prior.py'))['compare']
    save(out/'predeclared-inputs.json',manifest);(out/'native-build.json').write_bytes((build_root/'build.json').read_bytes());(out/'native-build.log.gz').write_bytes((build_root/'build.log.gz').read_bytes());report=dict(cases={},input_readbacks={},bias_readbacks={},comparisons={},whole_native_PPP_processEpoch_executed=True,practical_parity_accepted=False)
    def execute(stem,command,output_file=None,native_env=False):
        assert all(sha(repo/f)==h for f,h in manifest['source_sha256'].items())
        resource=out/(stem+'.resource.json');wrapper=[sys.executable,str(repo/'scripts/measure_ppp_process.py'),str(resource)];env={k:v for k,v in os.environ.items() if not (native_env and k.startswith('GNSS_'))};r=subprocess.run(wrapper+command,capture_output=True,cwd=repo,env=env)
        data=output_file.read_bytes() if output_file is not None and output_file.exists() else r.stdout
        retain(out,stem,command,data,r.stderr,r.returncode)
        if output_file is not None:
            (out/(stem+'.console-stdout.log.gz')).write_bytes(gzip.compress(r.stdout,mtime=0));e=json.loads((out/(stem+'.execution.json')).read_text());e['csv_from_explicit_output_file']=str(output_file);e['console_stdout_sha256']=hashlib.sha256(r.stdout).hexdigest();save(out/(stem+'.execution.json'),e)
        assert r.returncode==0,r.stderr.decode();return data
    def score(rr,refs,city,allowed):return scorer([dict(r,status='FLOAT') if r['status'] in allowed else dict(r,status='UNSCORED') for r in normalize(rr)],refs,city)
    for case,c in cases.items():
        ledger=work/(case+'.input.csv');export_data=execute(case+'-input-export',[str(binaries['export']),'observations',c['rover']]);ledger.write_bytes(export_data);epochs=[r for r in rows(export_data) if r['kind']=='epoch'];assert all(r['flag']=='0' and not r['receiver_offset_s'] for r in epochs)
        view=work/(case+'.native-code.bsx');audit=execute(case+'-bias-export',[str(binaries['export']),'native-bias',c['bias'],epochs[0]['week'],epochs[0]['tow'],epochs[-1]['tow'],str(view)]);(out/(case+'.native-code-view.bsx.gz')).write_bytes(gzip.compress(view.read_bytes(),mtime=0))
        with Path(c['reference']).open(newline='') as f:refs=[{k.strip():v for k,v in r.items()} for r in csv.DictReader(f)]
        references={(int(r['GPS Week']),Decimal(r['GPS TOW (s)'])):r for r in refs};outputs={}
        for profile,settings in profiles.items():
            stem=case+'-'+profile
            if settings['backend']=='rust':
                command=[str(binaries['rust'])]+[c[k] for k in ['rover','orbit','clock','bias']]+[antex,case,'1C2W','raw-float']+settings['arguments']+['-timing'];data=execute(stem,command);rr,metrics=baseline['inspect'](data);metrics['filter_processing']=timing(rr,'filter_ms');metrics['epoch_stage_processing']=timing(rr,'epoch_stage_ms')
                if profile=='rust_default':
                    before=repo/'docs/comparison/2026-10-21-ppc-ppp-first'/(case+'-1C2W-raw-float.stdout.csv.gz');assert project_rust(data)==gzip.decompress(before.read_bytes());metrics['all_previous_nontiming_CSV_bytes_identical']=True
            else:
                result=work/(stem+'.solution.csv');readback=work/(stem+'.readback.csv');bias_readback=work/(stem+'.bias-readback.csv');command=[str(binaries['native']),str(ledger),c['orbit'],c['clock'],str(view),str(Path(c['rover']).parent/'base.nav'),antex,str(result),str(readback),str(bias_readback),settings['seed']];data=execute(stem,command,result,True);rr,metrics=inspect_native(data)
                (out/(stem+'.input-readback.csv.gz')).write_bytes(gzip.compress(readback.read_bytes(),mtime=0));(out/(stem+'.bias-readback.csv.gz')).write_bytes(gzip.compress(bias_readback.read_bytes(),mtime=0));report['input_readbacks'][stem]=readback_audit(export_data,readback.read_bytes());report['bias_readbacks'][stem]=bias_audit(audit,bias_readback.read_bytes())
                metrics['SPP_reference_errors']=score(rr,references,c['city'],{'SPP'});metrics['all_valid_reference_errors']=score(rr,references,c['city'],{'SPP','PPP_FLOAT'})
            assert len(rr)==len(epochs);outputs[profile]=normalize(rr);metrics['PPP_reference_errors']=score(rr,references,c['city'],{'PPP_FLOAT'});metrics['resource_log']=stem+'.resource.json';metrics['runtime_memory_acceptance']=False;report['cases'].setdefault(case,{})[profile]=metrics;save(out/'report.json',report)
            print(case,profile,metrics['status_counts'],metrics['PPP_reference_errors']['hypotheses'].get('reference_is_antenna'),flush=True)
        for control,current in [('native_builtin','rust_default'),('native_builtin','rust_native_head'),('native_header_zero','rust_native_head'),('native_builtin','native_header_zero')]:
            key=case+'__'+control+'__'+current;report['comparisons'][key]=comparison(outputs[current],outputs[control],references,c['city'],scorer)
        save(out/'report.json',report)
    assert all(sha(f)==h for f,h in inputs.items()) and all(sha(f)==manifest['binary_sha256'][n] for n,f in binaries.items());print('Completed 24 whole-solver profiles with 12 native input/bias readbacks; no practical parity acceptance',flush=True)


if __name__=='__main__':main()
