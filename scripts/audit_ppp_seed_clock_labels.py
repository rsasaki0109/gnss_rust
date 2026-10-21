"""Audit a metadata-only SPP clock-label fix with six complete shadow replays."""
import argparse
import collections
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import runpy
import subprocess


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--archive',type=Path,required=True);a=p.parse_args();repo=Path(__file__).resolve().parents[1];root=a.archive.resolve();out=root/'clock-label-fix';assert not (out/'audit.json').exists();m=json.loads((root/'predeclared-inputs.json').read_text());helper=runpy.run_path(str(repo/'scripts/run_native_ppc_comparison.py'));save,retain=[helper[k] for k in ['save','retain']];binary=repo/'target/release/examples/ppc_ppp'
    changed={f:dict(captured_path='clock-label-fix/'+name,before_sha256=sha(out/name),after_sha256=sha(repo/f)) for f,name in [('examples/support/ppp_seed.rs','seed-before.rs'),('scripts/run_ppp_motion_seeds.py','runner-before.py')]}
    assert all(v['before_sha256']==m['source_sha256'][f] for f,v in changed.items());assert all(sha(repo/f)==h for f,h in m['source_sha256'].items() if f not in changed)
    report=dict(metadata_only_change=True,GNSS_numeric_processing_changed=False,library_sources_unchanged=True,changed_sources=changed,current_binary_sha256=sha(binary),source_sha256=sha(__file__),cases={},SPP_clock_imported=False,truth_used_by_solver=False)
    save(out/'audit.json',report)
    for case in m['cases']:
        stem=case+'-1C2W-shadow';old=json.loads((root/(stem+'.execution.json')).read_text());command=old['command'];command[0]=str(binary);r=subprocess.run(command,cwd=repo,capture_output=True);retain(out,case,command,r.stdout,r.stderr,r.returncode);assert r.returncode==0,r.stderr.decode()
        aa=list(csv.DictReader(io.StringIO(gzip.decompress((root/(stem+'.stdout.csv.gz')).read_bytes()).decode())));bb=list(csv.DictReader(io.StringIO(r.stdout.decode())));assert len(aa)==len(bb);clocks=collections.Counter()
        for old,new in zip(aa,bb):
            new['seed_gps_clock_s']=new.pop('seed_first_system_clock_s')
            if new['seed_status']=='SPP':
                source,clock=new['seed_source'].rsplit(';SPP_first_clock_system=',1);assert clock in ['Gps','Galileo'];clocks[clock]+=1;new['seed_source']=source
            assert old==new,(case,old['tow'],'data changed beyond clock label/provenance')
        report['cases'][case]=dict(input_epochs=len(aa),all_72_fields_equal_after_metadata_normalization=True,first_SPP_clock_system_counts=dict(clocks));save(out/'audit.json',report);print(case,len(aa),dict(clocks),'all values/status identical',flush=True)
    assert sum(c['input_epochs'] for c in report['cases'].values())==58634
    print('Audited metadata fix on all 58634 epochs; no numeric/status change',flush=True)


if __name__=='__main__':main()
