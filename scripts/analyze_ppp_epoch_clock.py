"""Post-hoc, read-only worst-point diagnostics; never an acceptance/profile selector."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import runpy


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def derive(root,repo):
    m=json.loads((root/'predeclared-inputs.json').read_text());targets=runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))['reference_targets'];records=[]
    for case,c in m['cases'].items():
        with Path(c['reference']).open(newline='') as f:refs=[{k.strip():v for k,v in r.items()} for r in csv.DictReader(f)]
        references={(int(r['GPS Week']),Decimal(r['GPS TOW (s)'])):r for r in refs}
        for pair in m['pairs']:
            for head in ['default','native-head']:
                def read(clock):return list(csv.DictReader(io.StringIO(gzip.decompress((root/(case+'-'+pair+'-'+head+'-'+clock+'.stdout.csv.gz')).read_bytes()).decode())))
                rows=read('independent');control=read('random-walk');assert len(rows)==len(control);worst={}
                for i,r in enumerate(rows):
                    if r['status']!='PPP_FLOAT':continue
                    key=(int(r['week']),Decimal(r['tow']).quantize(Decimal('.000001')));reference=references[key];hypotheses,_,_=targets(reference,c['city']);xyz=[float(r[k]) for k in ['x_m','y_m','z_m']]
                    for hypothesis,target in hypotheses.items():
                        delta=[x-y for x,y in zip(xyz,target)];error=math.sqrt(sum(x*x for x in delta))
                        if hypothesis in worst and worst[hypothesis]['error_3D_m']>=error:continue
                        previous=rows[i-1] if i else {};gap=None
                        if previous.get('committed_week'):gap=float((Decimal(r['week'])-Decimal(previous['committed_week']))*604800+Decimal(r['tow'])-Decimal(previous['committed_tow']))
                        worst[hypothesis]=dict(case=case,pair=pair,head_prior=head,hypothesis=hypothesis,row_index=i,solver_time=dict(week=r['week'],tow=r['tow']),reference_time=dict(week=reference['GPS Week'],tow=reference['GPS TOW (s)']),error_3D_m=error,error_ECEF_m=delta,prior_clock_control_status=control[i]['status'],common_success=control[i]['status']=='PPP_FLOAT',gap_from_previous_commit_s=gap,solver_fields=r,reference_ENU_speed_mps=math.sqrt(sum(float(reference[k])**2 for k in ['East Velocity (m/s)','North Velocity (m/s)','Up Velocity (m/s)'])),estimated_ECEF_speed_mps=math.sqrt(sum(float(r[k])**2 for k in ['vx_mps','vy_mps','vz_mps'])))
                records.extend(worst[k] for k in sorted(worst))
    return dict(post_hoc_diagnostic_only=True,used_to_select_solver_profile=False,thresholds_changed=False,source_sha256=sha(Path(__file__)),reference_targets_source_sha256=sha(repo/'scripts/run_ppc_rtk_float.py'),predeclared_inputs_sha256=sha(root/'predeclared-inputs.json'),records=records)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--archive',type=Path,required=True);p.add_argument('--verify',action='store_true');a=p.parse_args();repo=Path(__file__).resolve().parents[1];root=a.archive.resolve();result=derive(root,repo);path=root/'worst-points.json'
    if a.verify:assert json.loads(path.read_text())==result
    else:assert not path.exists();path.write_text(json.dumps(result,indent=2)+'\n')
    print(('Verified' if a.verify else 'Retained')+' '+str(len(result['records']))+' worst-origin diagnostic records; no solver selection')


if __name__=='__main__':main()
