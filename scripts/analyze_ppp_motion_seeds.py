"""Read-only worst-origin and frozen-native comparisons for all fixed seed modes."""
import argparse
import collections
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import runpy


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rows(p):return list(csv.DictReader(io.StringIO(gzip.decompress(p.read_bytes()).decode())))
def service(rr):return [dict(r,status=r['service_status'],**{k:r['service_'+k] for k in ['x_m','y_m','z_m','vx_mps','vy_mps','vz_mps']}) for r in rr]
def normalized(rr):return [dict(r,week=str(int(r['week'])),tow=str(Decimal(r['tow']).quantize(Decimal('.000001')))) for r in rr]


def common_valid(current,control,references,city,scorer):
    aa,bb=[{(r['week'],r['tow']):r for r in normalized(rr) if r['status'] in ['SPP','PPP_FLOAT']} for rr in [control,current]];common=aa.keys()&bb.keys()
    def score(rr):return scorer([dict(r,status='FLOAT') for r in rr],references,city)
    return dict(common_valid_epochs=len(common),current_only_valid=len(bb.keys()-aa.keys()),control_only_valid=len(aa.keys()-bb.keys()),common_current_status_counts=dict(collections.Counter(bb[k]['status'] for k in common)),common_control_status_counts=dict(collections.Counter(aa[k]['status'] for k in common)),common_current_reference_errors=score([bb[k] for k in sorted(common)]),common_control_reference_errors=score([aa[k] for k in sorted(common)]),current_only_reference_errors=score([bb[k] for k in sorted(bb.keys()-aa.keys())]),control_only_reference_errors=score([aa[k] for k in sorted(aa.keys()-bb.keys())]),SPP_relabelled_PPP=False)


def derive(root,repo):
    m=json.loads((root/'predeclared-inputs.json').read_text());targets=runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))['reference_targets'];scorer=runpy.run_path(str(repo/'scripts/run_ppc_rtk_float.py'))['score'];initial=runpy.run_path(str(repo/'scripts/run_ppc_ppp_initial_prior.py'));native=repo/'docs/comparison/2026-10-21-native-ppp-common-inputs';records=[];comparisons={};known_point={};service_comparisons={}
    for case,c in m['cases'].items():
        with Path(c['reference']).open(newline='') as f:refs=[{k.strip():v for k,v in r.items()} for r in csv.DictReader(f)]
        references={(int(r['GPS Week']),Decimal(r['GPS TOW (s)'])):r for r in refs};assert len(references)==len(refs)
        for pair in m['pairs']:
            outputs={profile:rows(root/(case+'-'+pair+'-'+profile+'.stdout.csv.gz')) for profile in m['profiles']}
            for profile,rr in outputs.items():
                if profile!='shadow':service_comparisons[case+'-'+pair+'-'+profile]=common_valid(service(rr),service(outputs['shadow']),references,c['city'],scorer)
                if pair=='1C2W':
                    for native_seed in ['native_builtin','native_header_zero']:
                        nr=rows(native/(case+'-'+native_seed+'.stdout.csv.gz'));key=case+'-'+profile+'__'+native_seed
                        comparisons[key]=dict(PPP_only=initial['compare'](normalized(rr),normalized(nr),references,c['city'],scorer),combined_valid=common_valid(service(rr),nr,references,c['city'],scorer))
                worst={}
                for i,r in enumerate(rr):
                    key=(int(r['week']),Decimal(r['tow']).quantize(Decimal('.000001')));ref=references[key];hypotheses,_,_=targets(ref,c['city'])
                    if case=='nagoya-run3' and pair=='1C2L' and key==(2325,Decimal('554542.200000')):known_point[profile]=dict(row_index=i,fields=r)
                    for output,status,prefix in [('PPP',r['status'],''),('combined',r['service_status'],'service_'),('SPP_all',r['seed_status'],'seed_')]:
                        if status not in ['PPP_FLOAT','SPP']:continue
                        xyz=[float(r[prefix+k+'_m']) for k in ['x','y','z']]
                        for hypothesis,target in hypotheses.items():
                            error=math.sqrt(sum((x-y)**2 for x,y in zip(xyz,target)));label=(output,hypothesis)
                            if label in worst and worst[label]['error_3D_m']>=error:continue
                            previous=rr[i-1] if i else {};gap=None
                            if previous.get('committed_week'):gap=float((Decimal(r['week'])-Decimal(previous['committed_week']))*604800+Decimal(r['tow'])-Decimal(previous['committed_tow']))
                            velocity=[r[prefix+k+'_mps'] for k in ['vx','vy','vz']]
                            worst[label]=dict(case=case,pair=pair,profile=profile,output=output,status=status,hypothesis=hypothesis,row_index=i,error_3D_m=error,reference_time=dict(week=ref['GPS Week'],tow=ref['GPS TOW (s)']),gap_from_previous_commit_s=gap,reference_speed_mps=math.sqrt(sum(float(ref[k])**2 for k in ['East Velocity (m/s)','North Velocity (m/s)','Up Velocity (m/s)'])),estimated_speed_mps=math.sqrt(sum(float(v)**2 for v in velocity)) if all(velocity) else None,fields=r,code_rows_minus_XYZ_clock_ZTD=int(r['satellites'])-5 if r['status']=='PPP_FLOAT' else None,PPP_to_SPP_distance_m=math.sqrt(sum((float(r[k+'_m'])-float(r['seed_'+k+'_m']))**2 for k in ['x','y','z'])) if r['status']=='PPP_FLOAT' and r['seed_status']=='SPP' else None)
                records.extend(worst[k] for k in sorted(worst))
    return dict(post_hoc_diagnostic_only=True,used_to_select_profile=False,source_sha256=sha(__file__),reference_targets_source_sha256=sha(repo/'scripts/run_ppc_rtk_float.py'),predeclared_inputs_sha256=sha(root/'predeclared-inputs.json'),worst_origin_records=records,previous_134m_epoch=known_point,service_vs_shadow=service_comparisons,frozen_native_archive=str(native.relative_to(repo)),frozen_native_evidence_sha256=sha(native/'evidence-sha256.json'),native_GPS_1C2W_comparisons=comparisons,native_feature_match=False,native_feature_difference='Native prior GPS common-input/raw-phase baseline versus Rust calibrated-phase GPS PPP and optional GPS/Galileo native-variance-QC SPP/Doppler seeds. Native SPP broadcast/no DCB; native precise effective clock random walk versus Rust IndependentEpoch; native built-in/header-zero initialization versus Rust native-head/PVA. No new native whole-PPP calls.',runtime_memory_acceptance=False,SPP_relabelled_PPP=False,practical_parity_accepted=False)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--archive',type=Path,required=True);p.add_argument('--verify',action='store_true');a=p.parse_args();repo=Path(__file__).resolve().parents[1];root=a.archive.resolve();result=derive(root,repo);path=root/'post-hoc-diagnostics.json'
    if a.verify:assert json.loads(path.read_text())==result
    else:assert not path.exists();path.write_text(json.dumps(result,indent=2)+'\n')
    print(('Verified' if a.verify else 'Retained')+' '+str(len(result['worst_origin_records']))+' worst-origin records and '+str(len(result['native_GPS_1C2W_comparisons']))+' frozen-native comparisons; no acceptance')


if __name__=='__main__':main()
