"""Freeze raw FLOAT gauge equivalence, default controls and real product admission."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import runpy
import subprocess
import tempfile


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def rows(data): return list(csv.DictReader(io.StringIO(data.decode())))
def compare(calibrated, raw, expected):
    assert len(calibrated) == len(raw) == 64
    position = ambiguity = 0.0
    for a, b in zip(calibrated, raw):
        assert (a['week'], a['tow'], a['status']) == (b['week'], b['tow'], 'PPP_FLOAT')
        assert a['phase_bias_mode'] == 'RequireSignalCorrections'
        assert b['phase_bias_mode'] == 'UncalibratedFloat' and int(b['applied_phase_OSBs']) == 0
        assert int(a['applied_phase_OSBs']) == 18
        position = max(position, sum((float(a[k])-float(b[k]))**2 for k in ['x_m','y_m','z_m'])**.5)
        aa, bb = [dict(x.split(':') for x in r['ambiguities_m'].split(';')) for r in [a,b]]
        assert set(aa) == set(bb) == set(expected)
        ambiguity = max(ambiguity, max(abs(float(bb[s])-float(aa[s])-expected[s]) for s in aa))
    assert position < 1e-6 and ambiguity < 1e-6
    return dict(epochs=64, maximum_position_difference_m=position, maximum_ambiguity_shift_error_m=ambiguity, bound_m=1e-6)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--code-root',type=Path,required=True)
    p.add_argument('--antex-root',type=Path,required=True)
    p.add_argument('--cas-root',type=Path,required=True)
    a=p.parse_args(); repo=Path(__file__).resolve().parents[1]; out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    assert not (out/'predeclared-inputs.json').exists()
    prior=repo/'docs/comparison/2026-10-21-kinematic-ppp'
    old=json.loads((prior/'predeclared-inputs.json').read_text())
    helper=runpy.run_path(str(repo/'scripts/run_native_ppc_comparison.py'));save,retain=[helper[k] for k in ['save','retain']]
    sources=list((repo/'src').rglob('*.rs'))+[repo/x for x in [
        'Cargo.toml','Cargo.lock','rust-toolchain.toml','scripts/run_ppp_raw_phase.py',
        'scripts/acquire_ppc_code_rapid.py','scripts/acquire_igs20_antex.py','scripts/select_gps_antex.py',
        'scripts/upstream_ppp_raw_phase.cpp','examples/ppp_kinematic.rs','examples/ppp.rs',
        'examples/ppc_input_audit.rs','examples/ppc_antex.rs','tests/ppp_kinematic.rs','tests/precise.rs','tests/bias.rs','docs/raw-phase-ppp.md']]+list((repo/'tests/fixtures').glob('*ppp*'))
    binaries={n:repo/'target/release/examples'/n for n in ['ppp_kinematic','ppp','ppc_input_audit','ppc_antex']}
    acquisition={}
    for name,root in [('code',a.code_root),('antex',a.antex_root)]:
        acquisition[name]=dict(root=str(root.resolve()),acquisition_sha256=sha(root/'acquisition.json'))
        (out/(name+'-acquisition.json')).write_bytes((root/'acquisition.json').read_bytes())
    for f in ['acquisition-initial-byte-comment-failure.json','acquisition-attempt-0.json','acquire-source-initial.py']:
        if (a.code_root/f).exists():(out/('code-'+f)).write_bytes((a.code_root/f).read_bytes())
    atx=a.antex_root/'gps-subset/gps-satellite-igs20.atx'
    (out/'gps-satellite-igs20.atx.gz').write_bytes(gzip.compress(atx.read_bytes(),mtime=0))
    (out/'antex-selection.json').write_bytes((a.antex_root/'gps-subset/selection.json').read_bytes())
    manifest=dict(base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),source_sha256={str(x.relative_to(repo)):sha(x) for x in sources},binary_sha256={n:sha(b) for n,b in binaries.items()},prior_archive=str(prior.relative_to(repo)),prior_evidence_sha256=sha(prior/'evidence-sha256.json'),acquisition=acquisition,position_and_ambiguity_bound_m=1e-6,default_mode_changed=False,phase_OSB_zero_fallback=False,thresholds_changed=False,whole_native_PPP_executed=False,real_PPC_PPP_executed=False,practical_parity_accepted=False)
    save(out/'predeclared-inputs.json',manifest)
    def execute(name,binary,args):
        assert all(sha(repo/f)==h for f,h in manifest['source_sha256'].items())
        command=[str(binaries[binary])]+list(map(str,args));r=subprocess.run(command,capture_output=True,cwd=repo)
        retain(out,name,command,r.stdout,r.stderr,r.returncode);assert r.returncode==0,r.stderr.decode();return r.stdout
    report=dict(gauge_profiles={},default_controls={},real_PPC_PPP_executed=False,practical_parity_accepted=False)
    expected={r['satellite']:float(r['unmodelled_phase_IF_bias_m']) for r in rows((repo/'tests/fixtures/upstream_ppp_raw_phase.csv').read_bytes())}
    for variant in ['clean','noisy']:
        calibrated=rows(execute('calibrated-'+variant,'ppp_kinematic',['kinematic',variant,'-phase-diagnostics']))
        for mode in ['raw-phase','raw-dsb']:
            raw=rows(execute(mode+'-'+variant,'ppp_kinematic',['kinematic',variant,'-'+mode]))
            report['gauge_profiles'][mode+'-'+variant]=compare(calibrated,raw,expected)
    for name,c in old['moving_profiles'].items():
        data=execute('control-'+name,'ppp_kinematic',c['arguments']);before=gzip.decompress((prior/(name+'.stdout.csv.gz')).read_bytes());assert data==before,name
        report['default_controls'][name]=dict(all_CSV_bytes_identical=True,csv_sha256=hashlib.sha256(data).hexdigest(),prior_capture=name+'.stdout.csv.gz')
    for name,args in old['static_profiles'].items():
        data=execute('control-static-'+name,'ppp',args);capture='static-new-'+name+'.stdout.csv.gz';assert data==gzip.decompress((prior/capture).read_bytes()),name
        report['default_controls']['static-'+name]=dict(all_CSV_bytes_identical=True,csv_sha256=hashlib.sha256(data).hexdigest(),prior_capture=capture)
    admission=rows(execute('code-products','ppc_input_audit',['code-products',a.code_root]));assert len(admission)==9
    assert all(r['GPS_IF_C1W_C2W_declared']==r['bias_APC_IGS20']=='true' and int(r['GPS_sample_states_with_clock'])==32 and int(r['GPS_sample_states_unavailable'])==int(r['gauge_pairs_missing'])==0 for r in admission)
    report['code_product_admission']=admission
    conflicts=[]
    for stamp,week,tow in [('20242020000',2323,518400),('20242050000',2324,172800),('20242160000',2325,518400)]:
        selected=rows(execute('antex-'+stamp,'ppc_antex',[atx,week,tow]));assert len(selected)==32 and all(r['status']=='AVAILABLE' for r in selected)
        svn={r['satellite']:r['svn'] for r in selected}
        for family,root,name in [('CODE',a.code_root,f'COD0OPSRAP_{stamp}_01D_01D_OSB.BIA'),('CAS',a.cas_root,f'CAS0OPSRAP_{stamp}_01D_01D_DCB.BIA')]:
            path=root/name; entries={}
            for line in path.read_bytes().splitlines():
                f=line.split()
                if len(f)>3 and f[0] in [b'OSB',b'DSB'] and len(f[1])==4 and f[1].startswith(b'G') and len(f[2])==3 and f[2].startswith(b'G'):
                    entries.setdefault(f[2].decode(),set()).add(f[1].decode())
            for sat in svn:
                if entries.get(sat)!={svn[sat]}:
                    item=dict(date=stamp,family=family,satellite=sat,ANTEX_SVN=svn[sat],bias_SVNs=sorted(entries.get(sat,set())),bias_sha256=sha(path));conflicts.append(item)
                    assert family=='CAS',item
    report['retained_SVN_conflicts']=conflicts
    save(out/'report.json',report)
    with tempfile.TemporaryDirectory(prefix='ppp-raw-native-',dir='/workspace/.cloud-setup') as d:
        root=Path(d);upstream=Path('/workspace/.cloud-setup/reference/gnssplusplus-library');binary=root/'generator'
        command=['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror','-I'+str(upstream/'include'),str(repo/'scripts/upstream_ppp_raw_phase.cpp'),'-o',str(binary)]
        r=subprocess.run(command,capture_output=True);(out/'native-build.log.gz').write_bytes(gzip.compress(r.stdout+r.stderr,mtime=0));assert r.returncode==0
        native=dict(command=command,returncode=r.returncode,compiler=subprocess.check_output(['g++','--version'],text=True).splitlines()[0],native_constants_header_sha256=sha(upstream/'include/libgnss++/core/constants.hpp'),generator_sha256=sha(repo/'scripts/upstream_ppp_raw_phase.cpp'),binary_sha256=sha(binary),whole_native_PPP_executed=False,repeat_runs=[])
        for i in range(2):
            r=subprocess.run([str(binary),str(root)],capture_output=True);assert r.returncode==0
            files={n:sha(root/n) for n in ['synthetic_ppp_raw_phase_dsb.bsx','upstream_ppp_raw_phase.csv']};assert all(h==sha(repo/'tests/fixtures'/n) for n,h in files.items());native['repeat_runs'].append(dict(returncode=r.returncode,fixture_sha256=files))
        save(out/'native-generator.json',native)
    assert all(sha(b)==manifest['binary_sha256'][n] for n,b in binaries.items())
    print(json.dumps(dict(gauge_profiles=report['gauge_profiles'],byte_identical_controls=len(report['default_controls']),CODE_days=len(admission),SVN_conflicts=conflicts),indent=2))


if __name__=='__main__':main()
