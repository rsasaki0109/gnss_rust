"""Reproduce actual pinned native precise IF ambiguity initialization traces twice."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess


def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args();repo=Path(__file__).resolve().parents[1];out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    prior=repo/'docs/comparison/2026-10-21-ppc-ppp-first/native-initial.json';old=json.loads(prior.read_text());libs=old['native_library_sha256'];assert all(sha(f)==h for f,h in libs.items())
    native=Path('/workspace/.cloud-setup/reference/gnssplusplus-library');source=repo/'scripts/upstream_ppp_if_ambiguity_prior.cpp';binary=out/'native'
    command=['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror','-I'+str(native/'include'),'-I'+str(native/'src/algorithms'),'-I/workspace/.cloud-setup/eigen-3.4.0',str(source),'-Wl,--start-group']+list(libs)+['-Wl,--end-group','-o',str(binary)]
    r=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);(out/'build.log.gz').write_bytes(gzip.compress(r.stdout,mtime=0))
    sources=[source,repo/'scripts/build_ppp_if_ambiguity_fixtures.py'];header=native/'src/algorithms/ppp_internal.hpp'
    manifest=dict(command=command,returncode=r.returncode,source_sha256={str(f.relative_to(repo)):sha(f) for f in sources},native_library_sha256=libs,native_private_header_sha256={str(header):sha(header)},pinned_upstream_commit=old['upstream_commit'],source_archive_sha256=old['source_archive_sha256'],prior_proof=str(prior.relative_to(repo)),prior_proof_sha256=sha(prior),upstream_source_patched=False,whole_native_PPP_executed=False,actual_initializeAmbiguityState_calls_per_repeat=32,repeat_runs=[])
    def save():(out/'native-generator.json').write_text(json.dumps(manifest,indent=2)+'\n')
    save();print('Build',r.returncode,flush=True)
    if r.returncode:print(r.stdout.decode());raise SystemExit(r.returncode)
    manifest['binary_sha256']=sha(binary);manifest['compiler']=subprocess.check_output(['g++','--version'],text=True).splitlines()[0];env={k:v for k,v in os.environ.items() if not k.startswith('GNSS_')};manifest['GNSS_environment_removed']=sorted(set(os.environ)-set(env))
    for repeat in range(2):
        target=out/str(repeat);target.mkdir();cmd=[str(binary),str(target)];r=subprocess.run(cmd,capture_output=True,env=env)
        for stream,data in [('stdout',r.stdout),('stderr',r.stderr)]:(out/f'run{repeat}.{stream}.log.gz').write_bytes(gzip.compress(data,mtime=0))
        assert r.returncode==0,r.stderr.decode();files={f.name:sha(f) for f in sorted(target.iterdir())};manifest['repeat_runs'].append(dict(command=cmd,returncode=0,fixture_sha256=files));save();print('Generation',repeat,'passed',flush=True)
    assert manifest['repeat_runs'][0]['fixture_sha256']==manifest['repeat_runs'][1]['fixture_sha256']
    for name,h in manifest['repeat_runs'][0]['fixture_sha256'].items():
        dest=repo/'tests/fixtures'/name
        if not dest.exists():dest.write_bytes((out/'0'/name).read_bytes())
        assert sha(dest)==h,name


if __name__=='__main__':main()
