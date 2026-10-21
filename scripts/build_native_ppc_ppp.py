"""Build a whole native PPP adapter against the pinned, unchanged reference libraries."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess


def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args();repo=Path(__file__).resolve().parents[1];out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True);assert not (out/'build.json').exists()
    proof=repo/'docs/comparison/2026-10-21-ppc-ppp-first/native-initial.json';prior=json.loads(proof.read_text());libs=prior['native_library_sha256'];assert all(sha(p)==h for p,h in libs.items())
    source=repo/'scripts/native_ppc_ppp.cpp';binary=out/'native_ppc_ppp';cmd=['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror','-I/workspace/.cloud-setup/reference/gnssplusplus-library/include','-I/workspace/.cloud-setup/eigen-3.4.0',str(source),'-Wl,--start-group']+list(libs)+['-Wl,--end-group','-o',str(binary)]
    r=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);log=out/'build.log.gz';log.write_bytes(gzip.compress(r.stdout,mtime=0));m=dict(command=cmd,returncode=r.returncode,log=log.name,log_sha256=sha(log),source=str(source.relative_to(repo)),source_sha256=sha(source),builder_source='scripts/build_native_ppc_ppp.py',builder_sha256=sha(Path(__file__)),compiler=subprocess.check_output(['g++','--version'],text=True).splitlines()[0],pinned_upstream_commit=prior['upstream_commit'],source_archive_sha256=prior['source_archive_sha256'],native_library_sha256=libs,earlier_native_build_proof=str(proof.relative_to(repo)),earlier_native_build_proof_sha256=sha(proof),upstream_source_patched=False,native_observation_RINEX_pipeline=False,shared_predecoded_input=True)
    if r.returncode:print(r.stdout.decode())
    else:m['binary_sha256']=sha(binary)
    (out/'build.json').write_text(json.dumps(m,indent=2)+'\n');print('Native PPP adapter build:',r.returncode,flush=True);raise SystemExit(r.returncode)


if __name__=='__main__':main()
