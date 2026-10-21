"""Build the committed static PPP in a temporary copy, without a Git checkout/worktree."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess


def sha(data):return hashlib.sha256(data).hexdigest()
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output-dir",type=Path,required=True);a=p.parse_args()
    repo=Path(__file__).resolve().parents[1];root=a.output_dir.resolve();root.mkdir(parents=True,exist_ok=False)
    commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=repo,text=True).strip()
    paths=subprocess.check_output(["git","ls-tree","-r","--name-only",commit,"src"],cwd=repo,text=True).splitlines()
    paths=[f for f in paths if f.endswith(".rs")]+["examples/ppp.rs"]
    example=subprocess.check_output(["git","show",commit+":examples/ppp.rs"],cwd=repo).decode()
    paths+=sorted({str((Path("examples")/f).as_posix()).replace("examples/../","") for f in re.findall(r'"(\.\./tests/fixtures/[^\"]+)"',example)})
    sources={}
    for name in paths:
        data=subprocess.check_output(["git","show",commit+":"+name],cwd=repo);path=root/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data);sources[name]=sha(data)
    commands=[['rustc','--edition=2024','--crate-name','gnss_rust','--crate-type','rlib','-C','opt-level=3','src/lib.rs','-o','libgnss_rust.rlib'],
        ['rustc','--edition=2024','--crate-name','legacy_ppp','-C','opt-level=3','examples/ppp.rs','--extern','gnss_rust=libgnss_rust.rlib','-o','ppp']]
    checks=[]
    for i,command in enumerate(commands):
        r=subprocess.run(command,cwd=root,capture_output=True);data=r.stdout+r.stderr;name=f"build-{i}.log.gz";(root/name).write_bytes(gzip.compress(data,mtime=0));checks.append(dict(command=command,returncode=r.returncode,log=name,uncompressed_sha256=sha(data)));assert r.returncode==0,data.decode()
    result=dict(base_commit=commit,source_sha256=sources,builder_sha256=sha(Path(__file__).read_bytes()),checks=checks,compiler=subprocess.check_output(['rustc','--version','--verbose'],text=True),binary_sha256=sha((root/'ppp').read_bytes()),library_sha256=sha((root/'libgnss_rust.rlib').read_bytes()),scope="committed static solver control only; no runtime acceptance")
    (root/'build.json').write_text(json.dumps(result,indent=2)+'\n');print(f"Built static control from {commit}; {len(sources)} retained sources/fixtures")


if __name__=="__main__":main()
