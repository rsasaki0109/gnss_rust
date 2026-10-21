"""Rebuild retained historical static sources without Git history or a checkout."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess


def sha(data):return hashlib.sha256(data).hexdigest()
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output-dir",type=Path,required=True);a=p.parse_args()
    archive=Path(__file__).resolve().parent;repo=archive.parents[2];root=a.output_dir.resolve();root.mkdir(parents=True,exist_ok=False)
    old=json.loads((archive/"old-static-build.json").read_text())
    for name,h in old["source_sha256"].items():
        data=gzip.decompress((archive/"old-static-ppp.rs.gz").read_bytes()) if name=="src/ppp.rs" else (repo/name).read_bytes()
        assert sha(data)==h,name
        target=root/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
    checks=[]
    for i,check in enumerate(old["checks"]):
        r=subprocess.run(check["command"],cwd=root,capture_output=True);data=r.stdout+r.stderr;name=f"build-{i}.log.gz";(root/name).write_bytes(gzip.compress(data,mtime=0));checks.append(dict(command=check["command"],returncode=r.returncode,log=name,uncompressed_sha256=sha(data)));assert r.returncode==0,data.decode()
    result=dict(old,checks=checks,builder_sha256=sha(Path(__file__).read_bytes()),builder_source=str(Path(__file__).resolve().relative_to(repo)),original_build_sha256=sha((archive/"old-static-build.json").read_bytes()),compiler=subprocess.check_output(["rustc","--version","--verbose"],text=True),binary_sha256=sha((root/"ppp").read_bytes()),library_sha256=sha((root/"libgnss_rust.rlib").read_bytes()))
    (root/"build.json").write_text(json.dumps(result,indent=2)+"\n");print("Rebuilt retained static reference",old["base_commit"])


if __name__=="__main__":main()
