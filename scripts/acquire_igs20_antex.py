"""Acquire public IGS20 ANTEX with source/transport hashes; do not apply calibrations."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output-dir",type=Path,required=True);a=p.parse_args();root=a.output_dir.resolve();root.mkdir(parents=True,exist_ok=False)
    sources=[('igs','https://files.igs.org/pub/station/general/igs20.atx'),('bkg','https://igs.bkg.bund.de/root_ftp/IGS/station/general/igs20.atx')]
    m=dict(source_sha256=sha(Path(__file__)),product='IGS20 ANTEX',downloaded=[],applied_to_PPP=False,receiver_calibration_identified=False,performance_accepted=False)
    for tag,url in sources:
        output=root/(tag+'-igs20.atx');headers=root/(tag+'-igs20.atx.headers')
        command=['curl','--location','--fail','--silent','--show-error','--connect-timeout','10','--max-time','60','--max-filesize','67108864','--dump-header',str(headers),'--output',str(output),url]
        r=subprocess.run(command,capture_output=True);item=dict(source=tag,url=url,command=command,returncode=r.returncode,stderr=r.stderr.decode(),acquired=False)
        if r.returncode==0:
            data=output.read_bytes();lines=data.decode('ascii').splitlines();assert 'ANTEX VERSION / SYST' in lines[0]
            item.update(acquired=True,filename=output.name,bytes=len(data),sha256=sha(output),headers_sha256=sha(headers),first_line=lines[0],antenna_blocks=sum('START OF ANTENNA' in s for s in lines))
        m['downloaded'].append(item);(root/'acquisition.json').write_text(json.dumps(m,indent=2)+'\n');print({k:item[k] for k in ['source','returncode','acquired']},flush=True)
    successes=[p for p in m['downloaded'] if p['acquired']];assert successes
    m['all_acquired_sources_identical']=len(successes)==2 and successes[0]['sha256']==successes[1]['sha256'];(root/'acquisition.json').write_text(json.dumps(m,indent=2)+'\n');print('Acquired IGS20 ANTEX; calibration admission remains pending')


if __name__=='__main__':main()
