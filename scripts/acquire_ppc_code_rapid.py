"""Acquire matching CODE rapid orbit/clock/OSB products on all nine PPC days."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import date,datetime,timedelta
import gzip
import hashlib
import json
from pathlib import Path
import subprocess


def sha(data):return hashlib.sha256(data).hexdigest()
def bias_lines(data):
    active=None;non_ascii=[];lines=[]
    for i,line in enumerate(data.splitlines(),1):
        trimmed=line.strip()
        for label in [b'FILE/REFERENCE',b'FILE/COMMENT',b'INPUT/ACKNOWLEDGMENTS']:
            if trimmed==b'+'+label:active=b'-'+label
        if not line.isascii():
            if active is None and not trimmed.startswith(b'*'):raise ValueError(f'non-ASCII technical BIA record at line {i}')
            non_ascii.append(i)
        lines.append(line.decode('ascii',errors='backslashreplace'))
        if active==trimmed:active=None
    return lines,non_ascii


def acquire(root,day,kind,reuse=False):
    stamp=day.strftime('%Y%j');week=(day-date(1980,1,6)).days//7
    suffix={'orbit':'01D_05M_ORB.SP3.gz','clock':'01D_30S_CLK.CLK.gz','osb':'01D_01D_OSB.BIA.gz'}[kind]
    name=f'COD0OPSRAP_{stamp}0000_{suffix}';url=f'https://igs.bkg.bund.de/root_ftp/IGSac/products/{week}/{name}'
    target=root/name;headers=root/(name+'.headers')
    command=['curl','--location','--fail','--silent','--show-error','--retry','2','--retry-delay','1','--connect-timeout','10','--max-time','60','--dump-header',str(headers),'--output',str(target),url]
    reused=reuse and target.exists() and target.read_bytes()[:2]==b'\x1f\x8b'
    if reused:code,stderr=0,''
    else:
        result=subprocess.run(command,capture_output=True);code,stderr=result.returncode,result.stderr.decode()
    item=dict(date=str(day),kind=kind,filename=name,url=url,command=None if reused else command,transport_reused=reused,returncode=code,stderr=stderr,acquired=False)
    if code:return item
    try:
        data=target.read_bytes()
        with gzip.open(target,'rb') as stream:expanded=stream.read(256*1024*1024+1)
        assert len(expanded)<=256*1024*1024
        if kind=='osb':
            lines,legacy_lines=bias_lines(expanded);item['non_ASCII_free_text_lines']=legacy_lines;item['free_text_encoding_guessed']=False
        else:lines=expanded.decode('ascii').splitlines()
        if kind=='orbit':
            assert lines[0].startswith('#') and any(s.strip()=='EOF' for s in lines)
            first=next(s for s in lines if s.startswith('*'));assert tuple(map(int,first.split()[1:4]))==(day.year,day.month,day.day)
            records=sum(s.startswith('PG') for s in lines)
        elif kind=='clock':
            assert 'RINEX VERSION / TYPE' in lines[0] and any('END OF HEADER' in s for s in lines)
            first=next(s for s in lines if s.startswith('AS G'));assert tuple(map(int,first.split()[2:5]))==(day.year,day.month,day.day)
            records=sum(s.startswith('AS G') for s in lines)
        else:
            assert lines[0].startswith('%=BIA') and any(s.startswith('-BIAS/SOLUTION') for s in lines)
            assert any(s.strip().split()==['BIAS_MODE','ABSOLUTE'] for s in lines)
            records=sum(s.startswith(' OSB') for s in lines)
            item['phase_OSB_rows']=sum(s.startswith(' OSB') and any(t.startswith('L') and len(t)==3 for t in s.split()[3:6]) for s in lines)
        assert records>0
        path=target.with_suffix('');path.write_bytes(expanded)
        item.update(acquired=True,compressed_bytes=len(data),compressed_sha256=sha(data),expanded_filename=path.name,expanded_bytes=len(expanded),expanded_sha256=sha(expanded),headers_sha256=sha(headers.read_bytes()),relevant_records=records,first_line=lines[0])
    except (AssertionError,ValueError,OSError,EOFError,StopIteration) as e:item['validation_error']=type(e).__name__+': '+str(e)
    return item


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--resume',action='store_true');a=p.parse_args();root=a.output_dir.resolve();root.mkdir(parents=True,exist_ok=a.resume)
    previous=None
    if a.resume:
        data=(root/'acquisition.json').read_bytes();number=len(list(root.glob('acquisition-attempt-*.json')));saved=root/f'acquisition-attempt-{number}.json';saved.write_bytes(data);previous=dict(filename=saved.name,sha256=sha(data))
    centres=[date(2024,7,20),date(2024,7,23),date(2024,8,3)];days=sorted({d+timedelta(days=n) for d in centres for n in [-1,0,1]})
    m=dict(source_sha256=sha(Path(__file__).read_bytes()),source='BKG public IGS analysis-center product mirror',product_family='matching COD0OPSRAP orbit/clock/absolute OSB',central_dates=list(map(str,centres)),padding_days=1,previous_attempt=previous,products=[],applied_to_PPP=False,clock_signal_reference_verified=False,phase_OSB_signal_coverage_verified=False,performance_accepted=False,started_utc=datetime.now().isoformat(timespec='seconds')+'Z')
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures=[pool.submit(acquire,root,day,kind,a.resume) for day in days for kind in ['orbit','clock','osb']]
        for future in as_completed(futures):
            item=future.result();m['products'].append(item);(root/'acquisition.json').write_text(json.dumps(m,indent=2)+'\n');print({k:item[k] for k in ['filename','returncode','acquired']},flush=True)
    m['products'].sort(key=lambda x:(x['date'],x['kind']));m['acquired_count']=sum(i['acquired'] for i in m['products']);m['completed_utc']=datetime.now().isoformat(timespec='seconds')+'Z';(root/'acquisition.json').write_text(json.dumps(m,indent=2)+'\n');print(f"Acquired {m['acquired_count']}/27 matching CODE rapid products");raise SystemExit(0 if m['acquired_count']==27 else 1)


if __name__=='__main__':main()
