"""Select complete GPS satellite ANTEX blocks, retaining every exclusion's provenance."""
import argparse
import hashlib
import json
from pathlib import Path
import re


def sha(data):return hashlib.sha256(data).hexdigest()
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--input",type=Path,required=True);p.add_argument("--output-dir",type=Path,required=True);a=p.parse_args();root=a.output_dir.resolve();root.mkdir(parents=True,exist_ok=False)
    source=a.input.resolve();data=source.read_bytes();lines=data.splitlines(keepends=True);header_end=next(i+1 for i,line in enumerate(lines) if line[60:].strip()==b'END OF HEADER')
    output=list(lines[:header_end]);selected=[];excluded=[];block=[];start=None
    for i,line in enumerate(lines[header_end:],header_end+1):
        label=line[60:].strip()
        if label==b'START OF ANTENNA':
            assert start is None;start=i;block=[]
        if start is not None:block.append(line)
        else:assert not line.strip(),(i,'nonempty record outside antenna')
        if label==b'END OF ANTENNA':
            assert start is not None
            types=[line for line in block if line[60:].strip()==b'TYPE / SERIAL NO'];assert len(types)==1
            kind=types[0][:20].decode('ascii').strip();serial=types[0][20:40].decode('ascii').strip();svn=types[0][40:50].decode('ascii').strip();raw=b''.join(block)
            row=dict(antenna_type=kind,serial=serial,svn=svn,source_start_line=start,source_end_line=i,block_sha256=sha(raw))
            if re.fullmatch(r'G[0-9]{2}',serial):row['excerpt_start_line']=len(output)+1;output.extend(block);selected.append(row)
            else:row['reason']='outside explicit GPS satellite scope; original source retained';excluded.append(row)
            start=None
    assert start is None and selected
    result=b''.join(output);target=root/'gps-satellite-igs20.atx';target.write_bytes(result)
    manifest=dict(source=str(source),source_sha256=sha(data),selector_sha256=sha(Path(__file__).read_bytes()),header_sha256=sha(b''.join(lines[:header_end])),output=target.name,output_sha256=sha(result),selected_blocks=selected,excluded_blocks=excluded,all_selected_block_bytes_unchanged=True,receiver_calibration_selected=False,applied_to_PPP=False,full_ANTEX_format_support_claimed=False)
    (root/'selection.json').write_text(json.dumps(manifest,indent=2)+'\n');print(f'Selected {len(selected)} whole GPS satellite blocks; retained {len(excluded)} exclusions')


if __name__=='__main__':main()
