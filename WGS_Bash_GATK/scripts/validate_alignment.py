#!/usr/bin/env python3
"""Non-destructive prealigned BAM/CRAM header compatibility validation."""
from __future__ import annotations
import subprocess
import sys
from pathlib import Path
from model import reference_dict_path


def validate(path, sample, reference):
    ref=Path(reference)
    index=ref.with_name(ref.name+'.fai')
    reflens={row.split('\t')[0]:int(row.split('\t')[1]) for row in index.read_text().splitlines() if row.strip()}
    dictionary={}
    for line in reference_dict_path(ref).read_text().splitlines():
        if line.startswith('@SQ\t'):
            tags=dict(field.split(':',1) for field in line.split('\t')[1:] if ':' in field)
            dictionary[tags['SN']]=tags
    r=subprocess.run(['samtools','view','-H','-T',str(ref),path],check=True,capture_output=True,text=True)
    sm=set();rg=0; sq={};sort=None
    for line in r.stdout.splitlines():
        fields=line.split('\t')
        tags=dict(field.split(':',1) for field in fields[1:] if ':' in field)
        if line.startswith('@HD\t'): sort=tags.get('SO')
        if line.startswith('@RG\t'):
            rg+=1
            sm.add(tags.get('SM',''))
        if line.startswith('@SQ\t'):
            name=tags.get('SN','')
            if name in sq: raise ValueError('Duplicate @SQ in input')
            sq[name]=tags
    if sort!='coordinate': raise ValueError(f'{path}: expected @HD SO:coordinate, got {sort}')
    if rg==0 or sm!={sample}: raise ValueError(f'{path}: @RG SM must equal {sample}, got {sm}')
    if set(sq)!=set(reflens):
        raise ValueError(f'{path}: @SQ contig set differs from selected reference')
    for name,tags in sq.items():
        if int(tags.get('LN','-1'))!=reflens[name]:
            raise ValueError(f'{path}: reference length mismatch for {name}')
        ref_md5=dictionary[name].get('M5')
        bam_md5=tags.get('M5')
        if ref_md5 and bam_md5 and ref_md5!=bam_md5:
            raise ValueError(f'{path}: reference MD5 differs for {name}')
    return True


if __name__=='__main__':
    try: validate(*sys.argv[1:4])
    except Exception as e: sys.exit(f'ERROR: {e}')
