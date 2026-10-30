#!/usr/bin/env python3
"""Prepare and download COMPLETE Illumina paired-end WGS for beluga/polar bear.

Two-step workflow: discover records a frozen selection and ENA MD5/bytes; download
retrieves ENTIRE FASTQ files, verifies them, then generates production config.
No FASTQ subsampling, reference modification, or interval restriction.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

import yaml

RUN_RE = re.compile(r'^[SED]RR\d+$')
ENA = 'https://www.ebi.ac.uk/ena/portal/api/filereport'
FIELDS = ('run_accession','sample_accession','scientific_name','library_strategy',
          'library_layout','instrument_platform','instrument_model','library_name',
          'experiment_accession','fastq_ftp','fastq_md5','fastq_bytes',
          'base_count','study_accession')
REFERENCE = {
 'Delphinapterus_leucas': ('GCA_059928025.1', 'Delphinapterus_leucas_GCA_059928025.1',21),
 'Ursus_maritimus': ('GCA_062004589.1','Ursus_maritimus_GCA_062004589.1',36),
}
DEFAULT_BEAR = ('SRR15170741','SRR15170742')


def fetch_rows(accession: str):
    q=urllib.parse.urlencode({'accession':accession,'result':'read_run','fields':','.join(FIELDS),'format':'tsv'})
    req=urllib.request.Request(ENA+'?'+q,headers={'User-Agent':'FullWGS-Downloader/1.0'})
    with urllib.request.urlopen(req,timeout=120) as r:
        body=r.read().decode('utf-8-sig')
    parsed=list(csv.DictReader(io.StringIO(body),delimiter='\t'))
    if not parsed or not parsed[0].get('run_accession'):
        raise ValueError(f'ENA returned no runs for {accession}: {body[:250]}')
    return parsed


def fastq_url(url):
    p=urllib.parse.urlparse(url if '://' in url else 'https://'+url)
    if p.hostname != 'ftp.sra.ebi.ac.uk' or not p.path.startswith('/vol1/') or '..' in p.path:
        raise ValueError(f'Unexpected ENA URL: {url}')
    return 'https://ftp.sra.ebi.ac.uk'+p.path


def normalize(row, species):
    if (not RUN_RE.fullmatch(row.get('run_accession','')) or
        row.get('library_strategy','').upper()!='WGS' or
        row.get('library_layout','').upper()!='PAIRED' or
        row.get('instrument_platform','').upper()!='ILLUMINA' or
        row.get('scientific_name','').strip().lower().replace(' ','_')!=species.lower() or
        not row.get('sample_accession')):
        return None
    paths=row.get('fastq_ftp','').split(';')
    checks=row.get('fastq_md5','').split(';')
    sizes=row.get('fastq_bytes','').split(';')
    if len(paths)!=2 or len(checks)!=2 or len(sizes)!=2:
        return None  # Reject paired runs with ambiguous or additional orphan FASTQs
    if any(not re.fullmatch('[0-9a-fA-F]{32}',x) for x in checks):
        return None
    if any(not x.isdecimal() or int(x)<=0 for x in sizes):
        return None
    try:
        bases=int(row.get('base_count') or 0)
    except ValueError:
        bases=0
    if bases<=0: return None
    urls=[fastq_url(x) for x in paths]
    if not urls[0].endswith('_1.fastq.gz') or not urls[1].endswith('_2.fastq.gz'):
        return None
    return {'species':species,'run':row['run_accession'],'biosample':row['sample_accession'],
            'experiment':row.get('experiment_accession',''),
            'library':row.get('library_name','') or row.get('experiment_accession','') or 'LIB1',
            'instrument':row.get('instrument_model',''),
            'bases':bases,'urls':urls,'md5':[x.lower() for x in checks],
            'bytes':[int(x) for x in sizes],
            'study':row.get('study_accession','')}


def choose_beluga(rows, n, min_bases):
    groups=defaultdict(list)
    for r in rows:
        v=normalize(r,'Delphinapterus_leucas')
        if v: groups[v['biosample']].append(v)
    eligible=[]
    for biosample,runs in groups.items():
        libs=set(x['library'] for x in runs)
        if len(libs)!=1: continue  # current pipeline permits one library per sample
        if len({x['run'] for x in runs})!=len(runs): continue
        total=sum(x['bases'] for x in runs)
        if total>=min_bases:
            eligible.append((sum(sum(x['bytes']) for x in runs),biosample,sorted(runs,key=lambda r:r['run'])))
    eligible.sort(key=lambda x:(x[0],x[1]))
    if len(eligible)<n:
        raise ValueError(f'Only {len(eligible)} distinct beluga BioSamples meet the requested filters; need {n}. '
                         'Review ENA metadata and choose another BioProject if necessary.')
    return [x[2] for x in eligible[:n]]


def discover(project_acc, beluga_n=2, bear_runs=DEFAULT_BEAR, min_bases=50_000_000_000):
    all_beluga=fetch_rows(project_acc)
    beluga_groups=choose_beluga(all_beluga,beluga_n,min_bases)
    bear_groups=[]
    biosamples=set()
    for run in bear_runs:
        rows=fetch_rows(run)
        rec=[normalize(x,'Ursus_maritimus') for x in rows if x['run_accession']==run]
        rec=[x for x in rec if x is not None]
        if len(rec)!=1: raise ValueError(f'{run} not verified as unique paired Illumina WGS bear run')
        if rec[0]['biosample'] in biosamples:
            raise ValueError('Bear runs share a BioSample; do not mistake separate runs for separate bears')
        biosamples.add(rec[0]['biosample'])
        bear_groups.append(rec)
    samples=[]
    for i,group in enumerate(beluga_groups,1):
        samples.append({'sample':f'BELUGA_WGS{i:02d}','species':'Delphinapterus_leucas',
                        'biosample':group[0]['biosample'],'runs':group})
    for i,group in enumerate(bear_groups,1):
        samples.append({'sample':f'BEAR_WGS{i:02d}','species':'Ursus_maritimus',
                        'biosample':group[0]['biosample'],'runs':group})
    return {'source':'ENA filereport','beluga_project':project_acc,'samples':samples}


def report(plan):
    total=0
    for s in plan['samples']:
        n=sum(sum(r['bytes']) for r in s['runs']); total+=n
        bases=sum(r['bases'] for r in s['runs'])
        print(f"{s['sample']:15s} {s['species']:24s} biosample={s['biosample']} "
              f"runs={','.join(r['run'] for r in s['runs'])} "
              f"FASTQ={n/1e9:.1f} GB; bases={bases/1e9:.1f} Gb")
    print(f'Total downloadable gzip FASTQ: {total/1e9:.1f} GB ({total:,} bytes)')


def fasta_info(refroot: Path, species: str):
    accession, folder, expected = REFERENCE[species]
    fadir=refroot/folder/'ncbi_dataset'/'data'/accession
    fas=list(fadir.glob('*_genomic.fna'))
    if len(fas)!=1 or not fas[0].is_file():
        raise ValueError(f'Exactly one reference FASTA required in {fadir}; found {len(fas)}')
    fasta=fas[0]
    fai=Path(str(fasta)+'.fai')
    if not fai.is_file():raise ValueError(f'Missing reference index {fai}')
    seq=fadir/'sequence_report.jsonl'
    if not seq.is_file():
        raise ValueError(f'Missing {seq}; retrieve NCBI sequence report for {accession} before production.')
    records=[]
    for l in seq.read_text(encoding='utf8').splitlines():
        if l.strip(): records.append(json.loads(l))
    faiseqs={}
    order=[]
    for l in fai.read_text().splitlines():
        cols=l.split('\t');faiseqs[cols[0]]=int(cols[1]);order.append(cols[0])
    chosen=set()
    for r in records:
        chr_name=str(r.get('chrName','')).strip()
        loc=str(r.get('assignedMoleculeLocationType','')).lower()
        role=str(r.get('role','')).lower()
        if re.fullmatch(r'(?:chr)?\d+',chr_name,re.IGNORECASE) and loc=='chromosome' and role=='assembled-molecule':
            acc=r.get('genbankAccession','') or r.get('refseqAccession','')
            if acc not in faiseqs:
                raise ValueError(f'Sequence report chromosome {chr_name}, accession {acc} not in FASTA')
            if faiseqs[acc]!=int(r.get('length',-1)):
                raise ValueError(f'Chromosome length discrepancy: {acc}')
            chosen.add(acc)
    selected=[c for c in order if c in chosen]
    if len(selected)!=expected:
        raise ValueError(f'{species}: selected {len(selected)} autosomes (expected {expected}); '
                         'DO NOT guess by length. Inspect sequence_report.jsonl and chromosome labels.')
    return str(fasta),selected


def md5_file(path):
    h=hashlib.md5()
    with path.open('rb') as r:
        for chunk in iter(lambda:r.read(16*1024*1024),b''):
            h.update(chunk)
    return h.hexdigest()


def downloader(url,dest,expected_md5,expected_bytes):
    dest.parent.mkdir(parents=True,exist_ok=True)
    if dest.exists():
        if dest.stat().st_size==expected_bytes and md5_file(dest)==expected_md5:
            print('Already verified:',dest,flush=True);return
        raise ValueError(f'Existing destination has wrong size or MD5: {dest}; investigate before deleting')
    partial=Path(str(dest)+'.part')
    if not shutil.which('curl'):
        raise RuntimeError('curl required for resumable full FASTQ transfer')
    print('Downloading:',dest, 'expected bytes',expected_bytes,flush=True)
    cmd=['curl','--fail','--location','--retry','8','--retry-delay','5',
         '--connect-timeout','30','--continue-at','-','--output',str(partial),url]
    subprocess.run(cmd,check=True)
    if partial.stat().st_size != expected_bytes:
        raise RuntimeError(f'Download size mismatch for {partial}: expected {expected_bytes}')
    if md5_file(partial)!=expected_md5:
        raise RuntimeError(f'MD5 mismatch for {partial}; do NOT accept corrupted FASTQ')
    partial.rename(dest)


def write_configuration(plan,project,refroot,reads,output,workdir):
    refs={}
    for species in REFERENCE:
        f,contigs=fasta_info(refroot,species)
        refs[species]={'fasta':f,'contigs':contigs}
    source=project/'config'/'config.yaml'
    if not source.is_file():raise ValueError(f'Not found: {source}')
    cfg=yaml.safe_load(source.read_text(encoding='utf8'))
    cfg['references']=refs
    cfg['workdir']=str(workdir.resolve())
    cfg['samples']=str((output/'full_samples.tsv').resolve())
    cfg.setdefault('preprocessing',{})['trim_fastq']=False
    output.mkdir(parents=True,exist_ok=True)
    (output/'full_wgs.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False,allow_unicode=True))
    cols=['sample','species','input_type','library','lane','r1','r2','bam','cram']
    with (output/'full_samples.tsv').open('w',newline='',encoding='utf8') as handle:
        writer=csv.DictWriter(handle,fieldnames=cols,delimiter='\t');writer.writeheader()
        for s in plan['samples']:
            for i,rec in enumerate(s['runs'],1):
                dest=reads/s['species']/rec['run']
                writer.writerow({'sample':s['sample'],'species':s['species'],'input_type':'fastq',
                                 'library':'LIB1','lane':f'L{i:03d}',
                                 'r1':str(dest/f"{rec['run']}_1.fastq.gz"),
                                 'r2':str(dest/f"{rec['run']}_2.fastq.gz"),
                                 'bam':'','cram':''})
    (output/'autosome_selection.tsv').write_text('species\tchromosome_accession\n'+''.join(
        f'{sp}\t{c}\n' for sp,r in refs.items() for c in r['contigs']))
    return output/'full_wgs.yaml'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('references','discover','download','config','report'))
    parser.add_argument('--plan',type=Path,required=True,help='JSON file freezing exact sample/run selection')
    parser.add_argument('--project',type=Path)
    parser.add_argument('--reference-root',type=Path)
    parser.add_argument('--reads',type=Path)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--workdir',type=Path)
    parser.add_argument('--beluga-project',default='PRJNA1328577')
    parser.add_argument('--beluga-samples',type=int,default=2)
    parser.add_argument('--min-beluga-bases',type=int,default=50_000_000_000)
    parser.add_argument('--bear-runs',default=','.join(DEFAULT_BEAR))
    a=parser.parse_args()
    a.plan=a.plan.expanduser().resolve()
    if a.action=='references':
        if a.reference_root is None:raise ValueError('--reference-root required')
        for species in REFERENCE:
            fasta, contigs=fasta_info(a.reference_root.expanduser().resolve(),species)
            print(f'{species}: {len(contigs)} autosomes; FASTA={fasta}; first={contigs[0]}; last={contigs[-1]}')
        return
    if a.action=='discover':
        if a.plan.exists():raise ValueError(f'Frozen selection already exists: {a.plan}; move it explicitly before rediscovery')
        runs=tuple(a.bear_runs.split(','))
        if a.beluga_samples<1 or a.min_beluga_bases<1:raise ValueError('Invalid selection counts')
        plan=discover(a.beluga_project,a.beluga_samples,runs,a.min_beluga_bases)
        a.plan.parent.mkdir(parents=True,exist_ok=True)
        a.plan.write_text(json.dumps(plan,ensure_ascii=False,indent=2))
        report(plan);print('Frozen plan:',a.plan)
        return
    if not a.plan.is_file():raise ValueError('Run discover first to freeze accessions/MD5')
    plan=json.loads(a.plan.read_text())
    report(plan)
    if a.action=='report':return
    for key in ('project','reference_root','reads','output','workdir'):
        if getattr(a,key) is None:raise ValueError(f'--{key.replace("_","-")} is required for {a.action}')
        setattr(a,key,getattr(a,key).expanduser().resolve())
    if a.action=='download':
        remaining=sum(r['bytes'][i] for s in plan['samples'] for r in s['runs'] for i in (0,1)
                      if not (a.reads/s['species']/r['run']/f"{r['run']}_{i+1}.fastq.gz").exists())
        a.reads.mkdir(parents=True,exist_ok=True)
        free=shutil.disk_usage(a.reads).free
        # This checks FASTQ download capacity; intermediate WGS files need additional space.
        if free < remaining*1.1:
            raise RuntimeError(f'Not enough free space for downloads: free {free/1e9:.1f} GB, '
                               f'FASTQ planned {remaining/1e9:.1f} GB. GATK workdir needs additional storage.')
        for s in plan['samples']:
            for r in s['runs']:
                d=a.reads/s['species']/r['run']
                for i in (0,1):
                    downloader(r['urls'][i],d/f"{r['run']}_{i+1}.fastq.gz",r['md5'][i],r['bytes'][i])
    # Both 'download' and 'config' require all R1/R2 files with correct size+MD5.
    for s in plan['samples']:
        for r in s['runs']:
            d=a.reads/s['species']/r['run']
            for i in (0,1):
                f=d/f"{r['run']}_{i+1}.fastq.gz"
                if (not f.is_file() or f.stat().st_size != r['bytes'][i] or
                    md5_file(f) != r['md5'][i]):
                    raise RuntimeError(f'FASTQ missing, incomplete or MD5 mismatch: {f}')
    config=write_configuration(plan,a.project,a.reference_root,a.reads,a.output,a.workdir)
    print('Ready:',config)
    print('Run: bash run_pipeline.sh validate --config',config)
    print('Run: bash run_pipeline.sh plan --config',config)
    print('Run: bash run_pipeline.sh submit --config',config)


if __name__=='__main__':
    try:main()
    except (ValueError,RuntimeError,OSError,subprocess.CalledProcessError) as e:sys.exit('ERROR: '+str(e))
