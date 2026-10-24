#!/usr/bin/env python3
"""WGS-only Bash/Python/SLURM orchestration (BWA-MEM + GATK4)."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path
import yaml

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PROJECT/'scripts'))
from model import abspath, check_references, load_manifest  # noqa E402

STAGES=('trim','align','markdup','import','qc','haplotypecaller','genomicsdb','genotype','concat','multiqc')
RESOURCE_DEFAULTS={  # CPUs, MB RAM, hours; conservative starting values, NOT profiling-based
    'trim':(4,8000,4),'align':(8,20000,16),'markdup':(8,20000,12),
    'import':(4,12000,8),'qc':(4,8000,8),'haplotypecaller':(4,24000,24),
    'genomicsdb':(2,32000,18),'genotype':(2,16000,12),
    'concat':(2,8000,4),'multiqc':(2,8000,2),
}
SCRIPT={'trim':'trim.sh','align':'align_lane.sh','markdup':'markdup.sh',
        'import':'import_aligned.sh','qc':'bam_qc.sh','haplotypecaller':'haplotypecaller.sh',
        'genomicsdb':'genomicsdb.sh','genotype':'genotype.sh','concat':'concat.sh',
        'multiqc':'multiqc.sh'}
SAFE_ID=re.compile(r'^[A-Za-z0-9][A-Za-z0-9_.-]*$')


def shell(stage,*args):
    return shlex.join(['bash',str(PROJECT/'scripts'/'stages'/SCRIPT[stage]),*(str(x) for x in args)])


def config_load(config):
    cpath=Path(config).expanduser()
    if not cpath.is_absolute(): cpath=PROJECT/cpath
    cpath=cpath.resolve()
    cfg=yaml.safe_load(cpath.read_text(encoding='utf8'))
    if not isinstance(cfg,dict): raise ValueError('YAML must be mapping')
    if not isinstance(cfg.get('workdir'),str) or not cfg['workdir'].strip(): raise ValueError('Set workdir in config')
    work=Path(abspath(cfg['workdir'],PROJECT))
    if work==Path('/') or work==PROJECT or PROJECT in work.parents:
        raise ValueError('workdir must be outside the project directory and cannot be /')
    if not isinstance(cfg.get('references'),dict): raise ValueError('Set references in config')
    cfg['samples']=abspath(cfg['samples'],PROJECT)
    for sp,r in cfg['references'].items():
        r['fasta']=abspath(r['fasta'],PROJECT)
    check_references(cfg['references'],PROJECT,check_files=False)
    return cpath,cfg,work


def validate(cfg,check_files=True):
    check_references(cfg['references'],PROJECT,check_files=check_files)
    units,samples,groups=load_manifest(cfg['samples'],cfg['references'],PROJECT,check_files=check_files)
    p=cfg.get('preprocessing',{})
    if not isinstance(p,dict) or not isinstance(p.get('trim_fastq',False),bool):
        raise ValueError('preprocessing.trim_fastq must be true or false')
    d=cfg.get('duplicate_marking',{})
    if not isinstance(d,dict) or int(d.get('optical_distance',0))<0: raise ValueError('Invalid duplicate_marking')
    g=cfg.get('gatk',{})
    for k in ('hc_heap_gb','db_heap_gb','genotype_heap_gb','db_batch_size','db_readers','hc_threads'):
        if not isinstance(g.get(k,1),int) or g.get(k,1)<1:
            raise ValueError(f'gatk.{k} must be positive integer')
    if check_files:
        for u in units:
            if u.input_type in ('bam','cram'):
                f=Path(u.source)
                indexpaths=([Path(str(f)+'.bai'),f.with_suffix('.bai'),Path(str(f)+'.csi')]
                            if u.input_type=='bam' else [Path(str(f)+'.crai'),f.with_suffix('.crai')])
                if not any(x.is_file() and x.stat().st_size for x in indexpaths):
                    raise ValueError(f'Missing index for external {u.input_type.upper()} {f}; index outside the pipeline')
    return units,samples,groups


def task(stage, name, args, outputs, deps=(), inputs=(), extra=''):
    if not SAFE_ID.fullmatch(name): raise ValueError(f'Unsafe job ID {name}')
    body='set -euo pipefail\n'+shell(stage,*args)+'\n'
    signature=hashlib.sha256(body.encode()+ (PROJECT/'scripts'/'stages'/SCRIPT[stage]).read_bytes()+str(extra).encode())
    if stage=='import':
        signature.update((PROJECT/'scripts'/'validate_alignment.py').read_bytes())
    return {'id':name,'stage':stage,'deps':list(dict.fromkeys(deps)),
            'inputs':list(map(str,inputs)),'outputs':list(map(str,outputs)),
            'body':body,'signature':signature.hexdigest()}


def make_tasks(cfg,work,units,samples,groups):
    refs=cfg['references']
    trim=cfg.get('preprocessing',{}).get('trim_fastq',False)
    optical=int(cfg.get('duplicate_marking',{}).get('optical_distance',0))
    gc=cfg.get('gatk',{})
    hc_threads=gc.get('hc_threads',4); hc_heap=gc.get('hc_heap_gb',16)
    db_heap=gc.get('db_heap_gb',24); gt_heap=gc.get('genotype_heap_gb',12)
    db_batch=gc.get('db_batch_size',50); readers=gc.get('db_readers',1)
    # Verify that fixed resource requests are not below hardcoded active threads/heap.
    overrides=cfg.get('slurm',{}).get('resources',{})
    hc_cpus=int(overrides.get('haplotypecaller',{}).get('cpus',4))
    if hc_threads>hc_cpus: raise ValueError('gatk.hc_threads exceeds Slurm cpus for haplotypecaller')
    align_cpus=int(overrides.get('align',{}).get('cpus',8))
    if align_cpus<8:
        raise ValueError('Alignment script requires at least 8 SLURM cpus')
    if int(overrides.get('markdup',{}).get('cpus',8))<8:
        raise ValueError('Markdup script requires at least 8 SLURM cpus')
    if int(overrides.get('haplotypecaller',{}).get('mem_mb',24000))<=hc_heap*1024:
        raise ValueError('HaplotypeCaller Java heap must be less than SLURM memory')
    if int(overrides.get('genomicsdb',{}).get('mem_mb',32000))<=db_heap*1024:
        raise ValueError('GenomicsDB Java heap must be less than SLURM memory')
    if int(overrides.get('genotype',{}).get('mem_mb',16000))<=gt_heap*1024:
        raise ValueError('GenotypeGVCFs Java heap must be less than SLURM memory')
    ts=[]; seen=set()
    def add(stage,id,args,outputs,deps=(),inputs=(),extra=''):
        if id in seen: raise ValueError(f'Duplicate job ID {id}')
        seen.add(id)
        t=task(stage,id,args,outputs,deps,inputs,extra)
        ts.append(t)
        return id
    lanes={}
    bysample={}
    for u in units:
        if u.input_type!='fastq': continue
        uid=f'{u.sample}.{u.library}.{u.lane}'
        ref=refs[u.species]['fasta']
        r1,r2=u.r1,u.r2
        deps=[]; inputs=[]
        if trim:
            pre=work/'trimmed'/u.sample
            r1=str(pre/f'{u.library}.{u.lane}_R1.fastq.gz')
            r2=str(pre/f'{u.library}.{u.lane}_R2.fastq.gz')
            q=work/'qc'/'fastp'/u.sample
            trimid=add('trim',f'trim.{uid}',[u.r1,u.r2,r1,r2,q/f'{u.library}.{u.lane}.json',q/f'{u.library}.{u.lane}.html'],
                       [r1,r2,q/f'{u.library}.{u.lane}.json',q/f'{u.library}.{u.lane}.html'], inputs=[u.r1,u.r2])
            deps=[trimid]
        else: inputs=[u.r1,u.r2]
        lane=work/'lanes'/u.sample/f'{u.library}.{u.lane}.bam'
        rid=add('align',f'align.{uid}',[ref,r1,r2,lane,u.sample,u.library,u.lane,align_cpus],
                [lane],deps=deps,inputs=[*inputs,ref,str(ref)+'.fai',*(str(ref)+'.'+x for x in ('amb','ann','bwt','pac','sa'))])
        lanes.setdefault(u.sample,[]).append((lane,rid))
    for sample,u in samples.items():
        ref=refs[u.species]['fasta']
        if u.input_type=='fastq':
            bam=work/'bam'/f'{sample}.bam'
            metrics=work/'qc'/'duplicates'/f'{sample}.json'
            lane_files, lane_deps=zip(*lanes[sample])
            deps=[*lane_deps]
            prep=add('markdup',f'markdup.{sample}',[bam,metrics,optical,8,*lane_files],
                     [bam,str(bam)+'.bai',metrics],deps=deps)
        elif u.input_type=='bam':
            bam=u.source
            report=work/'qc'/'input_validation'/f'{sample}.idxstats.tsv'
            prep=add('import',f'import.{sample}', ['bam',u.source,report,sample,ref],
                     [report],inputs=[u.source,ref,str(ref)+'.fai',refs[u.species]['dict']])
        else:
            bam=work/'bam'/f'{sample}.bam'
            prep=add('import',f'import.{sample}', ['cram',u.source,bam,sample,ref],
                     [bam,str(bam)+'.bai'],inputs=[u.source,ref,str(ref)+'.fai',refs[u.species]['dict']])
        bysample[sample]=(str(bam),prep)
        qc=work/'qc'/'bam'/sample
        qc_out=[str(qc)+e for e in ('.flagstat.txt','.stats.txt','.mosdepth.summary.txt','.regions.bed.gz','.thresholds.bed.gz')]
        add('qc',f'qc.{sample}',[bam,qc,4],qc_out,deps=[prep])
    maps={}
    for species,sample_ids in sorted(groups.items()):
        ref=refs[species]['fasta']
        for contig in refs[species]['contigs']:
            gvcfs={};hc_ids=[]
            for sample in sample_ids:
                bam,prep=bysample[sample]
                gvcf=work/'gvcf'/species/sample/f'{contig}.g.vcf.gz'
                tmp=work/'tmp'/'hc'/species/sample/contig
                ht=add('haplotypecaller',f'hc.{species}.{sample}.{contig}',
                       [ref,bam,contig,gvcf,hc_heap,hc_threads,tmp],
                       [gvcf],deps=[prep],inputs=[ref,str(ref)+'.fai',refs[species]['dict']])
                gvcfs[sample]=str(gvcf)
                hc_ids.append(ht)
            samplemap=''.join(f'{s}\t{gvcfs[s]}\n' for s in sorted(sample_ids))
            mapfile=work/'jobs'/'sample_maps'/species/f'{contig}.tsv'
            maps[str(mapfile)]=samplemap
            db=work/'genomicsdb'/species/contig
            d=add('genomicsdb',f'db.{species}.{contig}',
                  [ref,contig,mapfile,db,db_heap,db_batch,readers,work/'tmp'/'genomicsdb'/species/contig],
                  [db/'.completed'],deps=hc_ids,extra=samplemap)
            vcf=work/'joint'/species/f'{contig}.vcf.gz'
            v=add('genotype',f'genotype.{species}.{contig}',
                  [ref,contig,db,vcf,gt_heap,work/'tmp'/'genotype'/species/contig],
                  [vcf],deps=[d])
        vcfs=[work/'joint'/species/f'{c}.vcf.gz' for c in refs[species]['contigs']]
        out=work/'vcf'/species/'cohort.raw.vcf.gz'
        add('concat',f'concat.{species}',[out,*vcfs],[out,str(out)+'.csi'],
            deps=[f'genotype.{species}.{c}' for c in refs[species]['contigs']])
    qcout=work/'qc'/'multiqc'/'multiqc_report.html'
    add('multiqc','multiqc',[work/'qc',qcout.parent],[qcout],
        deps=[t['id'] for t in ts if t['stage'] in ('trim','qc')])
    ts.sort(key=lambda t:STAGES.index(t['stage']))
    return ts,maps


def done(t,work):
    marker=work/'checkpoints'/f"{t['id']}.json"
    try:
        record=json.loads(marker.read_text(encoding='utf8'))
        if record['signature']!=t['signature']: return False
        mtime=marker.stat().st_mtime
        if any(not Path(f).is_file() or Path(f).stat().st_size==0 for f in t['outputs']): return False
        if t['stage'] in ('haplotypecaller','genotype'):
            vcf=t['outputs'][0]
            if not any(Path(vcf+ext).is_file() and Path(vcf+ext).stat().st_size for ext in ('.tbi','.idx','.csi')): return False
        if any(not Path(f).is_file() or Path(f).stat().st_mtime>mtime for f in t['inputs']): return False
        return True
    except (OSError,ValueError,KeyError): return False


def pending(ts,work,force=False):
    needs={}
    for t in ts:
        needs[t['id']]=force or not done(t,work) or any(needs[d] for d in t['deps'])
    return needs


def save_plan(ts,maps,cfg_path,work):
    (work/'jobs'/'scripts').mkdir(parents=True,exist_ok=True)
    (work/'jobs'/'logs').mkdir(parents=True,exist_ok=True)
    (work/'checkpoints').mkdir(parents=True,exist_ok=True)
    for t in ts:
        for out in t['outputs']: Path(out).parent.mkdir(parents=True,exist_ok=True)
        p=work/'jobs'/'scripts'/f"{t['id']}.sh"
        p.write_text(t['body'],encoding='utf8')
        p.chmod(0o750)
    for filename,content in maps.items():
        p=Path(filename)
        p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text(content,encoding='utf8')
    path=work/'pipeline_plan.json'
    tmp=work/'pipeline_plan.json.tmp'
    tmp.write_text(json.dumps({'root':str(PROJECT),'workdir':str(work),
        'config':str(cfg_path),'tasks':ts},indent=2,ensure_ascii=False),encoding='utf8')
    tmp.replace(path)
    return path


def execute(plan,id):
    saved=json.loads(Path(plan).read_text(encoding='utf8'))
    work=Path(saved['workdir'])
    t=next((x for x in saved['tasks'] if x['id']==id),None)
    if not t: raise ValueError(f'Unknown task {id}')
    script=work/'jobs'/'scripts'/f'{id}.sh'
    if script.read_text(encoding='utf8')!=t['body']: raise ValueError('Task script differs from frozen plan')
    marker=work/'checkpoints'/f'{id}.json'
    marker.unlink(missing_ok=True)
    for f in t['outputs']: Path(f).unlink(missing_ok=True)
    print(f'RUN {id}',flush=True)
    subprocess.run(['bash',str(script)],cwd=PROJECT,check=True)
    if any(not Path(f).is_file() or Path(f).stat().st_size==0 for f in t['outputs']):
        raise RuntimeError(f'{id}: finished without expected nonempty outputs')
    if t['stage'] in ('haplotypecaller','genotype'):
        vcf=t['outputs'][0]
        if not any(Path(vcf+ext).is_file() and Path(vcf+ext).stat().st_size for ext in ('.tbi','.idx','.csi')):
            raise RuntimeError(f'{id}: missing indexed GATK VCF')
    tmp=marker.with_suffix('.tmp')
    tmp.write_text(json.dumps({'signature':t['signature'],'finished':time.time()}))
    tmp.replace(marker)
    print(f'DONE {id}',flush=True)


def active_jobs(work):
    record=work/'jobs'/'submission.json'
    if not record.exists() or not shutil.which('squeue'): return []
    ids=[str(n) for n in json.loads(record.read_text()).get('job_ids',{}).values()]
    if not ids: return []
    result=subprocess.run(['squeue','-h','-j',','.join(ids),'-o','%A'],capture_output=True,text=True,check=True)
    return sorted(set(result.stdout.split()))


def submit(ts,needs,work,plan,cfg):
    if not shutil.which('sbatch'): raise RuntimeError('sbatch not found (activate SLURM client environment)')
    sl=cfg.get('slurm',{})
    max_parallel=int(sl.get('max_array_parallel',4))
    if max_parallel<1: raise ValueError('max_array_parallel must be >=1')
    override=sl.get('resources',{})
    jobs={}
    for stage in STAGES:
        selection=[t for t in ts if needs[t['id']] and t['stage']==stage]
        if not selection: continue
        depends={jobs[prev['stage']] for t in selection for d in t['deps']
                 for prev in ts if prev['id']==d and needs[d] and prev['stage']!=stage}
        file=work/'jobs'/f'array_{stage}.txt'
        file.write_text('\n'.join(t['id'] for t in selection)+'\n')
        cpus,mem,hours=RESOURCE_DEFAULTS[stage]
        res=override.get(stage,{})
        cpus=int(res.get('cpus',cpus));mem=int(res.get('mem_mb',mem));hours=int(res.get('hours',hours))
        if min(cpus,mem,hours)<1: raise ValueError(f'{stage}: invalid SLURM resources')
        args=['sbatch','--parsable',f'--job-name=WGS_{stage}',
              f'--cpus-per-task={cpus}',f'--mem={mem}M',f'--time={hours}:00:00',
              f'--array=0-{len(selection)-1}%{max_parallel}',
              f'--output={work}/jobs/logs/{stage}_%A_%a.out',
              f'--error={work}/jobs/logs/{stage}_%A_%a.err',
              f'--export=ALL,PIPELINE_PYTHON={sys.executable}']
        if sl.get('partition'): args.append(f"--partition={sl['partition']}")
        if sl.get('account'): args.append(f"--account={sl['account']}")
        if sl.get('qos'): args.append(f"--qos={sl['qos']}")
        if depends: args.append('--dependency=afterok:'+':'.join(sorted(depends)))
        args.extend([str(PROJECT/'scripts'/'run_array.sh'),str(plan),str(file)])
        result=subprocess.run(args,text=True,capture_output=True)
        if result.returncode: raise RuntimeError(f'sbatch failed for {stage}: {result.stderr.strip()}')
        jobid=result.stdout.strip().split(';')[0]
        if not re.fullmatch(r'\d+',jobid): raise RuntimeError(f'Cannot parse sbatch jobid: {result.stdout}')
        jobs[stage]=jobid
        (work/'jobs'/'submission.json').write_text(json.dumps({'job_ids':jobs,'plan':str(plan)},indent=2))
        print(f'{stage:16s} {len(selection):5d} elements -> SLURM {jobid}',flush=True)
    return jobs


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    sp=ap.add_subparsers(dest='action',required=True)
    for cmdname in ('validate','plan','local','submit','status'):
        p=sp.add_parser(cmdname)
        p.add_argument('--config',default='config/config.yaml')
        if cmdname in ('plan','local','submit'):
            p.add_argument('--force',action='store_true',help='Recompute all tasks; never while jobs running')
        if cmdname=='plan':
            p.add_argument('--allow-missing-inputs',action='store_true',help='Preview before data/reference download')
    p=sp.add_parser('run-task');p.add_argument('--plan',required=True);p.add_argument('--task-id',required=True)
    args=ap.parse_args()
    if args.action=='run-task': execute(args.plan,args.task_id);return
    path,cfg,work=config_load(args.config)
    if args.action=='status':
        print('Active SLURM jobs:',active_jobs(work))
        return
    allow_missing=bool(getattr(args,'allow_missing_inputs',False))
    units,samples,groups=validate(cfg,check_files=not allow_missing)
    ts,maps=make_tasks(cfg,work,units,samples,groups)
    needs=pending(ts,work,force=getattr(args,'force',False))
    print('Samples:',len(samples),'Species:',', '.join(groups))
    for stage in STAGES:
        print(f'{stage:16s} pending {sum(needs[t["id"]] for t in ts if t["stage"]==stage):5d} / {sum(t["stage"]==stage for t in ts):5d}')
    if args.action in ('validate','plan'):
        print('Validation/plan only: no workdir files written and no tasks submitted.')
        return
    active=active_jobs(work)
    if active: raise RuntimeError(f'Existing jobs still in queue: {active}; avoid overwriting active plans')
    if not any(needs.values()): print('All expected products already completed');return
    plan=save_plan(ts,maps,path,work)
    if args.action=='local':
        for t in ts:
            if needs[t['id']]: execute(plan,t['id'])
    elif args.action=='submit': submit(ts,needs,work,plan,cfg)

if __name__=='__main__':
    try: main()
    except (ValueError,KeyError,RuntimeError,OSError,subprocess.CalledProcessError) as e:
        sys.exit(f'ERROR: {e}')
