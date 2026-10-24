import json
import shlex
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from model import check_references,load_manifest
import pipeline


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base=Path(self.tmp.name)
        self.root=self.base/'project'
        self.root.mkdir()
        self.ref=self.root/'test.fa'
        self.ref.write_text('>chr1\nACGTACGT\n>chr2\nTTTT\n')
        (self.root/'test.fa.fai').write_text('chr1\t8\t6\t8\t9\nchr2\t4\t21\t4\t5\n')
        (self.root/'test.dict').write_text('@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:8\n@SQ\tSN:chr2\tLN:4\n')
        for ext in ('amb','ann','bwt','pac','sa'):
            (self.root/f'test.fa.{ext}').write_text('index')
        self.refs={'Species_a':{'fasta':str(self.ref),'contigs':['chr1','chr2']}}
        self.tab=self.root/'samples.tsv'
        self.tab.write_text('sample\tspecies\tinput_type\tlibrary\tlane\tr1\tr2\tbam\tcram\n'
                            'A1\tSpecies_a\tfastq\tL1\tL001\t/data/A1_R1.fastq.gz\t/data/A1_R2.fastq.gz\t\t\n'
                            'A1\tSpecies_a\tfastq\tL1\tL002\t/data/A1L2_R1.fastq.gz\t/data/A1L2_R2.fastq.gz\t\t\n'
                            'A2\tSpecies_a\tfastq\tL1\tL001\t/data/A2_R1.fastq.gz\t/data/A2_R2.fastq.gz\t\t\n')
        self.work=self.base/'run'
        self.cfg={'references':self.refs,'samples':str(self.tab),'workdir':str(self.work),
                  'preprocessing':{'trim_fastq':False},'gatk':{'hc_threads':4},'slurm':{}}

    def build(self):
        check_references(self.refs,self.root)
        units,samples,groups=load_manifest(str(self.tab),self.refs,self.root,check_files=False)
        return pipeline.make_tasks(self.cfg,self.work,units,samples,groups)

    def test_existing_reference_is_never_rewritten(self):
        before=self.ref.stat().st_mtime_ns
        check_references(self.refs,self.root)
        self.assertEqual(self.ref.stat().st_mtime_ns,before)
        self.assertNotIn('ref',pipeline.STAGES)
        self.assertFalse(any('prepare_reference' in x for x in pipeline.SCRIPT.values()))

    def test_missing_reference_index(self):
        (self.root/'test.fa.sa').unlink()
        with self.assertRaisesRegex(ValueError,'Missing reference/index'):
            check_references(self.refs,self.root)

    def test_contig_order(self):
        self.refs['Species_a']['contigs']=['chr2','chr1']
        with self.assertRaisesRegex(ValueError,'order must match'):
            check_references(self.refs,self.root)

    def test_manifest_multilane_one_library(self):
        units,samples,groups=load_manifest(str(self.tab),self.refs,self.root,check_files=False)
        self.assertEqual(len(units),3)
        self.assertEqual(len(samples),2)
        self.assertEqual(groups['Species_a'],['A1','A2'])

    def test_manifest_no_multiple_libraries(self):
        self.tab.write_text(self.tab.read_text().replace('A1\tSpecies_a\tfastq\tL1\tL002','A1\tSpecies_a\tfastq\tL2\tL002'))
        with self.assertRaisesRegex(ValueError,'multiple independent libraries'):
            load_manifest(str(self.tab),self.refs,self.root,check_files=False)

    def test_no_lowpass(self):
        self.tab.write_text(self.tab.read_text().replace('sample\tspecies','coverage_type\tsample\tspecies'))
        with self.assertRaisesRegex(ValueError,'WGS only'):
            load_manifest(str(self.tab),self.refs,self.root,check_files=False)

    def test_hc_gvcf_and_joint_DAG(self):
        ts,maps=self.build()
        named={t['id']:t for t in ts}
        self.assertEqual(sum(t['stage']=='haplotypecaller' for t in ts),4)
        self.assertEqual(sum(t['stage']=='genomicsdb' for t in ts),2)
        self.assertIn('hc.Species_a.A1.chr1',named['db.Species_a.chr1']['deps'])
        self.assertIn('hc.Species_a.A2.chr1',named['db.Species_a.chr1']['deps'])
        self.assertIn('db.Species_a.chr1',named['genotype.Species_a.chr1']['deps'])
        self.assertEqual(len(named['concat.Species_a']['deps']),2)
        self.assertEqual(len(maps),2)
        self.assertTrue(all('A1\t' in s and 'A2\t' in s for s in maps.values()))
        self.assertNotIn('LowQual',named['concat.Species_a']['body'])

    def test_trim_disabled_by_default(self):
        ts,_=self.build()
        self.assertEqual(sum(t['stage']=='trim' for t in ts),0)
        a=next(t for t in ts if t['stage']=='align')
        self.assertIn('A1_R1.fastq.gz',a['body'])

    def test_trim_enabled(self):
        self.cfg['preprocessing']['trim_fastq']=True
        ts,_=self.build()
        self.assertEqual(sum(t['stage']=='trim' for t in ts),3)
        a=next(t for t in ts if t['stage']=='align')
        self.assertEqual(len(a['deps']),1)

    def test_sample_addition_changes_cohort_import_signature(self):
        old,_=self.build()
        oldsig=next(t['signature'] for t in old if t['id']=='db.Species_a.chr1')
        self.tab.write_text(self.tab.read_text()+'A3\tSpecies_a\tfastq\tL1\tL001\t/data/X1\t/data/X2\t\t\n')
        new,_=self.build()
        newsig=next(t['signature'] for t in new if t['id']=='db.Species_a.chr1')
        self.assertNotEqual(oldsig,newsig)

    def test_completion_recomputed_for_changed_upstream(self):
        a=pipeline.task('trim','trim.test',['R1','R2','O1','O2','J','H'],[self.work/'a.txt'])
        b=pipeline.task('align','align.test',['ref','R1','R2','O','S','L','01',8],
                        [self.work/'b.txt'],deps=[a['id']])
        self.assertTrue(all(pipeline.pending([a,b],self.work).values()))
        (self.work/'checkpoints').mkdir(parents=True)
        for t in (a,b):
            Path(t['outputs'][0]).write_text('data')
            (self.work/'checkpoints'/f"{t['id']}.json").write_text(json.dumps({'signature':t['signature']}))
        self.assertFalse(any(pipeline.pending([a,b],self.work).values()))
        a['signature']='changed'
        self.assertTrue(all(pipeline.pending([a,b],self.work).values()))

    def test_imported_bam_used_in_place(self):
        self.tab.write_text('sample\tspecies\tinput_type\tlibrary\tlane\tr1\tr2\tbam\tcram\n'
                            'A1\tSpecies_a\tbam\tL1\tL001\t\t\t/data/input.bam\t\n')
        ts,_=self.build(); named={t['id']:t for t in ts}
        self.assertIn('import.A1',named)
        self.assertNotIn('markdup.A1',named)
        self.assertIn('/data/input.bam',named['hc.Species_a.A1.chr1']['body'])

    def test_task_executor_creates_checkpoint(self):
        target=self.work/'hello.txt'
        t=pipeline.task('trim','trim.hello',['in1','in2','out1','out2','json','html'],[target])
        t['body']="set -euo pipefail\nprintf 'ready\\n' > "+shlex.quote(str(target))+"\n"
        plan=pipeline.save_plan([t],{},'config.yaml',self.work)
        pipeline.execute(str(plan),'trim.hello')
        self.assertTrue(pipeline.done(t,self.work))
        self.assertEqual(target.read_text().strip(),'ready')

    def test_failed_task_does_not_leave_checkpoint(self):
        target=self.work/'bad.txt'
        t=pipeline.task('trim','trim.bad',['in1','in2','out1','out2','json','html'],[target])
        t['body']='set -euo pipefail\nfalse\n'
        plan=pipeline.save_plan([t],{},'config.yaml',self.work)
        with self.assertRaises(subprocess.CalledProcessError):
            pipeline.execute(str(plan),'trim.bad')
        self.assertFalse((self.work/'checkpoints'/'trim.bad.json').exists())

    def test_slurm_stage_dependencies(self):
        ts,maps=self.build()
        needs={t['id']:True for t in ts}
        plan=pipeline.save_plan(ts,maps,'config.yaml',self.work)
        issued=[]
        def simulate(args,**kwargs):
            issued.append(args)
            return subprocess.CompletedProcess(args,0,stdout=str(100+len(issued))+'\n',stderr='')
        with patch.object(pipeline.shutil,'which',return_value='/usr/bin/sbatch'),\
             patch.object(pipeline.subprocess,'run',side_effect=simulate):
            pipeline.submit(ts,needs,self.work,plan,self.cfg)
        self.assertGreater(len(issued),5)
        self.assertTrue(any('--dependency=afterok:' in arg for row in issued for arg in row))
        self.assertTrue(any('--array=' in arg for row in issued for arg in row))


if __name__=='__main__': unittest.main()
