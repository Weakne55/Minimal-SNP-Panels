import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SOURCE=Path(__file__).resolve().parents[1]/'prepare_full_wgs.py'
spec=importlib.util.spec_from_file_location('prep',SOURCE)
p=importlib.util.module_from_spec(spec)
spec.loader.exec_module(p)


def row(run,sample,species='Delphinapterus_leucas',bases=65_000_000_000,lib='LIB1'):
    return dict(run_accession=run,sample_accession=sample,scientific_name=species.replace('_',' '),
                library_strategy='WGS',library_layout='PAIRED',instrument_platform='ILLUMINA',
                base_count=str(bases),instrument_model='NovaSeq 6000',library_name=lib,
                experiment_accession='SRX123',study_accession='PRJNA1328577',
                fastq_ftp=f'ftp.sra.ebi.ac.uk/vol1/fastq/{run}/{run}_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/{run}/{run}_2.fastq.gz',
                fastq_md5='a'*32+';'+'b'*32,fastq_bytes='2500;2300')


class FullWgsTests(unittest.TestCase):
    def test_reject_lc_wgs_and_unpaired(self):
        x=row('SRR123','SAM1');x['library_strategy']='RNA-Seq'
        self.assertIsNone(p.normalize(x,'Delphinapterus_leucas'))
        x=row('SRR123','SAM1');x['library_layout']='SINGLE'
        self.assertIsNone(p.normalize(x,'Delphinapterus_leucas'))
        x=row('SRR123','SAM1');x['fastq_md5']=''
        self.assertIsNone(p.normalize(x,'Delphinapterus_leucas'))

    def test_choose_distinct_beluga_samples_with_sufficient_bases(self):
        records=[row('SRR101','SAM1'),row('SRR102','SAM2'),row('SRR103','SAM3',bases=10_000_000_000),
                 row('SRR104','SAM1',bases=60_000_000_000)]
        groups=p.choose_beluga(records,2,50_000_000_000)
        self.assertEqual(len(groups),2)
        self.assertEqual(len({g[0]['biosample'] for g in groups}),2)
        self.assertTrue(any(len(g)==2 for g in groups))

    def test_reject_multiple_libraries_of_one_biosample(self):
        rows=[row('SRR101','SAM1',lib='A'),row('SRR102','SAM1',lib='B')]
        with self.assertRaisesRegex(ValueError,'Only 0'):
            p.choose_beluga(rows,1,50_000_000_000)

    def test_autosomes_via_ncbi_report_not_by_length(self):
        with tempfile.TemporaryDirectory() as folder:
            f=Path(folder)
            sp='Ursus_maritimus'
            acc, directory, expected=p.REFERENCE[sp]
            d=f/directory/'ncbi_dataset'/'data'/acc
            d.mkdir(parents=True)
            fasta=d/f'{acc}_genomic.fna';fasta.write_text('>test\nA\n')
            rows=[]
            seqs=[]
            # Deliberately X between chromosome 1 and 2.
            for i in range(1,expected+1):
                name=f'CM{i:06d}.1'
                rows.append(dict(chrName=str(i),role='assembled-molecule',assignedMoleculeLocationType='Chromosome',
                                 genbankAccession=name,length=1_000_000+i))
                seqs.append((name,1_000_000+i))
            rows.insert(1,dict(chrName='X',role='assembled-molecule',assignedMoleculeLocationType='Chromosome',
                               genbankAccession='CM999999.1',length=12_000_000))
            seqs.insert(1,('CM999999.1',12_000_000))
            (d/'sequence_report.jsonl').write_text('\n'.join(json.dumps(r) for r in rows)+'\n')
            Path(str(fasta)+'.fai').write_text(''.join(f'{n}\t{size}\t0\t80\t81\n' for n,size in seqs))
            _,selected=p.fasta_info(f,sp)
            self.assertEqual(len(selected),expected)
            self.assertNotIn('CM999999.1',selected)
            self.assertEqual(selected[0],'CM000001.1')
            self.assertEqual(selected[1],'CM000002.1')

    def test_fails_without_chromosome_annotations(self):
        with tempfile.TemporaryDirectory() as folder:
            f=Path(folder);sp='Delphinapterus_leucas'
            acc, directory, expected=p.REFERENCE[sp]
            d=f/directory/'ncbi_dataset'/'data'/acc;d.mkdir(parents=True)
            fasta=d/f'{acc}_genomic.fna';fasta.write_text('>A\nA\n')
            Path(str(fasta)+'.fai').write_text('A\t1\t3\t1\t2\n')
            (d/'sequence_report.jsonl').write_text(json.dumps(dict(chrName='Un',genbankAccession='A'))+'\n')
            with self.assertRaisesRegex(ValueError,'expected 21'):
                p.fasta_info(f,sp)


if __name__=='__main__':unittest.main()
