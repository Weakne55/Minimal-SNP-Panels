import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT/'setup_genomes.sh'
CHECKER = ROOT/'reference_summary.py'


class SetupTests(unittest.TestCase):
    def test_invalid_fai_and_dict_detected(self):
        with tempfile.TemporaryDirectory() as d:
            ref = Path(d)/'example.fa'
            ref.write_text('>chrA\nACGT\n>chrB\nAA\n')
            Path(str(ref)+'.fai').write_text('chrA\t4\t6\t4\t5\nchrB\t2\t17\t2\t3\n')
            ref.with_suffix('.dict').write_text('@HD\tVN:1.6\n@SQ\tSN:chrA\tLN:4\n@SQ\tSN:chrB\tLN:3\n')
            for ext in ('amb','ann','bwt','pac','sa'):
                Path(str(ref)+'.'+ext).write_text('index')
            result = subprocess.run([sys.executable,str(CHECKER),'check',str(ref)],text=True,capture_output=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('disagree',result.stderr)
            ref.with_suffix('.dict').write_text('@HD\tVN:1.6\n@SQ\tSN:chrA\tLN:4\n@SQ\tSN:chrB\tLN:2\n')
            ok = subprocess.run([sys.executable,str(CHECKER),'check',str(ref)],text=True,capture_output=True)
            self.assertEqual(ok.returncode,0,ok.stderr)

    def test_index_one_invokes_indexers_only_once(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            ref = d/'toy.fna'
            ref.write_text('>chr1\nACGT\n>chr2\nATCG\n')
            bin_dir = d/'bin'
            bin_dir.mkdir()
            calls = d/'calls.log'
            programs = {
                'samtools': '''#!/bin/bash
set -eu
echo samtools >> "$CALLS_FILE"
[[ "$1" == faidx ]] || exit 17
printf 'chr1\\t4\\t6\\t4\\t5\\nchr2\\t4\\t17\\t4\\t5\\n' > "$2.fai"
''',
                'bwa': '''#!/bin/bash
set -eu
echo bwa >> "$CALLS_FILE"
[[ "$1" == index ]] || exit 17
ref="$4"
for ext in amb ann bwt pac sa; do echo content > "$ref.$ext"; done
''',
                'gatk': '''#!/bin/bash
set -eu
echo gatk >> "$CALLS_FILE"
[[ "$1" == CreateSequenceDictionary ]] || exit 17
while (( $# )); do
  if [[ "$1" == -O ]]; then shift; dict="$1"; fi
  shift
done
printf '@HD\\tVN:1.6\\n@SQ\\tSN:chr1\\tLN:4\\n@SQ\\tSN:chr2\\tLN:4\\n' > "$dict"
''',
            }
            for name, script in programs.items():
                path = bin_dir/name
                path.write_text(script)
                path.chmod(0o755)
            env={**os.environ,'PATH':str(bin_dir)+os.pathsep+os.environ['PATH'],'CALLS_FILE':str(calls)}
            args=['bash',str(SCRIPT),'index-one',str(ref)]
            first=subprocess.run(args,env=env,text=True,capture_output=True)
            self.assertEqual(first.returncode,0,first.stderr)
            self.assertEqual(set(calls.read_text().splitlines()),{'samtools','gatk','bwa'})
            before=calls.read_text()
            again=subprocess.run(args,env=env,text=True,capture_output=True)
            self.assertEqual(again.returncode,0,again.stderr)
            self.assertEqual(calls.read_text(),before)

    def test_missing_assemblies_fails_fast(self):
        with tempfile.TemporaryDirectory() as d:
            x=subprocess.run(['bash',str(SCRIPT),'index','--root',d],capture_output=True,text=True)
            self.assertNotEqual(x.returncode,0)
            self.assertIn('Both assemblies',x.stderr)


if __name__=='__main__':
    unittest.main()


class NcbiDownloadTests(unittest.TestCase):
    def test_both_ncbi_accessions_requested(self):
        from zipfile import ZipFile
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            fake_bin=d/'fake-bin'; fake_bin.mkdir()
            fake_datasets=fake_bin/'datasets'
            fake_datasets.write_text('''#!/usr/bin/env python3
import sys,zipfile
args=sys.argv[1:]
assert args[:3] == ['download','genome','accession'], args
accession=args[3]
assert accession in ['GCA_059928025.1','GCA_062004589.1'], accession
assert args[args.index('--include')+1]=='genome,seq-report'
out=args[args.index('--filename')+1]
with zipfile.ZipFile(out,'w') as z:
    z.writestr(f'ncbi_dataset/data/{accession}/{accession}_toy_genomic.fna', '>chr1\\nACGT\\n')
''')
            fake_datasets.chmod(0o755)
            env={**os.environ, 'PATH':str(fake_bin)+os.pathsep+os.environ['PATH']}
            cp=subprocess.run(['bash',str(SCRIPT),'download','--root',str(d/'genomes')],
                              env=env,capture_output=True,text=True)
            self.assertEqual(cp.returncode,0,cp.stderr)
            for species,acc in [('Delphinapterus_leucas','GCA_059928025.1'),
                                ('Ursus_maritimus','GCA_062004589.1')]:
                fasta = d/'genomes'/f'{species}_{acc}'/'ncbi_dataset'/'data'/acc/f'{acc}_toy_genomic.fna'
                self.assertTrue(fasta.is_file(), str(fasta))
                self.assertIn('>chr1', fasta.read_text())
            # On a second download neither genome should be fetched again.
            cp2=subprocess.run(['bash',str(SCRIPT),'download','--root',str(d/'genomes')],
                               env=env,capture_output=True,text=True)
            self.assertEqual(cp2.returncode,0,cp2.stderr)
            self.assertIn('already present',cp2.stderr)
