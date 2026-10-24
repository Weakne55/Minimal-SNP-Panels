#!/usr/bin/env python3
"""List existing reference contigs without indexing or modifying the FASTA."""
import argparse
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('reference_fai', help='Path to existing reference.fa.fai')
a=p.parse_args()
f=Path(a.reference_fai)
if not f.is_file(): p.error('Missing .fai; this command never generates reference indexes')
print('contig\tlength')
for line in f.read_text().splitlines():
    if line.strip():
        cols=line.split('\t')
        print(f'{cols[0]}\t{cols[1]}')
