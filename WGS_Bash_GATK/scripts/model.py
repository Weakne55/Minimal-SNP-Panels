"""Input validation and DAG source models. No third-party imports except caller's PyYAML."""
from __future__ import annotations
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from collections import defaultdict

ID = re.compile(r'^[A-Za-z0-9][A-Za-z0-9_.-]*$')
COLUMNS = ('sample','species','input_type','library','lane','r1','r2','bam','cram')

@dataclass(frozen=True)
class Unit:
    sample: str
    species: str
    input_type: str
    library: str
    lane: str
    r1: str = ''
    r2: str = ''
    bam: str = ''
    cram: str = ''

    @property
    def key(self):
        return (self.sample, self.library, self.lane)
    @property
    def source(self):
        return self.bam or self.cram


def abspath(path: str, root: Path) -> str:
    p = Path(path).expanduser()
    return str((p if p.is_absolute() else root/p).resolve())


def load_manifest(path: str, refs: dict, root: Path, check_files: bool = True):
    f = Path(path)
    if not f.is_file():
        raise ValueError(f'Missing samples TSV: {f}')
    with f.open(newline='', encoding='utf8') as fd:
        reader = csv.DictReader(fd, delimiter='\t')
        if reader.fieldnames is None or not set(COLUMNS).issubset(reader.fieldnames):
            raise ValueError('samples.tsv requires columns: '+', '.join(COLUMNS))
        if 'coverage_type' in reader.fieldnames:
            raise ValueError('Legacy coverage_type is not supported (WGS only)')
        rows = list(reader)
    if not rows:
        raise ValueError('samples.tsv has no samples')
    units, samples = {}, {}
    libs = defaultdict(set)
    for line_no, row in enumerate(rows, 2):
        vals = {col: (row.get(col) or '').strip() for col in COLUMNS}
        for col in ('sample','species','library','lane'):
            if not ID.fullmatch(vals[col]):
                raise ValueError(f'Line {line_no}: invalid {col}: {vals[col]!r}')
        sp = vals['species']
        if sp not in refs:
            raise ValueError(f'Line {line_no}: species {sp} not in references')
        mode = vals['input_type']
        if mode not in ('fastq','bam','cram'):
            raise ValueError(f'Line {line_no}: invalid input_type {mode}')
        required = {'fastq': {'r1','r2'},'bam':{'bam'},'cram':{'cram'}}[mode]
        present = {x for x in ('r1','r2','bam','cram') if vals[x]}
        if present != required:
            raise ValueError(f'Line {line_no}: {mode} requires exactly {required}; got {present}')
        for field in required:
            vals[field] = abspath(vals[field], root)
            if check_files and (not Path(vals[field]).is_file() or Path(vals[field]).stat().st_size == 0):
                raise ValueError(f'Line {line_no}: missing or empty input: {vals[field]}')
        u = Unit(**vals)
        if u.key in units:
            raise ValueError(f'Line {line_no}: duplicated lane/library/sample {u.key}')
        if u.sample in samples:
            old = samples[u.sample]
            if (old.species,old.input_type) != (u.species,u.input_type):
                raise ValueError(f'Line {line_no}: inconsistent species/input_type for {u.sample}')
            if mode != 'fastq':
                raise ValueError(f'Line {line_no}: only one aligned input per sample is supported')
        samples[u.sample] = u
        units[u.key] = u
        if mode == 'fastq':
            libs[u.sample].add(u.library)
    for sample, libraries in libs.items():
        if len(libraries) != 1:
            raise ValueError(f'{sample}: multiple independent libraries require library-aware duplicate marking; currently unsupported')
    groups = defaultdict(list)
    for sample, u in samples.items():
        groups[u.species].append(sample)
    return list(units.values()),samples,{k:sorted(v) for k,v in groups.items()}


def reference_dict_path(fasta: Path) -> Path:
    name = fasta.name
    for ext in ('.fasta.gz','.fa.gz','.fna.gz','.fasta','.fna','.fa'):
        if name.endswith(ext):
            return fasta.with_name(name[:-len(ext)]+'.dict')
    return Path(str(fasta)+'.dict')


def check_references(refs: dict, root: Path, check_files=True):
    if not isinstance(refs,dict) or not refs:
        raise ValueError('References must be a nonempty mapping')
    for species, entry in refs.items():
        if not ID.fullmatch(species) or not isinstance(entry,dict):
            raise ValueError(f'Invalid reference configuration for {species}')
        contigs = entry.get('contigs')
        if not isinstance(contigs,list) or not contigs or any(not isinstance(c,str) or not ID.fullmatch(c) for c in contigs):
            raise ValueError(f'{species}: contigs must be a nonempty list of safe names')
        if len(contigs)!=len(set(contigs)):
            raise ValueError(f'{species}: duplicated contigs')
        ref = Path(abspath(str(entry.get('fasta','')),root))
        entry['fasta']=str(ref)
        d = reference_dict_path(ref)
        entry['dict']=str(d)
        if not check_files:
            continue
        req = [ref,Path(str(ref)+'.fai'),d]+[Path(str(ref)+'.'+ext) for ext in ('amb','ann','bwt','pac','sa')]
        for p in req:
            if not p.is_file() or p.stat().st_size==0:
                raise ValueError(f'Missing reference/index {p}; prepare reference indexes OUTSIDE the pipeline')
        fai = [line.rstrip('\n').split('\t') for line in Path(str(ref)+'.fai').read_text().splitlines() if line.strip()]
        seqs = [(a[0],int(a[1])) for a in fai]
        positions={name:i for i,(name,_) in enumerate(seqs)}
        if len(positions)!=len(seqs):
            raise ValueError(f'{ref}: duplicate FASTA contigs')
        if any(c not in positions for c in contigs):
            raise ValueError(f'{species}: configured contigs not found in {ref}.fai')
        indices=[positions[c] for c in contigs]
        if indices!=sorted(indices):
            raise ValueError(f'{species}: configured contig order must match FASTA index')
        dlines = [line for line in d.read_text().splitlines() if line.startswith('@SQ\t')]
        dict_seqs = {}
        for line in dlines:
            tags = dict(f.split(':',1) for f in line.split('\t')[1:] if ':' in f)
            dict_seqs[tags.get('SN','')]=int(tags.get('LN','0'))
        if dict_seqs != dict(seqs):
            raise ValueError(f'{species}: .dict and .fai disagree in reference names/lengths')
