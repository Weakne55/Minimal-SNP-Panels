#!/usr/bin/env python3
"""Fast reference checks and sequence previews using only .fai and .dict.

Does not scan or modify the gigabyte-scale FASTA.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path


def dict_path(ref: Path) -> Path:
    return ref.with_suffix('.dict')


def read_fai(path: Path):
    results = []
    names = set()
    for number, line in enumerate(path.read_text(encoding='utf-8').splitlines(), 1):
        if not line.strip():
            continue
        fields = line.split('\t')
        if len(fields) < 5:
            raise ValueError(f'{path}:{number}: invalid .fai line')
        name = fields[0]
        length = int(fields[1])
        if name in names or length < 1:
            raise ValueError(f'{path}:{number}: duplicate name or invalid length: {name}')
        names.add(name)
        results.append((name, length))
    if not results:
        raise ValueError(f'No sequences found in {path}')
    return results


def read_dict(path: Path):
    results = []
    names = set()
    for number, line in enumerate(path.read_text(encoding='utf-8').splitlines(), 1):
        if not line.startswith('@SQ\t'):
            continue
        properties = dict(field.split(':', 1) for field in line.split('\t')[1:] if ':' in field)
        name = properties.get('SN')
        length = int(properties.get('LN', '0'))
        if not name or name in names or length < 1:
            raise ValueError(f'{path}:{number}: malformed @SQ record')
        names.add(name)
        results.append((name, length))
    if not results:
        raise ValueError(f'No @SQ records found in {path}')
    return results


def validate(ref: Path):
    if not ref.is_file() or ref.stat().st_size == 0:
        raise ValueError(f'Reference FASTA missing/empty: {ref}')
    indexes = [Path(str(ref) + '.fai'), dict_path(ref)]
    indexes += [Path(str(ref) + '.' + e) for e in ('amb', 'ann', 'bwt', 'pac', 'sa')]
    for p in indexes:
        if not p.is_file() or p.stat().st_size == 0:
            raise ValueError(f'Required index missing/empty: {p}')
    fai = read_fai(indexes[0])
    dictionary = read_dict(indexes[1])
    if fai != dictionary:
        raise ValueError('Reference .fai and .dict disagree in sequence names, lengths, or order')
    return fai


def inspect(ref: Path, limit: int):
    seqs = validate(ref)
    bp = sum(length for _, length in seqs)
    print(f'\nREFERENCE: {ref}')
    print(f'CONTIGS:   {len(seqs)}')
    print(f'TOTAL_BP:  {bp:,}')
    print(f'FAI:       {ref}.fai')
    print(f'DICT:      {dict_path(ref)}')
    print('INDEX\tNAME\tLENGTH_BP  (largest sequences; INDEX is original FASTA order)')
    largest = sorted(enumerate(seqs), key=lambda x: -x[1][1])[:limit]
    for ordinal, (name, length) in largest:
        print(f'{ordinal + 1}\t{name}\t{length}')
    # Suggested small smoke test: one long sequence. Never silently classify sex vs autosomes.
    _, (top_name, _) = largest[0]
    print(f'TEST-ONLY contigs: [{top_name}]')
    print('DO NOT use the test contig as a finished autosome list.')
    print('For full runs choose autosomes using assembly sequence report / biological metadata,')
    print('then list them in the original .fai order (not in sorted-by-size order).')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['check', 'inspect'])
    parser.add_argument('fasta')
    parser.add_argument('--limit', type=int, default=30)
    args = parser.parse_args()
    if args.limit < 1:
        parser.error('--limit must be positive')
    ref = Path(args.fasta).resolve()
    try:
        if args.command == 'check':
            seqs = validate(ref)
            print(f'OK: {ref} | {len(seqs)} sequences | FASTA/.fai/.dict/BWA files present; .fai/.dict match')
        else:
            inspect(ref, args.limit)
    except (ValueError, OSError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
