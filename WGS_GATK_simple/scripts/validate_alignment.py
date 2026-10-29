#!/usr/bin/env python3
"""Validate a pre-aligned BAM/CRAM header against samples.tsv and reference.

Only the SAM header is checked: actual read order and alignments are not read.
Input BAM/CRAM and reference files are never modified.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from model import reference_dict_path


def _tags(line: str) -> dict[str, str]:
    """Return SAM header KEY:VALUE tags (keep colons within values)."""
    return dict(field.split(":", 1) for field in line.split("\t")[1:] if ":" in field)


def _reference_headers(ref: Path):
    reference_lengths = {
        row.split("\t")[0]: int(row.split("\t")[1])
        for row in Path(str(ref) + ".fai").read_text().splitlines()
        if row.strip()
    }
    dictionary = {}
    for line in reference_dict_path(ref).read_text().splitlines():
        if line.startswith("@SQ\t"):
            tags = _tags(line)
            dictionary[tags["SN"]] = tags
    return reference_lengths, dictionary


def _sam_header(path, ref: Path):
    result = subprocess.run(
        ["samtools", "view", "-H", "-T", str(ref), path],
        check=True, capture_output=True, text=True,
    )
    order = None
    read_groups = 0
    sample_names = set()
    sequences = {}

    for line in result.stdout.splitlines():
        tags = _tags(line)
        if line.startswith("@HD\t"):
            order = tags.get("SO")
        if line.startswith("@RG\t"):
            read_groups += 1
            sample_names.add(tags.get("SM", ""))
        if line.startswith("@SQ\t"):
            name = tags.get("SN", "")
            if name in sequences:
                raise ValueError("Duplicate @SQ in input")
            sequences[name] = tags
    return order, read_groups, sample_names, sequences


def validate(path, sample, reference):
    """Check coordinate sort, sample RGs, reference contigs, lengths and MD5.

    MD5 is compared only when it is specified in both headers.
    The required contig *sets* must match; order is not compared.
    """
    ref = Path(reference)
    reference_lengths, dictionary = _reference_headers(ref)
    order, read_groups, sample_names, sequences = _sam_header(path, ref)

    if order != "coordinate":
        raise ValueError(f"{path}: expected @HD SO:coordinate, got {order}")
    if read_groups == 0 or sample_names != {sample}:
        raise ValueError(f"{path}: @RG SM must equal {sample}, got {sample_names}")
    if set(sequences) != set(reference_lengths):
        raise ValueError(f"{path}: @SQ contig set differs from selected reference")

    for name, tags in sequences.items():
        if int(tags.get("LN", "-1")) != reference_lengths[name]:
            raise ValueError(f"{path}: reference length mismatch for {name}")
        ref_md5, bam_md5 = dictionary[name].get("M5"), tags.get("M5")
        if ref_md5 and bam_md5 and ref_md5 != bam_md5:
            raise ValueError(f"{path}: reference MD5 differs for {name}")
    return True


if __name__ == "__main__":
    try:
        validate(*sys.argv[1:4])
    except Exception as error:
        sys.exit(f"ERROR: {error}")
