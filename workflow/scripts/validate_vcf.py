#!/usr/bin/env python3
"""Validate a per-species cohort BCF and VCF; write compact QC evidence."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess


def run_capture(args):
    result = subprocess.run(args, check=True, capture_output=True, text=True)
    return result.stdout.strip()


def validate(bcf: Path, vcf: Path, samples: list[str], allowed: list[str], species: str) -> dict:
    bcf_ids = run_capture(["bcftools", "query", "-l", str(bcf)]).splitlines()
    vcf_ids = run_capture(["bcftools", "query", "-l", str(vcf)]).splitlines()
    expected = set(samples)
    if not expected or len(expected) != len(samples):
        raise ValueError("Empty or repeated sample IDs in FASTQ manifest")
    for name, observed in ((str(bcf), bcf_ids), (str(vcf), vcf_ids)):
        if len(observed) != len(samples) or set(observed) != expected:
            raise ValueError(f"Wrong cohort in {name}: {observed} != {samples}")
    if bcf_ids != vcf_ids:
        raise ValueError("BCF and VCF have a different sample order")
    counts = []
    for name in (bcf, vcf):
        total = int(run_capture(["bcftools", "index", "-n", str(name)]))
        if total <= 0:
            raise ValueError(f"No indexed variants in {name}")
        indexed_contigs = {}
        for line in run_capture(["bcftools", "index", "-s", str(name)]).splitlines():
            col = line.split("\t")
            if len(col) >= 3:
                indexed_contigs[col[0]] = int(col[2])
        if set(indexed_contigs) - set(allowed):
            raise ValueError(f"Unexpected non-autosomal contigs in {name}: {set(indexed_contigs)-set(allowed)}")
        counts.append({"format": "bcf" if name == bcf else "vcf.gz", "n_variants": total,
                       "variants_per_contig": indexed_contigs})
    if counts[0]["n_variants"] != counts[1]["n_variants"]:
        raise ValueError("Variant count changed during BCF -> VCF conversion")
    return {"species": species, "n_samples": len(samples), "sample_ids": bcf_ids,
            "assembly_contigs": allowed, "files": counts, "status": "validated"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bcf", type=Path, required=True)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--samples", nargs="+", required=True)
    parser.add_argument("--contigs", nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = validate(args.bcf, args.vcf, args.samples, args.contigs, args.species)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
