#!/usr/bin/env python3
"""Prepare a transparent manifest for a FUTURE population-structure study.

No VCF parsing, alignment, population inference, or cluster labels are produced.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


VCF_SUFFIXES = (".vcf", ".vcf.gz", ".bcf")


def _checked_variant(path_text: str) -> dict[str, object]:
    path = Path(path_text).expanduser().resolve()
    if not path.is_file():
        raise ValueError(f"VCF/BCF not found: {path}")
    if not any(path.name.lower().endswith(x) for x in VCF_SUFFIXES):
        raise ValueError(f"Unsupported variant filename: {path}")
    index_candidates = (
        [Path(str(path) + ".csi"), Path(str(path) + ".tbi")]
        if path.name.endswith(".vcf.gz") else
        [Path(str(path) + ".csi")]
        if path.name.endswith(".bcf") else []
    )
    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "index_present": any(index.is_file() for index in index_candidates),
    }


def make_plan(vcfs: list[str], species: str, assembly: str) -> dict[str, object]:
    """Validate only file existence/extension; document deferred scientific checks."""
    if not vcfs:
        raise ValueError("At least one --vcf is required")
    if not species.strip() or not assembly.strip():
        raise ValueError("Species and assembly must be non-empty")
    inputs = [_checked_variant(path) for path in vcfs]
    paths = [item["path"] for item in inputs]
    if len(paths) != len(set(paths)):
        raise ValueError("Duplicate VCF/BCF path")
    return {
        "status": "not_implemented",
        "purpose": "population_structure_research",
        "species_declared": species,
        "assembly_declared": assembly,
        "vcfs": inputs,
        "checks_completed": ["path_exists", "supported_extension", "duplicate_path"],
        "checks_pending": [
            "index_integrity", "bcftools_quickcheck", "assembly_and_contig_match",
            "sample_id_and_overlap", "record_order_and_variant_normalization",
            "GT_ploidy_quality_and_missingness", "variant_overlap_across_files",
        ],
        "planned_steps": [
            "harmonize_and_qc", "ld_prune", "pca", "ancestry_or_clustering",
            "interpret_and_validate_clusters", "review_and_write_population_labels_tsv",
        ],
        "population_labels_tsv": None,
        "warning": (
            "No genetic computation has been run. Declared species/assembly "
            "and VCF compatibility have NOT been checked. No labels inferred."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", action="append", required=True, help="Repeat for each VCF/BCF")
    parser.add_argument("--species", required=True, help="E.g. Delphinapterus_leucas")
    parser.add_argument("--assembly", required=True, help="Reference assembly accession")
    parser.add_argument("--out", type=Path, required=True, help="Output JSON plan file")
    args = parser.parse_args(argv)
    try:
        plan = make_plan(args.vcf, args.species, args.assembly)
        output = args.out.expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(plan, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    print(f"Wrote planning stub (NO analysis): {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
