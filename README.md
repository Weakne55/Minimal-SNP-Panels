# Minimal SNP Panels

Development and validation of minimal SNP panels for **subpopulation assignment** in:

- beluga (*Delphinapterus leucas*);
- polar bear (*Ursus maritimus*).

**Project deadline:** 25 December 2026.

## Goal

Build a reproducible Bash/Python pipeline that accepts:

```text
VCF/BCF + population_labels.tsv
```

and produces:

```text
QC
→ population statistics
→ SNP ranking
→ LD filtering
→ candidate panels
→ classification
→ cross-validation
→ minimal SNP panel
```

For the current study, WGS data are additionally processed to obtain the VCF/BCF files and to infer or validate genetic subpopulation labels.

## Current study

### Beluga
Publicly available WGS data, primarily from NCBI/SRA.

### Polar bear
Publicly available WGS data plus 12 in-house samples from the Kara group.

Raw and unpublished sequence data are **not stored in this repository**.

## Project architecture

### A. Study-specific data preparation

```text
WGS
→ mapping
→ BAM/CRAM
→ variant calling / joint genotyping
→ VCF/BCF
→ variant/sample QC
→ LD-pruned SNP set
→ population structure analysis
→ population_labels.tsv
```

### B. Reusable SNP-panel pipeline

```text
VCF/BCF + population_labels.tsv
→ QC
→ allele frequencies / population statistics
→ SNP ranking
→ LD filtering
→ candidate panels
→ classifier
→ cross-validation
→ minimal SNP panel
```

The reusable pipeline treats population labels as input and does not depend on how those labels were obtained.

## Repository layout

```text
.
├── .github/                GitHub issue templates and repository metadata
├── config/                 Example pipeline configuration
├── docs/                   Architecture, work plan and methodological decisions
├── literature/             Literature-methods table
├── metadata/               Small sample/label tables only
├── results/
│   └── final/              Small final tables suitable for version control
├── scripts/                Bash/Python analysis modules
├── slurm/                  SLURM job scripts
├── Dockerfile              Reproducible software environment
├── run_pipeline.sh         Main Bash entry point
└── README.md
```

## Data policy

Do **not** commit:

- FASTQ/FQ;
- BAM/CRAM and indexes;
- large VCF/BCF files;
- reference genomes;
- unpublished raw data;
- SLURM logs and temporary files.

Keep those on the HPC filesystem. Commit only code, configuration, metadata that may be shared, small summary tables, final marker lists and figures.

## Planned outputs

- standardized VCF/BCF for both species;
- population-structure results and `population_labels.tsv`;
- ranked SNP tables;
- candidate panels of different sizes;
- cross-validation metrics and confusion matrices;
- final minimal SNP panels for beluga and polar bear;
- comparison with published panels or SNP sets where reproducible;
- final report in the format of a short scientific article.

## Reproducibility

The main orchestration layer is Bash. Statistical and classification modules may be implemented in Python. Heavy computation is intended for SLURM. Dependencies are captured in the Dockerfile.

## Status

Active development. Methods are being fixed in parallel with literature review; see:

- [`docs/WORK_PLAN.md`](docs/WORK_PLAN.md)
- [`docs/METHODS_DECISIONS.md`](docs/METHODS_DECISIONS.md)
- [`literature/literature_methods.tsv`](literature/literature_methods.tsv)
