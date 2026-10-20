# Minimal SNP Panels — Beluga & Polar Bear

Development and validation of minimal SNP panels for **genetic subpopulation assignment** in:

- beluga (*Delphinapterus leucas*);
- polar bear (*Ursus maritimus*).

**Delivery deadline:** 25 December 2026.  
**Status of this README:** 20 October 2026.

## Objective and scope

The reusable Bash/Python workflow will take **VCF/BCF + `population_labels.tsv`** and produce quality-controlled informative SNP rankings, LD-filtered candidate panels, assignment results, cross-validation metrics and a minimal validated panel for each species.

The current study additionally requires obtaining suitable VCF/BCF from public WGS, checking genetic population structure and assigning/validating population labels. **Population-label inference is a separate study-specific analysis**, not a required operation of the reusable panel-design workflow.

## Current data inventory (20 October 2026)

| Species | Public WGS BioSamples | Selected SRA runs | Metadata status |
|---|---:|---:|---|
| Polar bear | 144 | 286 | Sample and run registries curated; independence and genome-level QC remain pending |
| Beluga | 106 | 106 | Sample and run registries curated; mapped depth and genotype-level QC remain pending |

- **Polar bear:** a confirmed cross-BioSample identity links two records for one animal; additional aliases/relatedness require review. Published population and locality fields are **not** finalized genetic labels.
- **Beluga:** 78 BioSamples are from `PRJNA1182041` and 28 from `PRJNA1328577`. There are 98 adults and 8 neonates. Estimated raw coverage is at least 15× for 98 samples, including 90 adults, and at least 20× for 60 samples, including 52 adults. These estimates are calculated from SRA base counts, **not** mapped coverage. In `PRJNA1182041`, the SRA assay field is `OTHER`; source publication/project information supports their inclusion as WGS candidates but should be reviewed before processing.
- **Additional in-house material:** 12 Kara polar bears are within the study scope but **are not included** in the 144 public-BioSample registry or stored here. Add their restricted metadata and genotypes to the analysis only through an approved access route; do not publish them without permission.
- The Montana et al. (2024) **low-coverage** beluga cohort is not part of this WGS registry.

All numbers above describe metadata records, **not** independently validated individuals or analysis-ready genomes. Reference assemblies, final inclusion thresholds and genetic labels are not fixed yet.

## Repository structure — implemented now

```text
.
├── README.md
├── LICENSE_NOTICE.md
├── PROJECT_SETUP.md
├── docs/
│   ├── ARCHITECTURE.md
│   ├── GITHUB_SETUP.md
│   ├── METHODS_DECISIONS.md
│   └── WORK_PLAN.md
├── literature/
│   └── literature_methods.tsv
└── data/
    └── metadata/
        ├── beluga/
        │   ├── README.md
        │   ├── metadata_column_classification.tsv
        │   ├── raw/
        │   └── curated/
        └── polar_bear/
            ├── README.md
            ├── metadata_column_classification.tsv
            ├── raw/
            └── curated/
```

This tree lists the **current tracked project content**, rather than directories planned for later implementation. The preliminary, untracked `sra_inventory_2026-10-13_v2/` directory is superseded by `data/metadata/` and can be removed locally. The GitHub project configuration is documented in `docs/GITHUB_SETUP.md`.

## Metadata: authoritative inputs

| Purpose | Polar bear | Beluga |
|---|---|---|
| Selected samples | `data/metadata/polar_bear/curated/polar_bear_samples_working_144.tsv` | `data/metadata/beluga/curated/beluga_samples_working_106.tsv` |
| Selected sequencing runs | `data/metadata/polar_bear/curated/polar_bear_runs_working_286.tsv` | `data/metadata/beluga/curated/beluga_runs_working_106.tsv` |
| Full metadata and provenance | `data/metadata/polar_bear/raw/` | `data/metadata/beluga/raw/` |
| Field definitions and caveats | `data/metadata/polar_bear/README.md` | `data/metadata/beluga/README.md` |

One curated **sample** record is keyed by `biosample`; one curated **run** record is keyed by `run_accession`. A BioSample is not automatically an independent animal. Use only the curated run manifests to acquire the selected public WGS runs; full `raw/` SRA tables are for traceability and may include non-cohort runs.

The sample tables contain source-locality or published group annotations, but **neither species currently has a finalized `population_labels.tsv`**. Do not feed `geographic_candidate_group`, `published_group` or `geographic_group_from_metadata` directly into the classifier as if they were verified genetic populations.

## Planned workflows

### A. Study-specific WGS and population-structure workflow

```text
Public WGS runs / approved local WGS
    -> FASTQ acquisition and integrity checks
    -> species-specific reference mapping
    -> sorted and indexed BAM/CRAM + QC
    -> variant calling / joint genotyping
    -> species-specific VCF/BCF + sample/variant QC
    -> LD-pruned SNPs for PCA and structure analysis
    -> checked population_labels.tsv
```

### B. Reusable SNP-panel workflow

```text
Species-specific VCF/BCF + checked population_labels.tsv
    -> QC and genotype harmonization
    -> train/test split by independent animal
    -> within-training SNP statistics, ranking and LD filtering
    -> candidate panels and classifier fitting
    -> held-out evaluation, with tuning confined to training data
    -> panel-size/quality curves, per-population metrics
    -> minimal validated SNP panel
```

Selection of markers, LD pruning for the predictive panel and model tuning must be **refitted inside validation folds**. For population labels inferred from the same WGS data, performance must be described as reconstruction of dataset-defined groups unless independently established labels/validation are available.

## Implementation status and next milestone

**Completed:** two-species public WGS metadata inventories, SRA run manifests and metadata documentation.

**Not yet confirmed as completed:** mapping prototype, reference assembly decision, joint-calling pipeline, VCF/BCF, SNP QC, population-label inference, informative SNP selection and classification.

**Next technical milestone:** run a reproducible FASTQ -> BAM/CRAM -> initial VCF test for at least one sample per species, then extend to joint genotyping and full-cohort SLURM processing. The formal task specification targets **8 November** for the WGS-to-VCF pipeline and **15 November** for population structure/labels; see `docs/WORK_PLAN.md`.

Planned implementation: Bash orchestration and SLURM scripts, Python for statistics and classifiers; reproducible environment specifications and a Dockerfile are deliverables **to be created**, not files already present. The pipeline may use pre-existing suitable VCF/BCF when their provenance, reference and sample metadata can be verified.

## Literature and reporting

`literature/literature_methods.tsv` is the **single working evidence table** for population-structure methods, informative SNP selection, published marker sets and validation/benchmarking. The current file contains only its header; a literature review has **not** yet been entered. See `docs/METHODS_DECISIONS.md` for criteria and `docs/WORK_PLAN.md` for milestones.

Research results should be produced every week while literature and final reporting proceed in parallel. The contractual task numbering is detailed in `docs/WORK_PLAN.md`.

## Data handling and licensing

Do **not** commit raw FASTQ, BAM/CRAM, large VCF/BCF, reference genomes, temporary analysis products, credentials, or unpublished in-house data. Put large files on the authorized HPC filesystem. Only share metadata whose redistribution is permitted; NCBI availability does not itself establish licensing for every associated annotation. See `LICENSE_NOTICE.md` before public redistribution or selecting a code license.
