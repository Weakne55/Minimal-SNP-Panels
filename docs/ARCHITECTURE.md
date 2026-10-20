# Architecture

**Updated:** 20 October 2026. This document describes the target architecture; only the metadata catalogues are currently implemented in the repository.

## Scope and separation of concerns

Two distinct workflows must remain separate:

1. **Study-specific WGS workflow:** obtain or re-use variants from sequenced belugas and polar bears, analyze structure and establish labels.
2. **Reusable SNP-panel workflow:** take a VCF/BCF and a separate `population_labels.tsv` and select/validate a minimal panel. The reusable module must not depend on PCA/ancestry software used to obtain labels.

Each species has its own reference assembly, VCF/BCF, QC and population labels. Do **not** merge variants of both species into one VCF.

## A. Study-specific WGS workflow

```text
curated sample and run manifests
            |
            v
public SRA / approved restricted WGS
            |
            v
FASTQ / pre-existing BAM/CRAM or suitable VCF/BCF
            |
            v
verify reference assembly, read/sample identity and checksums
            |
            v
mapping -> BAM/CRAM + mapping/coverage/contamination QC
            |
            v
variant calling + joint genotyping (per species)
            |
            v
normalized VCF/BCF + sample- and variant-level QC
            |
            v
structure-specific SNP filtering and LD pruning
            |
            v
PCA + complementary ancestry / clustering analysis
            |
            v
check outliers, duplicate animals, relatives, batch effects,
cluster stability and concordance with published geography
            |
            v
species-specific population_labels.tsv + label provenance
```

The sequencing workflow is implemented for SLURM. For an existing VCF/BCF, skip FASTQ processing **only** after verifying build, variants, samples, provenance and consistency with the intended cohort.

## B. Reusable SNP-panel workflow

**Required inputs:** `VCF/BCF` and `population_labels.tsv` with matching sample IDs. Additional configuration specifies panel sizes, QC and validation settings. Grouping of BioSamples by independent animal is required where repeats/duplicates are known.

```text
VCF/BCF + population_labels.tsv
              |
              v
sample matching, genotype harmonization, dataset-level QC
              |
              v
outer validation split (animal groups; stratified where feasible)
       _______|________________________
      |                                |
      v                                v
TRAIN only                         HELD-OUT TEST
      |                                |
allele frequencies / FST               |  untouched during
SNP ranking / LD filtering             |  SNP selection, tuning
candidate panel sizes                  |  and fitting
inner tuning/validation                |
classifier fitting                     |
      |                                |
      +---------- trained model -------+
                       |
                       v
              held-out predictions
                       |
                       v
         per-class metrics and confusion matrices
                       |
                       v
       aggregate panel-size performance and SNP stability
                       |
                       v
      final panel refit + independent evaluation if available
```

The diagram does **not** mean that the final deployed marker set is different in every use. During cross-validation a fresh set is selected within each training fold to estimate the full *selection procedure*. After choosing a size/hyperparameters, a final panel is refitted on the designated development set, documented, and tested on truly independent data when available.

## Data currently available

| Species | Source of current metadata | Samples | Selected runs | Critical caveat |
|---|---|---:|---:|---|
| Polar bear | `data/metadata/polar_bear/curated/` | 144 | 286 | Two BioSamples correspond to one reported animal; independent-animal grouping and other duplicate/kinship checks required |
| Beluga | `data/metadata/beluga/curated/` | 106 | 106 | 8 neonates; variable estimated raw depth; `OTHER` assay type for 78 runs |

The 12 in-house Kara bears are not in this public cohort and need an access-controlled, provenance-checked integration step. The WGS catalogues do **not** contain FASTQ, mapped BAM, VCF or validated population labels.

## Sample/variant QC must be purpose-specific

- **Before structure inference:** review mapping coverage, sex/contamination, missingness, close relatedness and technical batches; then create a high-quality, mostly independent, LD-pruned autosomal SNP set.
- **For panel design:** guard against SNP dropout, missing genotypes, batch artifacts, excessive linkage and unstable signals; select candidate markers only from the training partition during performance evaluation.
- **Depth estimates:** NCBI/SRA `bases / genome_size` is pre-mapping raw-depth approximation, not mean usable autosomal depth; calculate coverage metrics from aligned reads before final cohort selection.
- **Population labels:** geographic sampling location and publication group annotations are evidence for comparison, not automatic genetic truth. Record ambiguous assignments as such.

Do not impose precise QC cutoffs until the reference assemblies, depth distributions and empirical QC have been reviewed. Fix chosen cutoffs/configuration and report excluded samples and variants once established.

## Validation and leakage risks

- All feature ranking, LD-aware panel selection, imputations fitted on the cohort and classifier tuning must use **training-only** data in each validation split.
- Use grouped splits to keep repeated BioSamples, sequencing libraries and known same-animal records together. Assess other relatives and sites as potential dependence structures.
- Balance evaluation across populations: report balanced accuracy, macro-F1, recall by group, confusion matrices and uncertainty. Record class counts per fold.
- Where SNP-based labels were inferred from the **same individuals' genomes**, CV is not by itself proof of external population assignment: the outcome labels are data-derived. Document this circularity and use external assignments / held-out cohorts for stronger validation if possible.
- Technical metadata (`bioproject`, sequencer, read depth, locality) are for QC and interpretation, **not genotype features** fed to the classifier.

## Planned implementation modules (not yet present)

When code is added, retain a small, transparent layout: `scripts/` for Bash/Python processing, `slurm/` for job launchers, `config/` for explicit settings, and a top-level entry point (`run_pipeline.sh`) for the reusable SNP-panel workflow. Environment reproducibility (Conda lock/environment specification and Dockerfile) should be documented; cluster operation need not depend on privileged Docker/Singularity execution.

Do not add empty directories or placeholder scripts merely to match this target layout. Build them together with tested runnable components.

## Required scientific outputs

For **each species**, ultimately deliver: reference/build and cohort manifest; QC summary; species-specific VCF/BCF; structure plots and `population_labels.tsv` with label provenance; ranked SNPs; candidate and final panel lists; held-out CV predictions and metrics; panel-size curve; comparison with a reproducible published panel or closest published SNP set; limitations and a short-paper-style report.
