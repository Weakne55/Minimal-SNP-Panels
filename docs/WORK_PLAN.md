# Work plan and milestone status

**Updated:** 20 October 2026.  
**Contractual research period:** 1 October – 25 December 2026.  
**Delivery deadline:** **25 December 2026**.

Work proceeds in parallel across (a) data/computation/code (approximately 60–70%), (b) task-driven literature/methods (approximately 20–30%), and (c) documentation/reporting (approximately 10%). Each week should produce a tangible item: a verified manifest, tested script, BAM/CRAM, VCF/BCF, QC summary, figure, ranked SNPs or validation metrics.

## Formal specification: task numbers and dates

The authoritative task numbering below follows the supplied `ТЗ Карпович.docx`. These tasks **overlap**; they are not consecutive prerequisites.

| Task | Scope | Formal period | Required result |
|---|---|---|---|
| **1.1** | WGS cohort formation and technical metadata | 6–25 Oct | Beluga/polar-bear sample manifests, sources, FASTQ or available genome files |
| **1.2** | WGS -> VCF/BCF, Bash/Python + SLURM, Dockerfile | 12 Oct–8 Nov | Standardized VCF/BCF per species; runnable documented workflow |
| **1.3** | Population structure on WGS SNPs | 26 Oct–15 Nov | LD-pruned SNPs, PCA/ancestry analyses, `population_labels.tsv` per species |
| **1.4** | SNP informativeness, panels, cross-validation | 9 Nov–13 Dec | Ranked SNPs, candidate/final panels, held-out metrics |
| **1.5** | Published SNP panel / marker-set comparison | 30 Nov–20 Dec | Reproduced comparison where feasible; limitations otherwise |
| **1.6** | Data sufficiency, release and report | 14–25 Dec | Reproducible pipeline, final evaluation and short-paper-style report |

**Interpretation:** literature for task 1.3 covers PCA/ancestry and label construction; the review/comparison of published SNP panels is formally **task 1.5**, with method selection also supporting task 1.4.

## Verified progress as of 20 October

| Item | Status | Evidence / gap |
|---|---|---|
| GitHub repository and project documentation | **Available** | README, task docs, metadata directories committed |
| Public polar bear WGS registry | **Available, preliminary QC** | 144 BioSamples / 286 selected SRA runs; known cross-BioSample duplicate |
| Public beluga WGS registry | **Available, preliminary QC** | 106 BioSamples / 106 SRA runs; 98 adults and 8 neonates |
| 12 restricted Kara polar bears | **Not yet integrated** | Access-controlled cohort and metadata still required |
| `literature_methods.tsv` | **Header only** | No literature records entered yet |
| Reference assembly choices | **Not recorded as finalized** | Must be decided and documented |
| WGS download / mapping test | **Not demonstrated in this repository** | No reproducible FASTQ -> BAM/CRAM test artifact shown |
| Joint variant calling and species VCF/BCF | **Not demonstrated** | Code and test outputs not yet present |
| Population structure and labels | **Not done in repository** | Requires genotype data and QC |
| Ranking, panels, classification and CV | **Not done** | Depends on VCF/BCF + validated labels |

`data/metadata/` is the authoritative **metadata registry**, not a substitute for acquired FASTQ files, QC or genetic population labels. WGS count is a count of BioSamples, not verified independent animals.

## Immediate execution plan and measurable weekly deliverables

| Period | Computation / data (highest priority) | Literature / methods in parallel | Tangible checkable result |
|---|---|---|---|
| **20–25 Oct** | Clean old untracked inventory; verify `curated/` manifests; decide reference assemblies; acquire at least one tractable WGS per species; benchmark FASTQ -> aligned BAM/CRAM | WGS preprocessing, references, variant calling; start literature table | Updated docs, reference choices and successful small mapping test or explicit download blocker |
| **26 Oct–1 Nov** | Implement reproducible read QC, mapping and small-sample variant calling; begin batch-processing plan | Joint genotyping/QC, filtering parameters | Test BAM/CRAM QC and first callable species-specific variants |
| **2–8 Nov** | Extend to joint genotyping; submit and monitor selected-cohort jobs; produce VCF/BCF or document a reproducible tested production run with unfinished computation | Mapping/calling comparison and quality thresholds | **Formal 1.2 deadline:** runnable workflow, documented VCF/BCF status, reproducible logs/QC |
| **9–15 Nov** | Normalize/filter variants, LD prune for structure, PCA plus ancestry/clustering; review duplicates and batches; establish labels | PCA, ancestry, cluster selection and uncertainty | **Formal 1.3 deadline:** structure plots and supported `population_labels.tsv` or explicitly reported blockers |
| **16–22 Nov** | Population allele frequencies/FST; pilot SNP ranking within training partitions | Informative SNP selection, LD effects | Per-population SNP statistics and preliminary ranked table |
| **23–29 Nov** | Build and benchmark several LD-filtered panel sizes and a baseline classifier | Population assignment; performance metrics | Candidate SNP lists and first grouped CV estimates |
| **30 Nov–6 Dec** | Implement nested/grouped validation, held-out predictions, class-specific metrics | Leakage, relatedness and imbalanced classes | CV metrics, confusion matrices, panel-size curves |
| **7–13 Dec** | Repeat CV, assess marker stability and select preliminary minimal panel | Published panel evaluation criteria | **Formal 1.4 deadline:** documented preliminary final panels |
| **14–20 Dec** | Published-panel/SNP-set benchmark where technically possible, data sufficiency and robustness | Critical comparison and limitations | **Formal 1.5 deadline:** comparison summary; draft report |
| **21–25 Dec** | Clean release, test end-to-end entry point, verify configs and package results | Final report revision only | **Formal 1.6 deadline:** reproducible delivery; no major new methods |

The week-by-week entries are **operational checkpoints**, not guarantees of computation finishing by those dates. Work requiring network access/HPC capacity should be initiated early; record blockers rather than reporting an unrun pipeline as complete. Existing suitably documented public VCF/BCF may be considered to mitigate WGS calling delays.

## Current critical path

```text
Choose reference assemblies / sample subset
    -> acquire FASTQ or validated existing VCF/BCF
    -> mapping, sample QC, variant calling, joint genotyping
    -> per-species QC and LD-pruned SNPs
    -> PCA/ancestry and checked population labels
    -> leakage-aware SNP selection + assignment validation
    -> published comparison, sufficiency, release/report
```

Work on the reusable panel-design code can start with a synthetic or small example VCF **before** complete WGS processing. The software interface remains VCF/BCF + labels, independent of label-inference implementation. This mitigates delays without replacing the scientific WGS analysis.

## Documentation and reporting convention

For each weekly deliverable, record the Git commit, sample counts and input manifest version, reference/tool version, launch command, QC/results paths, unresolved limitations and the corresponding report figure/table. Keep large results on approved HPC storage; commit only shareable metadata, small tables, scripts and concise summaries.

The literature table is a continuous evidence log, **not** a deferred end-of-project task. First entries should cover beluga and polar-bear population structure, WGS/calling choices, informative SNP selection and published marker benchmarks, with DOI/URL and an explicit `what_we_use` field.
