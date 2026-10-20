# Methods decisions and open questions

**Updated:** 20 October 2026. Keep *agreed rules* separate from *choices still to be made*. This file records methodological decisions, not unverified claims that analyses have already run.

## 1. Fixed architecture — agreed

- Reusable entry point: **species-specific `VCF/BCF + population_labels.tsv`**.
- Population labels are input to the reusable panel pipeline. In this study they will be inferred/checked separately from the WGS cohort.
- Bash for workflow orchestration; Python for statistics/classification; SLURM for HPC processing. A Dockerfile is a formal deliverable, but a working cluster path must not depend on Docker privileges.
- Process beluga and polar bear against their own validated reference assemblies. No cross-species VCF mixing.

## 2. Cohorts and provenance — verified metadata, not finalized analysis sets

| Species | Public BioSamples | Selected runs | Current metadata issues |
|---|---:|---:|---|
| Polar bear | 144 | 286 | At least one confirmed same-animal pair under separate BioSamples; possible additional aliases; coverage and ancestry heterogeneity |
| Beluga | 106 | 106 | 78 `PRJNA1182041` runs with SRA assay `OTHER`; 8 neonates; raw-depth differences |

- Use `data/metadata/<species>/curated/*samples_working*.tsv` and `*runs_working*.tsv` for working inventories; `raw/` holds fuller provenance and audit data.
- BioSample, sequencing run and **animal** are distinct entities. Deduplicate or group known re-sequenced animals before genotype-based analyses and validation.
- SRA raw-depth estimate is *not* mapped/usable coverage. Define definitive thresholds only after alignment and QC.
- Beluga preliminary cohorts: 98 adults; 90 adults with estimated raw depth ≥15×; 52 adults with estimated raw depth ≥20×. These are useful sensitivity cohorts, **not yet final QC decisions**. Neonates remain documented and can be analyzed separately if justified.
- Twelve additional Kara polar bears are in the approved research scope but outside the public registry and must be integrated under the correct access and metadata policy; their presence/availability must not be assumed.
- Low-coverage beluga samples from Montana et al. (2024) are **not** part of the current WGS acquisition cohort.

## 3. WGS preprocessing and calling — to decide/test

- **Pending:** reference assembly/accession for each species, read trimming policy, aligner, duplicate marking policy, variant caller and joint-genotyping implementation.
- **Pending:** minimum mapped/autosomal depth, covered fraction, sample/variant missingness, quality-score thresholds, contamination/sex review and population-specific filtering exceptions.
- **Required:** record reference accession/version, coordinate normalization and tool versions in each workflow run.
- **Permitted:** use an independently sourced suitable VCF/BCF if sample identities, reference build, variant representation, quality and ascertainment are sufficiently documented.
- First demonstrate a FASTQ -> BAM/CRAM -> VCF test on a small sample subset before scaling to all WGS.

## 4. Population labels — agreed evaluation principles

- Geography and `published_group` annotations are **candidate evidence**, not proven genetic labels.
- Infer structure using PCA plus at least one complementary ancestry/clustering method on a QC-filtered, LD-pruned SNP set. Justify cluster-number choice and assess stability across SNP and sample subsets.
- Review kinship, duplicates, outliers, uneven sample sizes, depth and sequencing-project confounding before assigning labels.
- Maintain label provenance (method, confidence, ambiguity and any publication evidence). Samples with unstable/ambiguous ancestry must not silently become high-confidence training labels.
- Evaluate known caveats: two Southeast Greenland polar bears (`D24051`, `D24087`) have published genetic affinities atypical of their sampling region; do not automatically label these as typical SEG for training.

## 5. Informative SNP selection — baseline agreed, methods pending

- Baseline ranking candidates: population allele-frequency contrasts and **FST**; consider complementary population-informativeness measures/feature selection where justified by sample size and literature.
- LD-aware selection is required to reduce redundancy. Define chromosome/position constraints and window thresholds from the data and reference build.
- Start with panel sizes **10, 20, 50, 100, 200 SNP**; refine the grid if needed, rather than declaring any one count the minimum in advance.
- Candidate panel SNPs must be genotype-callable and portable across the intended sample set. Technical covariates or paper/geography labels must never be supplied as SNP predictors.

## 6. Validation — methodological requirements

- Evaluate the **entire selection procedure**: training-fold SNP QC fitting where applicable, frequency/FST estimates, ranking, LD filtering, panel selection, classifier fitting and parameter tuning.
- No information from held-out genotypes/labels may be used to select markers, tune hyperparameters or choose panel size for that fold. Use an **inner** split/validation for tuning when necessary and an **outer** held-out split for performance estimation.
- Split by independent animal, not run/BioSample; keep confirmed duplicates in one fold or retain one representation. Consider stronger family/locality/project hold-outs as sensitivity analyses where population sample sizes permit.
- Report per-population recall/precision, macro-F1, balanced accuracy, confusion matrix, number of evaluable individuals and variability across splits. Do not rely only on pooled accuracy.
- Repeated CV will measure prediction and marker selection stability; define feasibility before running on very small populations.
- If target labels are derived from clustering the same WGS individuals, internally cross-validated agreement tests recovery of these **dataset-derived** groups; it is not an independent accuracy estimate for objectively known population origin. Make this limitation explicit and seek external validation where feasible.

## 7. Literature and published-panel benchmark

The single table `literature/literature_methods.tsv` supports all relevant methods sections:

- **Formal task 1.3:** population-structure methods (PCA, ancestry/clustering, label assignment).
- **Formal task 1.4:** population-differentiation statistics, SNP ranking, LD-aware panel design and validation.
- **Formal task 1.5:** published fixed SNP panels or closest available SNP sets, comparison feasibility and limitations.

For each paper, record study species, sample sizes/populations, assembly, methods, selection, validation, panel size, main result and **what is reusable for our analysis**. Distinguish *a fixed reproducible SNP panel* from *a genome-wide marker set for population structure*. Do not assert that ready-to-test minimal panels exist without marker IDs, alleles, genomic coordinates/reference-build information and a usable comparison protocol.

**Current status:** `literature_methods.tsv` has only a header. No structured literature screening or paper-level method entries have yet been recorded in that file.

## 8. Open methodological decisions (resolve by measured tests)

| Decision | Evidence/experiment needed | When |
|---|---|---|
| Reference genome per species | assembly assessment, download availability, coordinate convention | Before mapping |
| WGS calling route | small FASTQ -> BAM/CRAM -> VCF benchmark, resource use | By 8 Nov |
| Final cohort/QC cutoffs | aligned coverage, mapping QC, relatedness, missingness | Before structure analysis |
| Cluster-number and labels | PCA/ancestry stability, geography concordance, problematic samples | By 15 Nov |
| Selection metric / LD settings | training-only pilot tests with fold-safe evaluation | Nov |
| Classifier / panel-size choice | balanced performance, uncertainty and stability | Nov–Dec |
| Literature comparators | accession/marker coordinates and method reproducibility | By 20 Dec |

## 9. Non-negotiable deliverables

By **25 December 2026**: documented reproducible Bash/Python + SLURM pipeline; Dockerfile and configs; both species' QC-ed VCF/BCF and population labels; candidate/final SNP lists; leakage-aware validation; published-method comparison or justified feasibility limit; data-sufficiency assessment; short-paper-style report. These are targets, **not claims of current completion**.
