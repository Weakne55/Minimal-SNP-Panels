# Контракт входных и выходных данных — VCF/BCF и метки популяций

## WGS → downstream

One species + one reference assembly per run. Typical WGS output:
`cohort.raw.vcf.gz` / `.bcf` with an index, genotype/sample IDs and a
recorded assembly accession (e.g., `GCA_059928025.1` for beluga or
`GCA_062004589.1` for polar bear). "Raw" does not imply SNP-panel QC.

WGS need not generate labels, SNP panels or population structure plots.
Downstream programs must accept *external* compatible VCF/BCF as well.

## VCF/BCF collection

An ordered list of input files is provided. Before SNP selection or population
inference, check reference accession and contigs, samples and duplicates,
variant `CHROM:POS:REF:ALT`, normalization (left alignment and ALT splitting),
ordering, genotype ploidy, missingness, indexing and representativeness.
Input files can contain (1) disjoint chromosomes but identical samples,
or (2) overlapping sites with distinct samples. These are different merge
problems (`bcftools concat` vs `bcftools merge`) and cannot be guessed from
filenames. Overlapping samples or sites require an explicit policy; reject
ambiguous mixing instead of silently merging.

## Metadata

Population labels: UTF-8, tab-delimited, exact column names:

```tsv
sample	population
sample_001	population_A
sample_002	population_B
```

`sample` is unique, non-empty and matches the VCF sample ID exactly.
Population labels are supplied to the *panel-selection* module, not WGS.
Population-structure inference (PCA, ADMIXTURE, clustering) is independent,
optional research; its reviewed final labels can later become this TSV.

## Minimal-panel selection

The user specifies candidate sizes and a threshold/metric, e.g. balanced
accuracy 0.95. "Minimal" is the smallest **tested** panel for which the
chosen validation rule supports the threshold, not necessarily a global
mathematical optimum. No qualifying size → explicit status `not_achieved`.
Nested CV (or strictly independent equivalent) is needed for tuning panel
size without optimistic bias. If samples are related, splits must account
for family/relatedness/biological grouping; class balance and uncertainty
must be reported. No random sampling of the same individual across folds.

## Intended outputs (not implemented)

`qc/`, `snp_statistics.tsv`, `ranked_snps.tsv`, `candidate_panels/`,
`validation/`, `panel_performance.tsv`, `final_panel.tsv` (only if qualified),
`run_metadata.yaml` with seed, version, input hashes and config.
