# Work plan

Work runs in parallel across:
1. data/pipeline;
2. literature/methods;
3. documentation/reporting.

| Period | Data / pipeline | Literature / methods | Measurable output |
|---|---|---|---|
| 6–11 Oct | Repository structure, WGS inventory, references | Key beluga/polar-bear population-genomics papers | sample inventory, repository, literature table |
| 12–18 Oct | Test WGS download, FASTQ → BAM/CRAM, SLURM, Docker | mapping and variant-calling approaches | working mapping prototype |
| 19–31 Oct | variant calling / joint genotyping | joint calling and variant QC | joint VCF/BCF for each species |
| 1–7 Nov | sample/SNP QC, LD-pruned data | MAF, missingness, LD and structure analysis | filtered VCF/BCF + QC summary |
| 8–15 Nov | PCA + population structure, fix labels | PCA, ancestry/clustering, number of clusters | plots + `population_labels.tsv` |
| 16–22 Nov | allele frequencies, FST, SNP ranking | informative-marker selection | SNP metrics + ranked SNPs |
| 23–29 Nov | LD filtering, candidate panel sizes, first classifier | population assignment methods | candidate panels + first classification |
| 30 Nov–6 Dec | correct CV with selection inside training folds | validation, leakage, imbalance, relatedness | CV metrics + confusion matrices |
| 7–13 Dec | repeated CV, SNP stability, preliminary final panels | published panels / SNP sets | preliminary final panels |
| 14–18 Dec | benchmark published sets, data sufficiency | limitations and interpretation | comparison table + limitations |
| 19–22 Dec | finalize pipeline, Docker, README, clean run | targeted reference checks | reproducible release candidate |
| 23–25 Dec | fixes and packaging only | no new major methods | final delivery |

## Critical checkpoints

- **31 Oct:** WGS → VCF/BCF works.
- **15 Nov:** QC datasets and population labels are fixed.
- **6 Dec:** SNP-selection + cross-validation pipeline works.
- **13 Dec:** preliminary final panels exist.
- **23–25 Dec:** no scope expansion.
