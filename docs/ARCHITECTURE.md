# Architecture

## Study-specific workflow

```text
NCBI / local WGS
        |
        v
FASTQ / BAM / CRAM
        |
        v
mapping + post-processing
        |
        v
variant calling / joint genotyping
        |
        v
VCF/BCF
        |
        v
sample + variant QC
        |
        v
LD-pruned SNP set
        |
        v
population structure analysis
        |
        v
population_labels.tsv
```

## Reusable panel-design workflow

```text
VCF/BCF + population_labels.tsv
        |
        v
QC / harmonization
        |
        v
allele frequencies + population statistics
        |
        v
SNP ranking
        |
        v
LD filtering
        |
        v
candidate panels (N SNP)
        |
        v
classification
        |
        v
cross-validation
        |
        v
minimal validated panel
```

## Methodological constraint

For unbiased performance estimation, SNP ranking/selection must be performed **inside the training data of each validation split**. A panel selected once on the full dataset and then evaluated on the same samples is not a valid estimate of out-of-sample assignment performance.
