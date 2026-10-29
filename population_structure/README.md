# Population structure — research stub

**Not implemented.** This is a separate research step for inferring and validating
population assignments from genome-wide SNPs; **it is not part of the WGS DAG or
the main SNP-panel DAG**.

The current `prepare.py` checks only the existence/extensions of one or more
VCF/BCF files and writes a JSON planning manifest. It does **not** validate
genotype contents, reference consistency, sample names, or variant overlaps.
It does not run PCA, LD pruning, ADMIXTURE or clustering. In particular,
it does **not** create arbitrary `population_labels.tsv`.

```bash
python population_structure/prepare.py \
  --vcf /path/to/cohort.raw.vcf.gz \
  --species Delphinapterus_leucas \
  --assembly GCA_059928025.1 \
  --out results/population_structure/beluga/plan.json
```

Later the analysis must validate input/reference, harmonize multiple VCFs,
perform sample QC and LD pruning, run PCA plus an independent structure method,
inspect ambiguous/admixed samples and only then produce reviewed labels:

```tsv
sample	population
sample_001	pop_A
sample_002	pop_B
```

This file is an **illustration**, not a generated label table. If labels were
obtained by clustering on the same WGS used to select SNP, classification
accuracy measures the ability to reproduce these cluster labels, not independent
biological confirmation of populations.
