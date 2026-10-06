# Methods decisions

This file records decisions that should not be silently changed.

## Fixed architecture

- Main reusable input: `VCF/BCF + population_labels.tsv`.
- Population labels are external input to the panel-design pipeline.
- For the current study, labels are inferred/validated separately from WGS population structure.

## Current baseline

- Bash for orchestration.
- Python for statistics/classification where appropriate.
- SLURM for heavy computation.
- Dockerfile for reproducible dependencies.

## Population structure

- Do not equate geographic sampling group with genetic population automatically.
- Use PCA plus at least one justified clustering/ancestry approach.
- The method for choosing the number of clusters must be documented.

## Panel selection

- FST is a baseline population-differentiation metric, not necessarily the only ranking criterion.
- Remove redundant markers with LD-aware filtering.
- Preliminary panel sizes: 10, 20, 50, 100, 200 SNP.

## Validation

- SNP selection occurs inside training folds.
- Report class-level metrics, not only global accuracy.
- Track confusion matrices.
- Examine imbalance and relatedness where relevant.
- Repeated validation should be used to assess SNP stability.

## Literature benchmark

- Do not assume that 2–3 fixed published panels exist for each species.
- Reproduce fixed panels when technically possible.
- Otherwise compare against the closest published SNP sets / population-assignment approaches.
