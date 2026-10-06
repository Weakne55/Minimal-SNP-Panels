# GitHub Project setup

## Project name

**Minimal SNP Panels — Beluga & Polar Bear**

## Short description

Development and validation of minimal SNP panels for subpopulation assignment from WGS data in beluga and polar bear. Deadline: 25 Dec 2026.

## Project README

This project tracks development of a reproducible Bash/Python workflow for minimal SNP-panel design and validation.

The work has three parallel streams:

1. **Data & pipeline** — WGS processing, VCF/BCF generation, QC, population structure, SNP selection and validation.
2. **Literature & methods** — population structure, assignment, informative SNP selection and cross-validation.
3. **Documentation & reporting** — weekly measurable outputs, README, methods decisions and final scientific report.

### Main workflow

```text
WGS
→ VCF/BCF
→ QC
→ population structure
→ population_labels.tsv
→ population statistics / SNP ranking
→ candidate panels
→ cross-validation
→ final minimal SNP panel
```

### Definition of project completion

By 25 Dec 2026 the project should contain:

- reproducible Bash/Python pipeline;
- Dockerfile and README;
- population labels for both study datasets;
- final SNP panels for beluga and polar bear;
- classification and cross-validation results;
- comparison with available published panels/SNP sets;
- assessment of data sufficiency;
- final short-paper-style report.

Raw sequence data and large intermediate files remain on the HPC filesystem and are not tracked in GitHub.
