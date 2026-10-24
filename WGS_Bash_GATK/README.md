# WGS → joint VCF: BWA-MEM + GATK HaplotypeCaller / GenomicsDB + Bash/Python/SLURM

**Version:** 2026-10-21, initial server-test candidate. **Status:** Python/SLURM DAG unit-tested; **NOT** yet end-to-end validated with real sequencing software/data.

Purpose: produce one unfiltered cohort variant VCF per species from short-read WGS for later, separate population structure / SNP panel projects. Supports independent species (examples: beluga `Delphinapterus_leucas` and polar bear `Ursus_maritimus`). Nothing in this pipeline assigns populations or selects SNPs.

**Scope agreed:** WGS only; BWA-MEM **not** BWA-MEM2; GATK HaplotypeCaller in gVCF mode and GenomicsDBImport+GenotypeGVCFs; no BQSR, QC gates, SNP filtering, imputation or SLURM resource profiling. BAM QC metrics are produced for inspection, but processing is not automatically stopped by poor QC.

## Data flow

```text
samples.tsv + config.yaml + EXISTING references / indexes
 |
 |-- FASTQ (per lane)  -> [optional fastp] -> BWA-MEM -> samtools collate/fixmate/sort
 |                  -> merge sorted lanes -> samtools markdup -> indexed BAM
 |
 |-- prealigned indexed BAM -> validate header/index -> use directly (no copy)
 |-- indexed CRAM -> validate -> convert to indexed BAM (no repeat markdup)
 |
 +-> samtools flagstat/stats + mosdepth -> MultiQC (report only)
 |
 +-> GATK HaplotypeCaller -ERC GVCF, per (sample, contig)
        -> gvcf/<species>/<sample>/<contig>.g.vcf.gz
     -> per-contig GenomicsDBImport (all samples of SAME species)
     -> per-contig GenotypeGVCFs
     -> bcftools concat (in FASTA contig order)
     -> vcf/<species>/cohort.raw.vcf.gz[.csi]
```

GATK parameters are intentionally simple (diploid by default). **Use only autosomes in the initial `contigs` configuration.** Sex chromosomes/mitochondria require explicit sex- and locus-dependent ploidy policy and will be handled separately.

## 1. Extract to shared filesystem

Upload ZIP to cluster, then on a login node:

```bash
unzip WGS_Bash_GATK_SLURM.zip
cd WGS_Bash_GATK
```

The project is self-contained. It has **no Snakemake dependency**. SLURM nodes must see the same paths to scripts, references, inputs, and workdir (shared filesystem). `run_pipeline.sh` resolves the project directory automatically.

## 2. Create software environment once

```bash
conda env create -f envs/wgs.yaml
conda activate wgs-gatk
bash scripts/check_stack.sh
python -m unittest discover -s tests -v
```

An equivalent `mamba env create -f envs/wgs.yaml` is fine. Package versions are specified in `envs/wgs.yaml` (BWA 0.7.18, GATK 4.6.2.0, samtools 1.21, etc.). Create an explicit lock/export **after** resolution on the cluster (e.g., `conda list --explicit > versions-lock.txt`). The supplied Dockerfile is optional and is **not** automatically run by SLURM; if the cluster permits Apptainer, integration is a separate step.

## 3. References are only EXTERNAL inputs

**No download, copying, indexing, symlinking, or modifying references is part of this pipeline.** In `config/config.yaml`, fill in for each species:

```yaml
references:
  Delphinapterus_leucas:
    fasta: /mnt/databases/beluga/beluga.fa
    contigs: [<EXACT_AUTOSOME_NAME_1>, <EXACT_AUTOSOME_NAME_2>]
  Ursus_maritimus:
    fasta: /mnt/databases/polar_bear/polar_bear.fa
    contigs: [<EXACT_AUTOSOME_NAME_1>, <EXACT_AUTOSOME_NAME_2>]
```

The bracketed contigs are placeholders: **use actual accession/contig IDs**, in the **same order** as the FASTA index. Include every autosome relevant to your analysis, not just two. The script checks names, lengths and order against both `.fai` and `.dict` quickly (without reading multi-GB FASTA contents).

Required files **already present** before starting:

```text
reference.fa
reference.fa.fai
reference.dict           # GATK sequence dictionary; note NOT reference.fa.dict
reference.fa.amb
reference.fa.ann
reference.fa.bwt
reference.fa.pac
reference.fa.sa
```

If any index is missing, stop and prepare it **once**, as a separate administrative operation outside the processing pipeline, or choose another preindexed reference. The pipeline itself refuses to proceed and never automatically creates it.

**Important:** BAM @SQ name/length equality does not guarantee that the exact reference sequence is identical. Source assembly IDs and reference checksum provenance should be checked when collecting WGS; the BAM-validator also compares `M5` MD5 values when both BAM and `.dict` headers supply them.

## 4. Sample table (WGS only)

Replace `config/samples.tsv` with real samples. It is initially header-only by design; see `config/samples.example.tsv` for fictitious rows.

Mandatory tab-separated columns:

```text
sample   species   input_type   library   lane   r1   r2   bam   cram
```

- `sample`: globally unique across all species, ASCII letters/digits/underscores/dots/dashes.
- `species`: matches the exact key under YAML `references`.
- `input_type`: `fastq`, `bam`, or `cram`; **no** low-pass or highcov branch.
- `library` and `lane`: used as BWA read-group IDs for FASTQ; provide identifiers also for aligned input rows (placeholders `NA` acceptable).
- FASTQ: one row per sequencing lane, `r1` + `r2` supplied. Multiple lanes from **one library** of the same sample are supported; multiple *independent* libraries in one sample are not (the validator fails rather than silently mis-marking duplicates).
- BAM: one row per sample, fill `bam` only; must be coordinate sorted, contain correct `@RG SM`, match reference, be indexed externally, and already have appropriate duplicate flags. Source file is used **in place**.
- CRAM: similarly, fill only `cram`; validate/reference-decode and **transcode once to BAM**, leaving original CRAM untouched. No repeat duplicate marking.

Paths can be absolute or relative to **project root**, not relative to the TSV file location. `workdir` should be an absolute path to a shared high-capacity directory **outside** the code repository.

## 5. Choose processing settings

`config/config.yaml`:

```yaml
preprocessing:
  trim_fastq: false   # default: directly align raw paired FASTQ

duplicate_marking:
  optical_distance: 0 # don't classify optical duplicates until instrument is known

gatk:
  hc_threads: 4
  hc_heap_gb: 16
  db_heap_gb: 24
  genotype_heap_gb: 12
  db_batch_size: 50
  db_readers: 1

slurm:
  max_array_parallel: 4
  # partition: main
  # account: users
  # qos: normal
```

`fastp` is optional. For some Illumina runs, it may be useful for adapters/quality; for already cleaned inputs it causes unnecessary full FASTQ rewrites. Default: **disabled**. Read QC can be examined before enabling it on a cohort.

SLURM uses fixed, conservative **starter requests per stage**, declared in `scripts/pipeline.py`. You can optionally override `cpus`, `mem_mb`, `hours` by stage using `slurm.resources`. **No profiling-driven SLURM optimization is included.** Avoid requesting fewer CPUs/RAM than GATK thread count/Java heap.

## 6. Validate, preview, submit

Before WGS files are downloaded, you can preview ONLY if you provide sample rows and the real contig names:

```bash
bash run_pipeline.sh plan --config config/config.yaml --allow-missing-inputs
```

This does **not** check files/indices; it is purely a task-planning preview.

After downloading data / selecting existing references and filling inputs:

```bash
conda activate wgs-gatk
bash run_pipeline.sh validate --config config/config.yaml
bash run_pipeline.sh plan     --config config/config.yaml
bash run_pipeline.sh submit   --config config/config.yaml
```

`submit` builds a frozen task plan and submits one SLURM **job array per stage** with `--dependency=afterok:<jobid>` between stages. A downstream stage waits until required upstream *arrays* complete, not just a single finished sample. `--array %N` limits simultaneous jobs per stage.

To observe:

```bash
squeue -u "$USER"
bash run_pipeline.sh status --config config/config.yaml
# stderr / stdout under <workdir>/jobs/logs/
```

**Resume after failures:** once jobs from the preceding attempt are no longer pending/running, rerun `submit`. Successful stages are skipped using checkpoint signatures and expected output files; unfinished and downstream tasks are resubmitted. If one task in a SLURM array fails, dependent arrays remain blocked; diagnose failed job logs before resubmitting. To restart everything from scratch: `submit --force` **only after** ensuring there are no active jobs.

For an extremely small *test* data set, `local` runs stages sequentially, without submitting to SLURM. Do not run multi-GB WGS on a shared login node; test `local` inside an interactive allocated compute session.

## 7. Outputs and working layout

```text
workdir/
├── bam/<sample>.bam[.bai]
├── lanes/<sample>/<library>.<lane>.bam
├── gvcf/<species>/<sample>/<contig>.g.vcf.gz
├── genomicsdb/<species>/<contig>/
├── joint/<species>/<contig>.vcf.gz
├── vcf/<species>/cohort.raw.vcf.gz
├── vcf/<species>/cohort.raw.vcf.gz.csi
├── qc/fastp/ (if trimming enabled)
├── qc/duplicates/
├── qc/bam/<sample>.{flagstat.txt,stats.txt,...}
├── qc/input_validation/
├── qc/multiqc/multiqc_report.html
├── jobs/{scripts,logs,sample_maps}/
├── checkpoints/
└── pipeline_plan.json
```

`cohort.raw.vcf.gz` is a **raw joint-genotyped variants VCF**. It is **not filtered SNP data**. GenotypeGVCFs can emit SNPs and INDELs, and non-PASS sites are not currently removed. Only the contigs in the configured list are included. No separate `population_labels.tsv` is generated here.

## 8. Incremental cohorts and reruns

New sample added to `samples.tsv`:

- A new per-lane alignment and per-sample duplicate-marked BAM is created (or an external alignment validated).
- New per-sample × contig GVCF is created.
- Existing gVCFs are reused if their sources, command signature and outputs are unchanged.
- GenomicsDB and downstream per-contig VCFs **are rebuilt for the new cohort**, as expected. It is not a fully incremental GenomicsDB *update* implementation; this is simpler and less error-prone for the initial study.

**Important:** task signatures cover generated shell commands, source script code and GenomicsDB sample-map content. Versions of tools/environment and checksums of all data are not captured automatically; archive `envs/wgs.yaml`, resolved package list, FASTA assembly/accession and relevant checksums in your analysis provenance.

## 9. Limitations before production

1. The framework is tested by Python unit tests and shell syntax checks but **not** by a full GATK/BWA run in this container. Perform a small real FASTQ → BAM → gVCF → joint VCF smoke test on the cluster first.
2. For large cohorts, whole-chromosome GenomicsDB workspaces and GATK Java heaps may require manual resource adjustments. Stage arrays deliberately use simple dependencies, not a profiling-based scheduler.
3. Input BAMs and CRAMs are expected to have been correctly prepared and duplicate-marked already. Cross-dataset library chemistry, sequencing platform, reference assembly and duplicate policy need checking during data collection.
4. Current pipeline does not perform BQSR, QC-based sample exclusion, VQSR/hard filtering, sample contamination analysis, genotype/SNP QC, PCA, LD pruning or SNP-panel selection.
5. Diploid HaplotypeCaller model is appropriate for initial **autosomal** WGS. Add sex-chromosome/ploidy strategies separately for each species and sex.
6. Gene annotations, downloaded FASTQ, external references and index preparation are **outside** the scope of this code; no downloading or indexing happens automatically.

### Tool documentation

- [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/30332006386459-HaplotypeCaller)
- [GATK GenomicsDBImport](https://gatk.broadinstitute.org/hc/en-us/articles/360036889212-GenomicsDBImport)
- [GATK GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360036731211-GenotypeGVCFs)
- [Samtools duplicate marking algorithm](https://www.htslib.org/algorithms/duplicate.html)

The first server exercise should use **1–2 WGS from a single species, only 1 short representative autosome** in the test config, keeping a separate temporary workdir; after the smoke test, use the full autosomal contig list and the entire species cohort. This prevents a full-genome run being launched before scripts and dependencies are confirmed on the specific SLURM cluster.
