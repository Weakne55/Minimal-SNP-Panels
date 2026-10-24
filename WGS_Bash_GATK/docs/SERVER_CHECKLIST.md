# Server acceptance checklist

## Stage A — before downloading WGS

- [ ] Unzip into a shared path accessible from SLURM compute nodes.
- [ ] Create and activate Conda environment; `bash scripts/check_stack.sh`.
- [ ] Run `python -m unittest discover -s tests -v`.
- [ ] Locate existing reference `.fa/.fai/.dict` and BWA-MEM indexes; do not duplicate them.
- [ ] Choose exactly matched reference assemblies per species and record accession and reference checksum.
- [ ] Create a complete autosome list based on `.fai` (start with 1 autosome for smoke test).
- [ ] Set workdir outside repository, on shared storage; confirm quota/free space.
- [ ] Configure SLURM partition/account if required. No resource profiling work is planned.

## Stage B — once representative FASTQ/CRAM/BAM are available

- [ ] Fill `config/samples.tsv` with a small subset (1–2 samples, one species).
- [ ] Check correct data type, library/lane IDs, paired FASTQs and absence of duplicate sample IDs.
- [ ] Validate: `bash run_pipeline.sh validate --config config/config.yaml`.
- [ ] Preview: `bash run_pipeline.sh plan --config config/config.yaml`.
- [ ] Submit: `bash run_pipeline.sh submit --config config/config.yaml`.
- [ ] Inspect `squeue -u "$USER"`, `workdir/jobs/logs`, checkpoint state.
- [ ] Verify BAM/BAI, duplicate stats, per-sample gVCF/index, GenomicsDB, per-contig VCF, cohort raw VCF/index.
- [ ] Confirm sample count in cohort VCF, VCF headers, reference contig name, nonempty variant sites, and expected coverage.
- [ ] Record versions, date, reference assembly and first test outputs in a reproducibility log.

## Stage C — after smoke test passes

- [ ] Expand `contigs` to all selected autosomes; delete/replace test-only config or use separate workdir.
- [ ] Add full sample manifests for both species.
- [ ] Run species WGS → gVCF → joint VCF; no QC-gate/SNP QC at this stage.
- [ ] Add separate downstream modules for population structure and SNP-panel selection later.

## If a SLURM array fails

1. Examine `jobs/logs/<stage>_<JOBID>_<array-index>.err`.
2. Fix missing indexes, incompatible ref, GATK memory or bad inputs.
3. Ensure no prior SLURM jobs from the run are still active (`squeue`).
4. Reissue `submit` (does not overwrite successful upstream products).
5. Avoid `--force` unless a complete rerun is deliberately required.

Contig-list helper (read-only):

```bash
python scripts/list_reference_contigs.py /mnt/databases/beluga/beluga.fa.fai
```

Copy selected **autosomal** contig names into `config/config.yaml` in `.fai` order. This is not a reference preparation/indexing command.
