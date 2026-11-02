# Full WGS acquisition — historical helper

This directory contains the earlier `prepare_full_wgs.py` downloader for a small initial, full-read beluga/polar-bear cohort. **It is not the orchestrator of the current 20-sample pilot** and is not a prerequisite for the reusable VCF/BCF → SNP-panel analysis.

The ongoing pilot `pilot_20261030` was submitted through `wgs_pilot.py` (36 verified runs; 10 belugas and 10 polar bears). Its generated SLURM scripts still invoke the preserved server-local path `tools/wgs_pilot.py`. **Do not run `prepare_full_wgs.py download` or any other downloader for the same files while the pilot is active.** Do not rename or remove the local `tools/` until those jobs finish.

Historical `prepare_full_wgs.py` uses an ENA `discover` → frozen `ena_selection.json` → `download` model with size/MD5 checking. It was written for the previous `WGS_Bash_GATK` directory layout and `wgs-gatk` environment. Its old generated configs **are not** the current `src/minsnp` input contract. Treat the script as an archival utility that needs a separate review before reuse.

For the active pipeline, read [repository README](../../README.md) and [technical guide](../../docs/PROJECT_GUIDE.md). The current FASTQ TSV contract is `sample_id`, `species`, `run_id`, `library`, `fastq_r1`, `fastq_r2`. Real example manifests live in `inputs/`.
