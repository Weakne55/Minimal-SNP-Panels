#!/usr/bin/env bash
# Usage: bash scripts/pull_sif.sh docker://ghcr.io/weakne55/minimal-snp-panels-wgs:snakemake [output.sif]
set -euo pipefail
image=${1:?Usage: pull_sif.sh OCI_URI [OUTPUT.sif]}
output=${2:-containers/wgs-gatk-mem3.sif}
if command -v apptainer >/dev/null; then runtime=apptainer
elif command -v singularity >/dev/null; then runtime=singularity
else echo 'Не найден apptainer или singularity' >&2; exit 1
fi
mkdir -p "$(dirname "$output")"
# pull of a pre-built OCI image: no fakeroot / setuid root build is necessary.
"$runtime" pull --force "$output" "$image"
"$runtime" exec "$output" sh -c 'set -e; command -v bwa-mem3; command -v gatk; command -v bcftools; command -v samtools'
