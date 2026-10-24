#!/usr/bin/env bash
set -euo pipefail
for program in python3 bwa samtools bcftools gatk mosdepth multiqc; do
  command -v "$program" >/dev/null || { echo "MISSING $program" >&2; exit 1; }
  printf '%-10s %s\n' "$program" "$(command -v "$program")"
done
python3 -c 'import yaml; print("PyYAML", yaml.__version__)'
command -v sbatch >/dev/null && echo 'SLURM sbatch found' || echo 'sbatch not found (local tests only)'
java -version 2>&1 | head -2
