#!/usr/bin/env bash
# Назначение: собрать метрики качества выравнивания и глубины покрытия.
# Вход: окончательный BAM, префикс файлов QC, число потоков.
# Выход: flagstat, samtools stats, mosdepth summary / dist.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: bam_qc.sh INPUT.bam OUTPUT_PREFIX THREADS

Метрики выравнивания: samtools flagstat, samtools stats и mosdepth.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then usage; exit 0; fi
if (( $# != 3 )); then usage >&2; exit 2; fi

input_bam=$1
output_prefix=$2
threads=$3

mkdir -p "$(dirname "$output_prefix")"

samtools flagstat -@ "$threads" "$input_bam" > "${output_prefix}.flagstat.txt"
samtools stats -@ "$threads" "$input_bam" > "${output_prefix}.stats.txt"
mosdepth -t "$threads" -n "$output_prefix" "$input_bam"
