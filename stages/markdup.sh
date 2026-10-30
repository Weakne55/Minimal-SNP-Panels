#!/usr/bin/env bash
# Назначение: отметить (не удалять) дубликаты после объединения sequencing runs.
# Usage: markdup.sh INPUT.bam OUT.bam METRICS.txt THREADS
# Требует BAM после samtools fixmate -m и coordinate sort.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/resource_lib.bash"
if [[ "${1:-}" == '-h' || "${1:-}" == '--help' ]]; then
    echo 'Usage: markdup.sh INPUT.bam OUT.bam METRICS.txt THREADS'; exit 0
fi
if (( $# != 4 )); then echo 'Usage: markdup.sh INPUT.bam OUT.bam METRICS.txt THREADS' >&2; exit 2; fi
input_bam=$1
output_bam=$2
metrics_file=$3
threads=$4
minsnp_check_threads "$threads"
[[ ! -e "$output_bam" ]] || { echo "BAM уже существует: $output_bam" >&2; exit 2; }
samtools quickcheck -v "$input_bam"
mkdir -p "$(dirname "$output_bam")" "$(dirname "$metrics_file")"
tmp_bam="${output_bam}.incomplete.bam"
trap 'rm -f -- "$tmp_bam" "${tmp_bam}.bai"' EXIT
workers=$(minsnp_samtools_workers "$threads")
samtools markdup -@ "$workers" -f "$metrics_file" "$input_bam" "$tmp_bam"
samtools quickcheck -v "$tmp_bam"
samtools index -@ "$workers" -b "$tmp_bam"
# Публикуем оба файла только после успешной проверки.
mv -- "$tmp_bam" "$output_bam"
mv -- "${tmp_bam}.bai" "${output_bam}.bai"
trap - EXIT
