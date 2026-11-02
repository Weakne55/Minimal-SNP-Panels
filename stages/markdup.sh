#!/usr/bin/env bash
# DupBlaster требует QNAME-grouped BAM; merged BAM уже coordinate-sorted.
# Usage: markdup.sh INPUT.bam OUT.bam METRICS_PREFIX THREADS
# Дубликаты отмечаем (FLAG 0x400), но не удаляем.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/resource_lib.bash"
if [[ "${1:-}" == '-h' || "${1:-}" == '--help' ]]; then
    echo 'Usage: markdup.sh INPUT.bam OUT.bam METRICS_PREFIX THREADS'; exit 0
fi
if (( $# != 4 )); then echo 'Usage: markdup.sh INPUT.bam OUT.bam METRICS_PREFIX THREADS' >&2; exit 2; fi
input_bam=$1
output_bam=$2
metrics_prefix=$3
threads=$4
minsnp_check_threads "$threads"
[[ ! -e "$output_bam" ]] || { echo "BAM уже существует: $output_bam" >&2; exit 2; }
samtools quickcheck -v "$input_bam"
mkdir -p "$(dirname "$output_bam")" "$(dirname "$metrics_prefix")"
memory_mb=$(minsnp_memory_mb)
# collate (главный поток) + dupblaster (worker и IO) + samtools sort.
# Ограничиваем sort потоками оставшегося бюджета; слабый CPU делится по времени.
sort_workers=$((threads > 5 ? threads - 5 : 0))
(( sort_workers > 16 )) && sort_workers=16
sort_memory_mb=$((memory_mb / 6 / (sort_workers + 1)))
(( sort_memory_mb > 768 )) && sort_memory_mb=768
(( sort_memory_mb < 64 )) && sort_memory_mb=64
scratch_root=${SLURM_TMPDIR:-${TMPDIR:-$(dirname "$output_bam")}}
[[ -d "$scratch_root" ]] || { echo "Нет временного каталога: $scratch_root" >&2; exit 2; }
workdir=$(mktemp -d "${scratch_root%/}/minsnp-dupblaster.XXXXXXXX")
tmp_bam="${output_bam}.incomplete.bam"
trap 'rm -rf -- "$workdir"; rm -f -- "$tmp_bam" "${tmp_bam}.bai"' EXIT
# Optical duplicate classification выключена: SRA может изменять формат QNAME.
# Маркировка PCR/library duplicates остаётся включённой.
samtools collate -O -u "$input_bam" "$workdir/collate" \
    | dupblaster --sequencing-duplicate-detection off \
        --single-end-strategy picard-exact \
        --metrics-prefix "$metrics_prefix" -o - \
    | samtools sort -@ "$sort_workers" -m "${sort_memory_mb}M" \
        -T "$workdir/sort" -o "$tmp_bam" -
samtools quickcheck -v "$tmp_bam"
[[ -s "${metrics_prefix}.duplicate-metrics.tsv" ]] || {
    echo "DupBlaster не создал метрики: ${metrics_prefix}.duplicate-metrics.tsv" >&2
    exit 1
}
samtools index -@ "$sort_workers" -b "$tmp_bam"
[[ -s "${tmp_bam}.bai" ]] || { echo "Нет индекса BAM" >&2; exit 1; }
mv -- "$tmp_bam" "$output_bam"
mv -- "${tmp_bam}.bai" "${output_bam}.bai"
trap - EXIT
rm -rf -- "$workdir"
