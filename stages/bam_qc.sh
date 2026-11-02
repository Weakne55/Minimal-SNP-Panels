#!/usr/bin/env bash
# BAM QC: samtools flagstat/stats и PanDepth по всем contigs.
# Usage: bam_qc.sh INPUT.bam OUTPUT_PREFIX THREADS
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/resource_lib.bash"
if [[ "${1:-}" == '-h' || "${1:-}" == '--help' ]]; then
    echo 'Usage: bam_qc.sh INPUT.bam OUTPUT_PREFIX THREADS'; exit 0
fi
if (( $# != 3 )); then echo 'Usage: bam_qc.sh INPUT.bam OUTPUT_PREFIX THREADS' >&2; exit 2; fi
input_bam=$1
prefix=$2
threads=$3
minsnp_check_threads "$threads"
samtools quickcheck -v "$input_bam"
mkdir -p "$(dirname "$prefix")"
workers=$(minsnp_samtools_workers "$threads")
# Сначала пишем служебные метрики; результат PanDepth проверяем отдельно.
samtools flagstat -@ "$workers" "$input_bam" > "${prefix}.flagstat.txt"
samtools stats -@ "$workers" "$input_bam" > "${prefix}.stats.txt"
pandepth -i "$input_bam" -o "$prefix" -t "$threads"
[[ -s "${prefix}.chr.stat.gz" ]] || {
    echo "PanDepth не создал сводку покрытия: ${prefix}.chr.stat.gz" >&2
    exit 1
}
gzip -t "${prefix}.chr.stat.gz"
