#!/usr/bin/env bash
# Назначение: BAM QC (flagstat, stats и mosdepth).
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
mkdir -p "$(dirname "$prefix")"
workers=$(minsnp_samtools_workers "$threads")
samtools flagstat -@ "$workers" "$input_bam" > "${prefix}.flagstat.txt"
samtools stats -@ "$workers" "$input_bam" > "${prefix}.stats.txt"
mosdepth -t "$threads" -n "$prefix" "$input_bam"
