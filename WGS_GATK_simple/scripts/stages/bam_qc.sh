#!/usr/bin/env bash
set -euo pipefail
[[ $# == 3 ]] || { echo 'Usage: bam_qc.sh BAM PREFIX THREADS' >&2; exit 2; }
bam=$1 prefix=$2 threads=$3
mkdir -p "$(dirname "$prefix")"
samtools flagstat -@ "$threads" "$bam" > "${prefix}.flagstat.txt"
samtools stats -@ "$threads" "$bam" > "${prefix}.stats.txt"
mosdepth --threads "$threads" --mapq 20 -n --by 1000000 --thresholds 1,5,10,20 "$prefix" "$bam"
