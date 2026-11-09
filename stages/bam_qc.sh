#!/usr/bin/env bash
set -euo pipefail

bam=$1; prefix=$2; threads=$3
mkdir -p "$(dirname "$prefix")"
workers=$(( threads > 1 ? threads - 1 : 0 ))

samtools flagstat -@ "$workers" "$bam" > "${prefix}.flagstat.txt"
samtools stats -@ "$workers" "$bam" > "${prefix}.stats.txt"
mosdepth -t "$threads" -n "$prefix" "$bam"
