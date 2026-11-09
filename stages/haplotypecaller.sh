#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/resource_lib.bash"

ref=$1; bam=$2; out=$3; threads=$4; shift 4
heap=$(minsnp_heap_mb 70 16384)
intervals=()
for contig in "$@"; do intervals+=(-L "$contig"); done
mkdir -p "$(dirname "$out")"
[[ ! -e "$out" ]] || { echo "Output exists: $out" >&2; exit 2; }

gatk --java-options "-Xmx${heap}m -XX:ActiveProcessorCount=${threads}" \
    HaplotypeCaller -R "$ref" -I "$bam" "${intervals[@]}" \
    -ERC GVCF --native-pair-hmm-threads "$((threads < 8 ? threads : 8))" -O "$out"
