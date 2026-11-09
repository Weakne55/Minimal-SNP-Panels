#!/usr/bin/env bash
set -euo pipefail

bam=$1; out=$2; metrics=$3; threads=$4
workers=$(( threads > 1 ? threads - 1 : 0 ))
mkdir -p "$(dirname "$out")" "$(dirname "$metrics")"
[[ ! -e "$out" ]] || { echo "Output exists: $out" >&2; exit 2; }

tmp="${out}.incomplete.bam"
trap 'rm -f -- "$tmp" "${tmp}.bai"' EXIT

samtools markdup -@ "$workers" -f "$metrics" "$bam" "$tmp"
samtools quickcheck -v "$tmp"
samtools index -@ "$workers" -b "$tmp"

mv -- "$tmp" "$out"
mv -- "${tmp}.bai" "${out}.bai"
