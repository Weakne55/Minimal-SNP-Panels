#!/usr/bin/env bash
set -euo pipefail

r1=$1; r2=$2; out1=$3; out2=$4; threads=$5
mkdir -p "$(dirname "$out1")" "$(dirname "$out2")"

fastp -i "$r1" -I "$r2" -o "$out1" -O "$out2" \
    -w "$((threads < 16 ? threads : 16))" \
    -h "${out1%.fastq.gz}.html" -j "${out1%.fastq.gz}.json"
