#!/usr/bin/env bash
set -euo pipefail
(( $# >= 2 )) || { echo 'Usage: concat.sh OUTPUT_VCF INPUT_CHR_VCF [INPUT_CHR_VCF...]' >&2; exit 2; }
out=$1; shift
mkdir -p "$(dirname "$out")"
if (( $# == 1 )); then
  bcftools view -Oz -o "$out" "$1"
else
  bcftools concat -Oz -o "$out" "$@"
fi
bcftools index -f -c "$out"
