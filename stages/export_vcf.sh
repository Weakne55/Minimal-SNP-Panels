#!/usr/bin/env bash
set -euo pipefail

bcf=$1; out=$2
threads=${MIN_SNP_THREADS:?}
mkdir -p "$(dirname "$out")"
[[ ! -e "$out" ]] || { echo "Output exists: $out" >&2; exit 2; }

tmp="${out}.incomplete.vcf.gz"
trap 'rm -f -- "$tmp" "${tmp}.csi"' EXIT

bcftools view --threads "$((threads > 1 ? threads - 1 : 0))" -Oz -o "$tmp" "$bcf"
bcftools index -c "$tmp"

mv -- "$tmp" "$out"
mv -- "${tmp}.csi" "${out}.csi"
