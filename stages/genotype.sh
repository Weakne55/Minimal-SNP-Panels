#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/resource_lib.bash"

ref=$1; db=$2; out=$3; shift 3
threads=${MIN_SNP_THREADS:?}
heap=$(minsnp_heap_mb 75 16384)

intervals=()
for contig in "$@"; do intervals+=(-L "$contig"); done
mkdir -p "$(dirname "$out")"
[[ ! -e "$out" ]] || { echo "Output exists: $out" >&2; exit 2; }

tmp_vcf="${out}.incomplete.vcf.gz"
tmp_bcf="${out}.incomplete.bcf"
trap 'rm -f -- "$tmp_vcf" "${tmp_vcf}.tbi" "${tmp_vcf}.csi" "$tmp_bcf" "${tmp_bcf}.csi"' EXIT

gatk --java-options "-Xmx${heap}m -XX:ActiveProcessorCount=${threads}" \
    GenotypeGVCFs -R "$ref" -V "gendb://${db}" "${intervals[@]}" -O "$tmp_vcf"

bcftools view -Ob -o "$tmp_bcf" "$tmp_vcf"
bcftools index -c "$tmp_bcf"

mv -- "$tmp_bcf" "$out"
mv -- "${tmp_bcf}.csi" "${out}.csi"
