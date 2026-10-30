#!/usr/bin/env bash
# Назначение: экспорт cohort.raw.bcf в индексированный bgzip VCF.
# Usage: export_vcf.sh INPUT.bcf OUTPUT.vcf.gz
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/resource_lib.bash"
if [[ "${1:-}" == '-h' || "${1:-}" == '--help' ]]; then
    echo 'Usage: export_vcf.sh INPUT.bcf OUTPUT.vcf.gz'; exit 0
fi
if (( $# != 2 )); then echo 'Usage: export_vcf.sh INPUT.bcf OUTPUT.vcf.gz' >&2; exit 2; fi
input_bcf=$1
output_vcf=$2
[[ ! -e "$output_vcf" && ! -e "${output_vcf}.csi" ]] || {
    echo "Выходной VCF уже существует: $output_vcf" >&2; exit 2;
}
mkdir -p "$(dirname "$output_vcf")"
threads=${MIN_SNP_THREADS:-1}
minsnp_check_threads "$threads"
tmp_vcf="${output_vcf}.incomplete.vcf.gz"
trap 'rm -f -- "$tmp_vcf" "${tmp_vcf}.csi"' EXIT
compression_threads=$(minsnp_samtools_workers "$threads")
(( compression_threads > 16 )) && compression_threads=16
bcftools view --threads "$compression_threads" -Oz -o "$tmp_vcf" "$input_bcf"
bcftools index -c "$tmp_vcf"
bcftools index -n "$tmp_vcf" >/dev/null
mv -- "$tmp_vcf" "$output_vcf"
mv -- "${tmp_vcf}.csi" "${output_vcf}.csi"
trap - EXIT
