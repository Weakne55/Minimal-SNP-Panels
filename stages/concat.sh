#!/usr/bin/env bash
# Назначение: конкатенация VCF разных contigs с одинаковыми sample columns.
# Usage: concat.sh OUT.bcf IN1.vcf.gz [IN2.vcf.gz ...]
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/resource_lib.bash"
if [[ "${1:-}" == '-h' || "${1:-}" == '--help' ]]; then
    echo 'Usage: concat.sh OUT.bcf IN1.vcf.gz [IN2.vcf.gz ...]'; exit 0
fi
if (( $# < 2 )); then echo 'Usage: concat.sh OUT.bcf IN1.vcf.gz [IN2.vcf.gz ...]' >&2; exit 2; fi
output_bcf=$1
shift
[[ ! -e "$output_bcf" && ! -e "${output_bcf}.csi" ]] || {
    echo "Выходной BCF уже существует: $output_bcf" >&2; exit 2;
}
mkdir -p "$(dirname "$output_bcf")"
first_samples=$(bcftools query -l "$1")
for part in "$@"; do
    [[ "$(bcftools query -l "$part")" == "$first_samples" ]] || {
        echo "Несовпадающий состав/порядок образцов в $part" >&2; exit 2;
    }
done
threads=${MIN_SNP_THREADS:-1}
minsnp_check_threads "$threads"
tmp_bcf="${output_bcf}.incomplete.bcf"
trap 'rm -f -- "$tmp_bcf" "${tmp_bcf}.csi"' EXIT
compression_threads=$(minsnp_samtools_workers "$threads")
(( compression_threads > 16 )) && compression_threads=16
bcftools concat --threads "$compression_threads" -Ob -o "$tmp_bcf" "$@"
bcftools index -f "$tmp_bcf"
bcftools index -n "$tmp_bcf" >/dev/null
mv -- "$tmp_bcf" "$output_bcf"
mv -- "${tmp_bcf}.csi" "${output_bcf}.csi"
trap - EXIT
