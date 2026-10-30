#!/usr/bin/env bash
# Назначение: экспортировать совместный BCF в сжатый VCF.gz.
# Вход: индексированный BCF после joint genotyping.
# Выход: VCF.gz с индексом CSI.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: export_vcf.sh INPUT.bcf OUTPUT.vcf.gz

Экспорт bcftools view -Oz; исходный BCF остаётся без изменений.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then usage; exit 0; fi
if (( $# != 2 )); then usage >&2; exit 2; fi

input_bcf=$1
output_vcf=$2

if [[ -e "$output_vcf" || -e "${output_vcf}.csi" ]]; then
  echo "Выходной VCF уже существует: $output_vcf" >&2
  exit 2
fi
mkdir -p "$(dirname "$output_vcf")"

# Временный файл и индекс не остаются после неуспешного экспорта.
temporary_vcf="${output_vcf}.incomplete.vcf.gz"
trap 'rm -f -- "$temporary_vcf" "${temporary_vcf}.csi"' EXIT

bcftools view -Oz -o "$temporary_vcf" "$input_bcf"
bcftools index -c "$temporary_vcf"
bcftools index -n "$temporary_vcf" >/dev/null
mv -- "$temporary_vcf" "$output_vcf"
mv -- "${temporary_vcf}.csi" "${output_vcf}.csi"
trap - EXIT
