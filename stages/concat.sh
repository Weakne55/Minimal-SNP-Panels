#!/usr/bin/env bash
# Назначение: объединить непересекающиеся хромосомные multi-sample VCF.
# Вход: VCF/BCF всех contigs одного вида и сборки, в порядке референса.
# Выход: общий multi-sample BCF с индексом CSI.
# Важно: состав и порядок образцов должны совпадать во всех VCF.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: concat.sh OUT.bcf CONTIG_1.vcf.gz [CONTIG_2.vcf.gz ...]

bcftools concat объединяет genomic regions, но НЕ разные наборы образцов.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then usage; exit 0; fi
if (( $# < 2 )); then usage >&2; exit 2; fi

output_bcf=$1
shift
input_vcfs=("$@")

if [[ -e "$output_bcf" || -e "${output_bcf}.csi" ]]; then
  echo "Результат уже существует: $output_bcf" >&2
  exit 2
fi
mkdir -p "$(dirname "$output_bcf")"

# Проверяем sample columns до объединения файлов.
first_samples=$(bcftools query -l "${input_vcfs[0]}")
for input_vcf in "${input_vcfs[@]}"; do
  current_samples=$(bcftools query -l "$input_vcf")
  if [[ "$current_samples" != "$first_samples" ]]; then
    echo "Разный состав или порядок образцов: $input_vcf" >&2
    exit 2
  fi
done

# Сначала создаём и проверяем временный BCF, затем публикуем результат.
temporary_bcf="${output_bcf}.incomplete.bcf"
trap 'rm -f -- "$temporary_bcf" "${temporary_bcf}.csi"' EXIT

bcftools concat -Ob -o "$temporary_bcf" "${input_vcfs[@]}"
bcftools index -f "$temporary_bcf"
bcftools index -n "$temporary_bcf" >/dev/null
mv -- "$temporary_bcf" "$output_bcf"
mv -- "${temporary_bcf}.csi" "${output_bcf}.csi"
trap - EXIT
