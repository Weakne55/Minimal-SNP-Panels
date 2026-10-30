#!/usr/bin/env bash
# Назначение: нормализовать multi-sample VCF и отфильтровать варианты.
# Вход: FASTA, VCF/BCF, пороги генотипов и SNP, аутосомные contigs.
# Выход: нормализованный BCF и BCF после технического QC.
# Это технический QC, не отбор population-informative SNP.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: qc.sh REF.fa INPUT.vcf.gz NORMALIZED.bcf FILTERED.bcf MIN_DP MIN_GQ MAX_MISSING MIN_MAF CONTIGS

CONTIGS: имена хромосом через запятую, например CM187741.1,CM187742.1.
MIN_DP / MIN_GQ: пороги для маскирования низкокачественных GT.
MAX_MISSING / MIN_MAF: предельная доля пропусков / минимальная MAF.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then usage; exit 0; fi
if (( $# != 9 )); then usage >&2; exit 2; fi

reference_fasta=$1
input_vcf=$2
normalized_bcf=$3
filtered_bcf=$4
min_dp=$5
min_gq=$6
max_missing=$7
min_maf=$8
contigs=$9

mkdir -p "$(dirname "$normalized_bcf")" "$(dirname "$filtered_bcf")"

# 1. Оставляем указанные аутосомы и приводим REF/ALT к FASTA.
autosomes_bcf="${normalized_bcf}.autosomes.bcf"
bcftools view -r "$contigs" -Ob -o "$autosomes_bcf" "$input_vcf"
bcftools norm -f "$reference_fasta" -m -any \
  -Ob -o "${normalized_bcf}.tmp.bcf" "$autosomes_bcf"
mv -- "${normalized_bcf}.tmp.bcf" "$normalized_bcf"
rm -f -- "$autosomes_bcf"
bcftools index -f "$normalized_bcf"

# 2. Если в заголовке есть FORMAT/DP и FORMAT/GQ, маскируем плохие GT.
# Не требуем отсутствующих полей: некоторые VCF их не содержат.
header=$(bcftools view -h "$normalized_bcf")
mask_expression=""
if grep -q '^##FORMAT=<ID=DP,' <<< "$header"; then
  mask_expression="FMT/DP<${min_dp}"
fi
if grep -q '^##FORMAT=<ID=GQ,' <<< "$header"; then
  if [[ -n "$mask_expression" ]]; then
    mask_expression+=" || FMT/GQ<${min_gq}"
  else
    mask_expression="FMT/GQ<${min_gq}"
  fi
fi

source_bcf=$normalized_bcf
if [[ -n "$mask_expression" ]]; then
  source_bcf="${filtered_bcf}.masked.bcf"
  bcftools +setGT "$normalized_bcf" -Ob -o "$source_bcf" -- \
    -t q -n . -i "$mask_expression"
fi

# 3. После маскирования GT обновляем частоты аллелей (AC/AN/AF/MAF).
tagged_bcf="${filtered_bcf}.tags.bcf"
bcftools +fill-tags "$source_bcf" -Ob -o "$tagged_bcf" -- -t AC,AN,AF,MAF

# 4. Оставляем качественные биаллельные SNP с нужной MAF и долей пропусков.
bcftools view -m2 -M2 -v snps -f PASS,. \
  -i "F_MISSING<=${max_missing} && MAF>=${min_maf}" \
  -Ob -o "${filtered_bcf}.tmp.bcf" "$tagged_bcf"
mv -- "${filtered_bcf}.tmp.bcf" "$filtered_bcf"
bcftools index -f "$filtered_bcf"

rm -f -- "$tagged_bcf" "${filtered_bcf}.masked.bcf"
