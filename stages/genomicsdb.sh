#!/usr/bin/env bash
# Назначение: импортировать gVCF всех образцов для одной хромосомы.
# Вход: FASTA, contig, двухколоночный GATK sample-name-map (без заголовка).
# Выход: новый GenomicsDB workspace и маркер .completed.
# Важно: готовую GenomicsDB НЕ удаляем при повторном запуске.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: genomicsdb.sh REF.fa CONTIG SAMPLE_MAP.tsv NEW_DB_DIR

GATK GenomicsDBImport: один workspace для одного contig.
SAMPLE_MAP: SAMPLE_ID<TAB>/absolute/path/sample.contig.g.vcf.gz
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then usage; exit 0; fi
if (( $# != 4 )); then usage >&2; exit 2; fi

reference_fasta=$1
contig=$2
sample_map=$3
output_db=$4

if [[ -e "$output_db" ]]; then
  echo "GenomicsDB уже существует: $output_db" >&2
  exit 2
fi
mkdir -p "$(dirname "$output_db")"

gatk --java-options '-Xmx24g' GenomicsDBImport \
  -R "$reference_fasta" \
  -L "$contig" \
  --sample-name-map "$sample_map" \
  --validate-sample-name-map true \
  --genomicsdb-workspace-path "$output_db"

test -d "$output_db"
printf 'complete\n' > "${output_db}/.completed"
