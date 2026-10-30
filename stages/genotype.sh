#!/usr/bin/env bash
# Назначение: получить совместные генотипы всех образцов для одного contig.
# Вход: FASTA, готовый GenomicsDB workspace.
# Выход: multi-sample VCF.gz по одной хромосоме (с индексом GATK).
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: genotype.sh REF.fa GENOMICSDB_DIR OUT.vcf.gz

GATK GenotypeGVCFs: gendb:// → обычный хромосомный VCF.gz.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then usage; exit 0; fi
if (( $# != 3 )); then usage >&2; exit 2; fi

reference_fasta=$1
input_db=$2
output_vcf=$3

mkdir -p "$(dirname "$output_vcf")"

gatk --java-options '-Xmx16g' GenotypeGVCFs \
  -R "$reference_fasta" \
  -V "gendb://${input_db}" \
  -O "$output_vcf"

test -s "$output_vcf"
