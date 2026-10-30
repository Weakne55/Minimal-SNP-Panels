#!/usr/bin/env bash
# Назначение: получить gVCF одного образца на одной хромосоме.
# Вход: референс FASTA, final markdup BAM, contig, потоки.
# Выход: индексируемый хромосомный gVCF (.g.vcf.gz).
# Не является совместным генотипированием: это вход для GenomicsDBImport.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: haplotypecaller.sh REF.fa INPUT.bam CONTIG OUT.g.vcf.gz THREADS

GATK HaplotypeCaller -ERC GVCF для одного sample × contig.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then usage; exit 0; fi
if (( $# != 5 )); then usage >&2; exit 2; fi

reference_fasta=$1
input_bam=$2
contig=$3
output_gvcf=$4
threads=$5

mkdir -p "$(dirname "$output_gvcf")"

gatk --java-options "-Xmx8g -XX:ParallelGCThreads=${threads}" HaplotypeCaller \
  -R "$reference_fasta" \
  -I "$input_bam" \
  -L "$contig" \
  -ERC GVCF \
  --native-pair-hmm-threads "$threads" \
  -O "$output_gvcf"

test -s "$output_gvcf"
