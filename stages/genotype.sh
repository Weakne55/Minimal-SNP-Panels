#!/usr/bin/env bash
# Назначение: joint genotyping одного contig из GenomicsDB.
# Usage: genotype.sh REF.fa GENOMICS_DB_DIR OUT.vcf.gz
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/resource_lib.bash"
if [[ "${1:-}" == '-h' || "${1:-}" == '--help' ]]; then
    echo 'Usage: genotype.sh REF.fa GENOMICS_DB_DIR OUT.vcf.gz'; exit 0
fi
if (( $# != 3 )); then echo 'Usage: genotype.sh REF.fa GENOMICS_DB_DIR OUT.vcf.gz' >&2; exit 2; fi
reference_fasta=$1
db_dir=$2
output_vcf=$3
threads=${MIN_SNP_THREADS:-1}
minsnp_check_threads "$threads"
heap_mb=$(minsnp_heap_mb 75 16384)
gc_threads=$((threads < 8 ? threads : 8))
mkdir -p "$(dirname "$output_vcf")"
gatk --java-options "-Xmx${heap_mb}m -XX:ActiveProcessorCount=${threads} -XX:ParallelGCThreads=${gc_threads}" \
    GenotypeGVCFs -R "$reference_fasta" -V "gendb://${db_dir}" -O "$output_vcf"
test -s "$output_vcf"
