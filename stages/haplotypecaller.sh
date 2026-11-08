#!/usr/bin/env bash
# Один запуск HaplotypeCaller на sample: все выбранные аутосомы, один gVCF.
# Usage: haplotypecaller.sh REF.fa INPUT.bam OUT.g.vcf.gz THREADS CONTIG [CONTIG ...]
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/resource_lib.bash"
if [[ "${1:-}" == '-h' || "${1:-}" == '--help' ]]; then
    echo 'Usage: haplotypecaller.sh REF.fa INPUT.bam OUT.g.vcf.gz THREADS CONTIG [CONTIG ...]'; exit 0
fi
if (( $# < 5 )); then echo 'Expected REF BAM OUTPUT THREADS and at least one CONTIG' >&2; exit 2; fi
reference_fasta=$1
input_bam=$2
output_gvcf=$3
threads=$4
shift 4
minsnp_check_threads "$threads"
heap_mb=$(minsnp_heap_mb 70 16384)
gc_threads=$((threads < 8 ? threads : 8))
pairhmm_threads=$((threads < 8 ? threads : 8))
mkdir -p "$(dirname "$output_gvcf")"
[[ ! -e "$output_gvcf" ]] || { echo "gVCF уже существует: $output_gvcf" >&2; exit 2; }
intervals=()
for contig in "$@"; do intervals+=(-L "$contig"); done
# Один invocation для всех аутосом (не per-chromosome HaplotypeCaller).
gatk --java-options "-Xmx${heap_mb}m -XX:ActiveProcessorCount=${threads} -XX:ParallelGCThreads=${gc_threads}" \
    HaplotypeCaller -R "$reference_fasta" -I "$input_bam" "${intervals[@]}" \
    -ERC GVCF --native-pair-hmm-threads "$pairhmm_threads" -O "$output_gvcf"
test -s "$output_gvcf"
