#!/usr/bin/env bash
# Назначение: HaplotypeCaller (-ERC GVCF) для одного sample × contig.
# Usage: haplotypecaller.sh REF.fa INPUT.bam CONTIG OUT.g.vcf.gz THREADS
# Java heap зависит от RAM; PairHMM threads от доступных CPU.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/resource_lib.bash"
if [[ "${1:-}" == '-h' || "${1:-}" == '--help' ]]; then
    echo 'Usage: haplotypecaller.sh REF.fa INPUT.bam CONTIG OUT.g.vcf.gz THREADS'; exit 0
fi
if (( $# != 5 )); then echo 'Usage: haplotypecaller.sh REF.fa INPUT.bam CONTIG OUT.g.vcf.gz THREADS' >&2; exit 2; fi
reference_fasta=$1
input_bam=$2
contig=$3
output_gvcf=$4
threads=$5
minsnp_check_threads "$threads"
heap_mb=$(minsnp_heap_mb 70 16384)
gc_threads=$((threads < 8 ? threads : 8))
pairhmm_threads=$((threads < 8 ? threads : 8))
mkdir -p "$(dirname "$output_gvcf")"
gatk --java-options "-Xmx${heap_mb}m -XX:ActiveProcessorCount=${threads} -XX:ParallelGCThreads=${gc_threads}" \
    HaplotypeCaller -R "$reference_fasta" -I "$input_bam" -L "$contig" \
    -ERC GVCF --native-pair-hmm-threads "$pairhmm_threads" -O "$output_gvcf"
test -s "$output_gvcf"
