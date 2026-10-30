#!/usr/bin/env bash
set -euo pipefail
[[ $# == 7 ]] || { echo 'Usage: haplotypecaller.sh REF BAM CONTIG GVCF HEAP_GB THREADS TMPDIR' >&2; exit 2; }
ref=$1 bam=$2 contig=$3 out=$4 heap=$5 threads=$6 tmp=$7
mkdir -p "$(dirname "$out")" "$tmp"
gatk --java-options "-Xmx${heap}g -Djava.io.tmpdir=${tmp}" HaplotypeCaller \
  -R "$ref" -I "$bam" -L "$contig" -ERC GVCF \
  --native-pair-hmm-threads "$threads" -O "$out"
# Most GATK builds create the compressed VCF index automatically.
if [[ ! -s "${out}.tbi" && ! -s "${out}.idx" && ! -s "${out}.csi" ]]; then
  gatk --java-options "-Xmx${heap}g" IndexFeatureFile -I "$out"
fi
[[ -s "${out}.tbi" || -s "${out}.idx" || -s "${out}.csi" ]]
