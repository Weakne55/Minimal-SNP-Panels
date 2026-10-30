#!/usr/bin/env bash
set -euo pipefail
[[ $# == 6 ]] || { echo 'Usage: genotype.sh REF CONTIG DB_DIR OUT_VCF HEAP_GB TMPDIR' >&2; exit 2; }
ref=$1 contig=$2 db=$3 out=$4 heap=$5 tmp=$6
[[ -s "$db/.completed" ]] || { echo "Incomplete GenomicsDB: $db" >&2; exit 1; }
mkdir -p "$(dirname "$out")" "$tmp"
gatk --java-options "-Xmx${heap}g -Djava.io.tmpdir=${tmp}" GenotypeGVCFs \
  -R "$ref" -L "$contig" -V "gendb://$db" -O "$out" --tmp-dir "$tmp"
if [[ ! -s "${out}.tbi" && ! -s "${out}.idx" && ! -s "${out}.csi" ]]; then
  gatk --java-options "-Xmx${heap}g" IndexFeatureFile -I "$out"
fi
[[ -s "${out}.tbi" || -s "${out}.idx" || -s "${out}.csi" ]]
