#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/resource_lib.bash"

ref=$1; map=$2; db=$3; shift 3
threads=${MIN_SNP_THREADS:?}
heap=$(minsnp_heap_mb 60 24576)
batch=$(wc -l < "$map")
batch=$(( batch > 50 ? 50 : batch ))

intervals=()
for contig in "$@"; do intervals+=(-L "$contig"); done
workdir="$(dirname "$db")/genomicsdb_tmp"
mkdir -p "$workdir"
[[ ! -e "$db" ]] || { echo "Database exists: $db" >&2; exit 2; }

gatk --java-options "-Xmx${heap}m -XX:ActiveProcessorCount=${threads}" \
    GenomicsDBImport -R "$ref" "${intervals[@]}" \
    --sample-name-map "$map" --validate-sample-name-map true \
    --reader-threads 1 --batch-size "$batch" --tmp-dir "$workdir" \
    --genomicsdb-workspace-path "$db"

printf 'complete\n' > "$db/.completed"
