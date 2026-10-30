#!/usr/bin/env bash
set -euo pipefail
[[ $# == 8 ]] || { echo 'Usage: genomicsdb.sh REF CONTIG SAMPLE_MAP DB_DIR HEAP_GB BATCH_SIZE READERS TMPDIR' >&2; exit 2; }
ref=$1 contig=$2 samplemap=$3 db=$4 heap=$5 batch=$6 readers=$7 tmp=$8
[[ $db == */genomicsdb/* ]] || { echo "Unsafe GenomicsDB destination: $db" >&2; exit 2; }
[[ $db != '/' && -n $db ]] || exit 2
mkdir -p "$(dirname "$db")" "$tmp"
# Only the pipeline-owned single-contig workspace is replaced during a rerun.
if [[ -e $db ]]; then rm -rf -- "$db"; fi
gatk --java-options "-Xmx${heap}g -Djava.io.tmpdir=${tmp}" GenomicsDBImport \
  -R "$ref" -L "$contig" --sample-name-map "$samplemap" \
  --genomicsdb-workspace-path "$db" --batch-size "$batch" \
  --reader-threads "$readers" --validate-sample-name-map true --tmp-dir "$tmp"
[[ -d $db ]] || { echo 'GenomicsDBImport did not create workspace' >&2; exit 1; }
printf "completed\n" > "$db/.completed"
