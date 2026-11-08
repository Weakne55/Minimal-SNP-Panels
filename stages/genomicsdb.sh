#!/usr/bin/env bash
# Один GenomicsDBImport на вид: все выбранные аутосомы передаются как -L.
# Usage: genomicsdb.sh REF.fa SAMPLE_MAP.tsv DB_DIR AUTOSOME [AUTOSOME ...]
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/resource_lib.bash"
if [[ "${1:-}" == '-h' || "${1:-}" == '--help' ]]; then
    echo 'Usage: genomicsdb.sh REF.fa SAMPLE_MAP.tsv DB_DIR AUTOSOME [AUTOSOME ...]'; exit 0
fi
if (( $# < 4 )); then echo 'Usage: genomicsdb.sh REF.fa SAMPLE_MAP.tsv DB_DIR AUTOSOME [AUTOSOME ...]' >&2; exit 2; fi
reference_fasta=$1
sample_map=$2
db_dir=$3
shift 3
threads=${MIN_SNP_THREADS:-1}
minsnp_check_threads "$threads"
[[ ! -e "$db_dir" ]] || { echo "GenomicsDB уже существует: $db_dir" >&2; exit 2; }
[[ -s "$sample_map" ]] || { echo "Пустая sample-map: $sample_map" >&2; exit 2; }
mkdir -p "$(dirname "$db_dir")"
# Важно: при нескольких -L GenomicsDBImport не поддерживает reader-threads > 1.
# Оставляем свободную RAM для нативной памяти TileDB.
heap_mb=$(minsnp_heap_mb 60 24576)
gc_threads=$((threads < 8 ? threads : 8))
batch_size=$(wc -l < "$sample_map")
(( batch_size > 50 )) && batch_size=50
(( batch_size < 1 )) && batch_size=1
intervals=()
for contig in "$@"; do intervals+=(-L "$contig"); done
# Рабочее место для крупного временного хранилища импорта на общем диске.
import_tmp="$(dirname "$db_dir")/genomicsdb_tmp"
mkdir -p "$import_tmp"
gatk --java-options "-Xmx${heap_mb}m -XX:ActiveProcessorCount=${threads} -XX:ParallelGCThreads=${gc_threads}" \
    GenomicsDBImport -R "$reference_fasta" "${intervals[@]}" \
    --sample-name-map "$sample_map" --validate-sample-name-map true \
    --reader-threads 1 --batch-size "$batch_size" --tmp-dir "$import_tmp" \
    --genomicsdb-workspace-path "$db_dir"
test -d "$db_dir"
printf 'complete\n' > "$db_dir/.completed"
