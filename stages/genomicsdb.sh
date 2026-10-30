#!/usr/bin/env bash
# Назначение: импорт хромосомных gVCF в новый GenomicsDB workspace.
# Usage: genomicsdb.sh REF.fa CONTIG SAMPLE_MAP.tsv DB_DIR
# CPU: MIN_SNP_THREADS; память: MIN_SNP_MEMORY_MB; batch-size ограничен 50.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/resource_lib.bash"
if [[ "${1:-}" == '-h' || "${1:-}" == '--help' ]]; then
    echo 'Usage: genomicsdb.sh REF.fa CONTIG SAMPLE_MAP.tsv DB_DIR'; exit 0
fi
if (( $# != 4 )); then echo 'Usage: genomicsdb.sh REF.fa CONTIG SAMPLE_MAP.tsv DB_DIR' >&2; exit 2; fi
reference_fasta=$1
contig=$2
sample_map=$3
db_dir=$4
threads=${MIN_SNP_THREADS:-1}
minsnp_check_threads "$threads"
[[ ! -e "$db_dir" ]] || { echo "GenomicsDB уже существует: $db_dir" >&2; exit 2; }
[[ -s "$sample_map" ]] || { echo "Пустая sample-map: $sample_map" >&2; exit 2; }
mkdir -p "$(dirname "$db_dir")"
heap_mb=$(minsnp_heap_mb 60 24576)
reader_threads=$((threads < 8 ? threads : 8))
gc_threads=$((threads < 8 ? threads : 8))
batch_size=$(wc -l < "$sample_map")
(( batch_size > 50 )) && batch_size=50
(( batch_size < 1 )) && batch_size=1
(( reader_threads > batch_size )) && reader_threads=$batch_size

# Не заполняем Xmx всей SLURM RAM: GenomicsDB расходует native memory.
gatk --java-options "-Xmx${heap_mb}m -XX:ActiveProcessorCount=${threads} -XX:ParallelGCThreads=${gc_threads}" \
    GenomicsDBImport -R "$reference_fasta" -L "$contig" \
    --sample-name-map "$sample_map" --validate-sample-name-map true \
    --reader-threads "$reader_threads" --batch-size "$batch_size" \
    --genomicsdb-workspace-path "$db_dir"
test -d "$db_dir"
printf 'complete\n' > "$db_dir/.completed"
