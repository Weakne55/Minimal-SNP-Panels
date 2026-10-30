#!/usr/bin/env bash
# Назначение: одна пара FASTQ -> name-collate -> fixmate -> coordinate-sorted BAM.
# Usage: align.sh REF.fa R1.fastq.gz R2.fastq.gz OUT.bam SAMPLE LIBRARY RUN THREADS
# Память: MIN_SNP_MEMORY_MB; sort ограничивает память на поток через -m.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/resource_lib.bash"

usage() {
    echo 'Usage: align.sh REF.fa R1.fastq.gz R2.fastq.gz OUT.bam SAMPLE LIBRARY RUN THREADS'
}
if [[ "${1:-}" == '-h' || "${1:-}" == '--help' ]]; then usage; exit 0; fi
if (( $# != 8 )); then usage >&2; exit 2; fi

reference_fasta=$1
input_r1=$2
input_r2=$3
output_bam=$4
sample_id=$5
library_id=$6
run_id=$7
threads=$8
minsnp_check_threads "$threads"
for input in "$reference_fasta" "$input_r1" "$input_r2"; do
    [[ -s "$input" ]] || { echo "Нет входного файла: $input" >&2; exit 2; }
done
for extension in .amb .ann .pac .bwt.2bit.64; do
    [[ -s "${reference_fasta}${extension}" ]] || {
        echo "Нет индекса BWA-MEM3: ${reference_fasta}${extension}" >&2; exit 2;
    }
done
[[ ! -e "$output_bam" ]] || { echo "BAM уже существует: $output_bam" >&2; exit 2; }
mkdir -p "$(dirname "$output_bam")"

read_group="@RG\\tID:${sample_id}.${library_id}.${run_id}\\tSM:${sample_id}\\tLB:${library_id}\\tPL:ILLUMINA\\tPU:${library_id}.${run_id}"
memory_mb=$(minsnp_memory_mb)

# В pipe одновременно работают BWA, collate, fixmate, sort.
# Малый CPU-бюджет допускает конкурентные однопоточные процессы (time sharing).
# Дополнительные потоки sort выделяем только при CPU >= 5.
sort_workers=$((threads >= 5 ? threads / 4 : 0))
# RAM ограничивает число потоков сортировки при очень большом CPU-запросе.
max_sort_workers=$((memory_mb / 1024))
(( max_sort_workers > 1 )) && max_sort_workers=$((max_sort_workers - 1))
(( sort_workers > max_sort_workers )) && sort_workers=$max_sort_workers
bwa_threads=$((threads > sort_workers + 3 ? threads - sort_workers - 3 : 1))
sort_memory_mb=$((memory_mb / 6 / (sort_workers + 1)))
(( sort_memory_mb > 768 )) && sort_memory_mb=768
(( sort_memory_mb < 64 )) && sort_memory_mb=64

scratch_root=${SLURM_TMPDIR:-${TMPDIR:-$(dirname "$output_bam")}}
[[ -d "$scratch_root" ]] || { echo "Нет временного каталога: $scratch_root" >&2; exit 2; }
workdir=$(mktemp -d "${scratch_root%/}/minsnp-align.XXXXXXXX")
tmp_bam="${output_bam}.incomplete.bam"
trap 'rm -rf -- "$workdir"; rm -f -- "$tmp_bam"' EXIT

printf 'align: CPU=%s BWA=%s samtools-sort-workers=%s sort-m=%sM RAM-limit=%sMiB\n' \
    "$threads" "$bwa_threads" "$sort_workers" "$sort_memory_mb" "$memory_mb" >&2
bwa-mem3 mem --compat=bwa-mem -t "$bwa_threads" -R "$read_group" \
    "$reference_fasta" "$input_r1" "$input_r2" |
    samtools collate -O -u - "$workdir/collate" |
    samtools fixmate -m -u - - |
    samtools sort -@ "$sort_workers" -m "${sort_memory_mb}M" \
        -T "$workdir/sort" -o "$tmp_bam" -
samtools quickcheck -v "$tmp_bam"
mv -- "$tmp_bam" "$output_bam"
trap - EXIT
rm -rf -- "$workdir"
