#!/usr/bin/env bash
# Назначение: выравнивание одной FASTQ-пары на референс с BWA-MEM3.
# Вход: FASTA, R1/R2, sample/library/run IDs, число CPU.
# Выход: координатно отсортированный BAM с RG и тегами fixmate.
# Замечание: markdup выполняется ПОСЛЕ объединения BAM одного образца.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: align.sh REF.fa R1.fastq.gz R2.fastq.gz OUT.bam SAMPLE LIBRARY RUN THREADS

Выровнять одну FASTQ-пару: BWA-MEM3 → collate → fixmate → sort.
SAMPLE, LIBRARY, RUN используются для записи read group (@RG).
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then usage; exit 0; fi
if (( $# != 8 )); then usage >&2; exit 2; fi

reference_fasta=$1
input_r1=$2
input_r2=$3
output_bam=$4
sample_id=$5
library_id=$6
run_id=$7
threads=$8

mkdir -p "$(dirname "$output_bam")"

# Read groups сохраняют происхождение reads в объединённом BAM.
read_group="@RG\\tID:${sample_id}.${library_id}.${run_id}\\tSM:${sample_id}\\tLB:${library_id}\\tPL:ILLUMINA\\tPU:${library_id}.${run_id}"

# Временные файлы размещаются в локальном scratch, если он задан.
scratch_root=${SLURM_TMPDIR:-${TMPDIR:-$(dirname "$output_bam")}}
workdir=$(mktemp -d "${scratch_root}/minsnp-align.XXXXXXXX")
trap 'rm -rf -- "$workdir"' EXIT

# Сохраняем прежнее простое распределение CPU, не оптимизируя SLURM.
bwa_threads=$((threads > 2 ? threads - 2 : 1))
sort_threads=$((threads > 2 ? 2 : 0))

# fixmate требует name-collated BAM; sort возвращает порядок по координатам.
bwa-mem3 mem --compat=bwa-mem \
  -t "$bwa_threads" \
  -R "$read_group" \
  "$reference_fasta" "$input_r1" "$input_r2" |
  samtools collate -O -u - "$workdir/collate" |
  samtools fixmate -m -u - - |
  samtools sort -@ "$sort_threads" -T "$workdir/sort" \
    -o "${output_bam}.tmp" -

mv -- "${output_bam}.tmp" "$output_bam"
samtools quickcheck -v "$output_bam"
