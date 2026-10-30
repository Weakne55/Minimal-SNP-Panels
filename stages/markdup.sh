#!/usr/bin/env bash
# Назначение: отметить дубликаты на уровне всего образца.
# Вход: coordinate-sorted BAM с mate-тегами от samtools fixmate.
# Выход: BAM, индекс BAI и статистика duplicate marking.
# Важно: запускать после объединения всех runs одного sample.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: markdup.sh INPUT.bam OUT.bam METRICS.txt THREADS

Отметить дубликаты (не удалять reads) программой samtools markdup.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then usage; exit 0; fi
if (( $# != 4 )); then usage >&2; exit 2; fi

input_bam=$1
output_bam=$2
metrics_file=$3
threads=$4

mkdir -p "$(dirname "$output_bam")" "$(dirname "$metrics_file")"

samtools markdup -@ "$threads" -f "$metrics_file" \
  "$input_bam" "${output_bam}.tmp.bam"
mv -- "${output_bam}.tmp.bam" "$output_bam"
samtools index -@ "$threads" "$output_bam"
samtools quickcheck -v "$output_bam"
