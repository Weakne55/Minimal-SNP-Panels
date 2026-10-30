#!/usr/bin/env bash
# Назначение: очистка парных FASTQ с помощью fastp.
# Вход: парные FASTQ.gz, пути для очищенных FASTQ.gz, число потоков.
# Выход: два FASTQ.gz и HTML/JSON-отчёты fastp.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: trim.sh R1.fastq.gz R2.fastq.gz OUT_R1.fastq.gz OUT_R2.fastq.gz THREADS

Очистить парные FASTQ программой fastp.
Отчёты: OUT_R1 без суффикса .fastq.gz + .html / .json.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then usage; exit 0; fi
if (( $# != 5 )); then usage >&2; exit 2; fi

input_r1=$1
input_r2=$2
output_r1=$3
output_r2=$4
threads=$5

mkdir -p "$(dirname "$output_r1")" "$(dirname "$output_r2")"
report_prefix=${output_r1%.fastq.gz}

fastp \
  -i "$input_r1" -I "$input_r2" \
  -o "$output_r1" -O "$output_r2" \
  -w "$threads" \
  -h "${report_prefix}.html" \
  -j "${report_prefix}.json"
