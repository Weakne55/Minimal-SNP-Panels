#!/usr/bin/env bash
# Назначение: опциональная фильтрация/обрезка paired FASTQ через fastp.
# Usage: trim.sh R1.fastq.gz R2.fastq.gz OUT_R1.gz OUT_R2.gz THREADS
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/resource_lib.bash"
if [[ "${1:-}" == '-h' || "${1:-}" == '--help' ]]; then
    echo 'Usage: trim.sh R1.fastq.gz R2.fastq.gz OUT_R1.gz OUT_R2.gz THREADS'; exit 0
fi
if (( $# != 5 )); then echo 'Usage: trim.sh R1.fastq.gz R2.fastq.gz OUT_R1.gz OUT_R2.gz THREADS' >&2; exit 2; fi
input_r1=$1
input_r2=$2
output_r1=$3
output_r2=$4
threads=$5
minsnp_check_threads "$threads"
mkdir -p "$(dirname "$output_r1")" "$(dirname "$output_r2")"
# fastp имеет ограниченную пользу от очень высокого числа потоков; не запрашиваем лишние.
fastp_threads=$((threads < 16 ? threads : 16))
fastp -i "$input_r1" -I "$input_r2" -o "$output_r1" -O "$output_r2" \
    -w "$fastp_threads" -h "${output_r1%.fastq.gz}.html" -j "${output_r1%.fastq.gz}.json"
