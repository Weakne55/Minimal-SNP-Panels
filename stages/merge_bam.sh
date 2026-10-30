#!/usr/bin/env bash
# Назначение: объединить BAM нескольких sequencing runs одного образца.
# Usage: merge_bam.sh OUT.bam THREADS IN1.bam [IN2.bam ...]
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/resource_lib.bash"
if [[ "${1:-}" == '-h' || "${1:-}" == '--help' ]]; then
    echo 'Usage: merge_bam.sh OUT.bam THREADS IN1.bam [IN2.bam ...]'; exit 0
fi
if (( $# < 3 )); then echo 'Usage: merge_bam.sh OUT.bam THREADS IN1.bam [IN2.bam ...]' >&2; exit 2; fi
output_bam=$1
threads=$2
shift 2
minsnp_check_threads "$threads"
[[ ! -e "$output_bam" ]] || { echo "BAM уже существует: $output_bam" >&2; exit 2; }
mkdir -p "$(dirname "$output_bam")"
for input_bam in "$@"; do
    samtools quickcheck -v "$input_bam"
done
if (( $# == 1 )); then
    ln -s "$(realpath "$1")" "$output_bam"
else
    tmp_bam="${output_bam}.incomplete.bam"
    trap 'rm -f -- "$tmp_bam"' EXIT
    samtools merge -@ "$(minsnp_samtools_workers "$threads")" -o "$tmp_bam" "$@"
    samtools quickcheck -v "$tmp_bam"
    mv -- "$tmp_bam" "$output_bam"
    trap - EXIT
fi
samtools quickcheck -v "$output_bam"
