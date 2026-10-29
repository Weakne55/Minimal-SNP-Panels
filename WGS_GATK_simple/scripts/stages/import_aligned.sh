#!/usr/bin/env bash
set -euo pipefail
[[ $# == 5 ]] || { echo 'Usage: import_aligned.sh MODE INPUT OUTPUT_OR_REPORT SAMPLE REF' >&2; exit 2; }
mode=$1 src=$2 out=$3 sample=$4 ref=$5
python3 "$(dirname "$0")/../validate_alignment.py" "$src" "$sample" "$ref"
samtools quickcheck -v "$src"
if [[ $mode == bam ]]; then
  # Caller uses the original sorted/indexed BAM in place. Produce small report only.
  samtools idxstats "$src" > "$out"
elif [[ $mode == cram ]]; then
  # Conservative CRAM -> BAM conversion; no duplicate remarking.
  samtools view -@ 4 -T "$ref" -b -o "$out" "$src"
  samtools index -@ 2 "$out"
else
  echo "Unknown aligned input: $mode" >&2; exit 2
fi
