#!/usr/bin/env bash
set -euo pipefail
(( $# >= 5 )) || { echo 'Usage: markdup.sh OUT_BAM JSON_METRICS OPTICAL_DISTANCE THREADS LANE_BAM...' >&2; exit 2; }
out=$1 metrics=$2 optical=$3 threads=$4; shift 4
mkdir -p "$(dirname "$out")" "$(dirname "$metrics")"
# A lane BAM already has MC/ms tags and coordinate sorting.
if (( $# == 1 )); then
  samtools markdup -@ "$((threads-2))" --json -s -d "$optical" -f "$metrics" "$1" "$out"
else
  samtools merge -u -@ 2 - "$@" \
    | samtools markdup -@ "$((threads-4))" --json -s -d "$optical" -f "$metrics" - "$out"
fi
samtools index -@ 2 "$out"
