#!/usr/bin/env bash
set -euo pipefail
[[ $# == 8 ]] || { echo 'Usage: align_lane.sh REF R1 R2 OUT_BAM SAMPLE LIBRARY LANE THREADS' >&2; exit 2; }
ref=$1 r1=$2 r2=$3 out=$4 sample=$5 lib=$6 lane=$7 threads=$8
(( threads >= 6 )) || { echo 'Alignment needs at least 6 CPUs' >&2; exit 2; }
mkdir -p "$(dirname "$out")"
tmp_root="${SLURM_TMPDIR:-${TMPDIR:-$(dirname "$out")}}"
mkdir -p "$tmp_root"
tmpdir="$(mktemp -d "$tmp_root/wgs-align.XXXXXXXX")"
trap 'rm -rf -- "$tmpdir"' EXIT
printf -v rg '@RG\tID:%s.%s.%s\tSM:%s\tLB:%s\tPL:ILLUMINA\tPU:%s.%s' \
  "$sample" "$lib" "$lane" "$sample" "$lib" "$lib" "$lane"
# Use literal backslash-t in BWA -R (BWA decodes it to tab).
rg=${rg//$'\t'/\\t}
bwa-mem3 mem --compat=bwa-mem -t "$((threads-4))" -R "$rg" "$ref" "$r1" "$r2" \
  | samtools collate -O -u - "$tmpdir/collate" \
  | samtools fixmate -m -u - - \
  | samtools sort -@ 2 -m 1G -T "$tmpdir/sort" -o "$out" -
