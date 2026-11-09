#!/usr/bin/env bash
set -euo pipefail

ref=$1; r1=$2; r2=$3; out=$4
sample=$5; library=$6; run=$7; threads=$8
mkdir -p "$(dirname "$out")"
[[ ! -e "$out" ]] || { echo "Output exists: $out" >&2; exit 2; }

# Часть CPU оставляем samtools; остальные выделяем BWA-MEM3.
sort_threads=$(( threads >= 5 ? threads / 4 : 0 ))
bwa_threads=$(( threads > sort_threads + 3 ? threads - sort_threads - 3 : 1 ))
sort_mem=$(( ${MIN_SNP_MEMORY_MB:?} / (6 * (sort_threads + 1)) ))
sort_mem=$(( sort_mem < 64 ? 64 : sort_mem > 768 ? 768 : sort_mem ))

rg="@RG\\tID:${sample}.${library}.${run}\\tSM:${sample}\\tLB:${library}\\tPL:ILLUMINA\\tPU:${library}.${run}"

workdir=$(mktemp -d "${SLURM_TMPDIR:-${TMPDIR:-$(dirname "$out")}}/minsnp-align.XXXXXXXX")
tmp="${out}.incomplete.bam"
trap 'rm -rf -- "$workdir"; rm -f -- "$tmp"' EXIT

bwa-mem3 mem --compat=bwa-mem -t "$bwa_threads" -R "$rg" "$ref" "$r1" "$r2" |
    samtools collate -O -u - "$workdir/collate" |
    samtools fixmate -m -u - - |
    samtools sort -@ "$sort_threads" -m "${sort_mem}M" -T "$workdir/sort" -o "$tmp" -

samtools quickcheck -v "$tmp"
mv -- "$tmp" "$out"
