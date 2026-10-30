#!/usr/bin/env bash
set -euo pipefail

# Calculate alignment and coverage QC for one BAM file.
#
# Arguments:
#   BAM      - input coordinate-sorted and indexed BAM
#   PREFIX   - prefix for all QC output files
#   THREADS  - number of allocated CPUs (at least 1)

if (( $# != 3 )); then
    echo "Usage: bam_qc.sh BAM PREFIX THREADS" >&2
    exit 2
fi

bam=$1
prefix=$2
threads=$3

# THREADS=1 is valid.
# Reject only zero, negative and non-integer values.
if [[ ! $threads =~ ^[1-9][0-9]*$ ]]; then
    echo "THREADS must be a positive integer" >&2
    exit 2
fi

# SAMtools -@ and mosdepth --threads specify additional
# decompression threads, not the total CPU allocation.
extra_threads=$((threads - 1))

# More than 4 decompression threads usually brings little
# benefit for mosdepth.
mosdepth_threads=$extra_threads
if (( mosdepth_threads > 4 )); then
    mosdepth_threads=4
fi

# Create directory for QC reports.
mkdir -p "$(dirname "$prefix")"

# 1. Basic alignment statistics:
# mapped, properly paired, duplicates, supplementary, etc.
samtools flagstat \
    -@ "$extra_threads" \
    "$bam" > "${prefix}.flagstat.txt"

# 2. Detailed alignment statistics:
# read lengths, base qualities, insert sizes, mapping qualities, etc.
samtools stats \
    -@ "$extra_threads" \
    "$bam" > "${prefix}.stats.txt"

# 3. Genome coverage statistics.
#
# --mapq 20           Ignore alignments with MAPQ < 20.
# -n                  Do not write large per-base coverage output.
# --by 1000000        Calculate coverage in 1 Mb windows.
# --thresholds        Count bases with coverage >= 1, 5, 10, 20.
#
# Input BAM must be coordinate-sorted and indexed.
mosdepth \
    --threads "$mosdepth_threads" \
    --mapq 20 \
    --no-per-base \
    --by 1000000 \
    --thresholds 1,5,10,20 \
    "$prefix" "$bam"
