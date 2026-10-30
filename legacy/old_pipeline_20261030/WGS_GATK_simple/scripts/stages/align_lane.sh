#!/usr/bin/env bash
set -euo pipefail

# Align one paired-end FASTQ lane to the reference.
# Output: coordinate-sorted BAM with mate tags for markdup.

if (( $# != 8 )); then
    echo "Usage: align_lane.sh REF R1 R2 OUT_BAM SAMPLE LIBRARY LANE THREADS" >&2
    exit 2
fi

# Input arguments
reference=$1
read1=$2
read2=$3
output_bam=$4
sample=$5
library=$6
lane=$7
threads=$8

# Accept any positive number of threads, including 1.
if [[ ! $threads =~ ^[1-9][0-9]*$ ]]; then
    echo "THREADS must be a positive integer" >&2
    exit 2
fi

# samtools sort -@ specifies ADDITIONAL worker threads.
# With 1 CPU, all commands use their single-threaded mode.
sort_extra_threads=0

if (( threads >= 7 )); then
    sort_extra_threads=2
elif (( threads >= 4 )); then
    sort_extra_threads=1
fi

bwa_threads=$((threads - sort_extra_threads))

# Create output directory.
mkdir -p "$(dirname "$output_bam")"

# Prefer local scratch space on the SLURM compute node.
scratch_root=${SLURM_TMPDIR:-${TMPDIR:-$(dirname "$output_bam")}}
mkdir -p "$scratch_root"

workdir=$(mktemp -d "$scratch_root/wgs-align.XXXXXXXX")

# Remove temporary files when the script exits.
trap 'rm -rf -- "$workdir"' EXIT

# Read-group metadata used to identify sample, library and lane.
#
# ID = read-group identifier
# SM = sample name
# LB = library
# PL = sequencing platform
# PU = platform unit
#
# BWA expects literal \t separators.
printf -v read_group '@RG\\tID:%s.%s.%s\\tSM:%s\\tLB:%s\\tPL:ILLUMINA\\tPU:%s.%s' \
    "$sample" "$library" "$lane" \
    "$sample" "$library" "$library" "$lane"

# Alignment and BAM preparation:
#
# 1. bwa-mem3: align reads to the reference.
# 2. collate: group paired reads together.
# 3. fixmate: add mate information and tags required by markdup.
# 4. sort: sort BAM by genomic coordinates.
#
# Intermediate data are streamed through pipes.

bwa-mem3 mem \
    --compat=bwa-mem \
    -t "$bwa_threads" \
    -R "$read_group" \
    "$reference" "$read1" "$read2" \
  | samtools collate -O -u - "$workdir/collate" \
  | samtools fixmate -m -u - - \
  | samtools sort \
        -@ "$sort_extra_threads" \
        -m 1G \
        -T "$workdir/sort" \
        -o "$output_bam" -
