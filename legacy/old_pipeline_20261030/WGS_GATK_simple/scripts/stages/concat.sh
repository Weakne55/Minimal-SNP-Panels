#!/usr/bin/env bash
set -euo pipefail

# Combine chromosome VCF files into one compressed VCF.
#
# Arguments:
#   OUTPUT_VCF  - final .vcf.gz file
#   INPUT_VCFs  - chromosome VCFs in reference contig order
#
# All input files must contain the same samples in the
# same order. Chromosome regions must not overlap.

if (( $# < 2 )); then
    echo "Usage: concat.sh OUTPUT_VCF INPUT_CHR_VCF [INPUT_CHR_VCF...]" >&2
    exit 2
fi

output_vcf=$1
shift

# Check input files before starting.
for input_vcf in "$@"; do
    if [[ ! -s "$input_vcf" ]]; then
        echo "Missing or empty VCF: $input_vcf" >&2
        exit 1
    fi

    if [[ "$input_vcf" == "$output_vcf" ]]; then
        echo "Input and output VCF must differ" >&2
        exit 2
    fi
done

mkdir -p "$(dirname "$output_vcf")"

if (( $# == 1 )); then
    # One input: normalize output compression to BGZF.
    bcftools view \
        -Oz \
        -o "$output_vcf" \
        "$1"
else
    # Multiple chromosome VCFs: concatenate in supplied order.
    bcftools concat \
        -Oz \
        -o "$output_vcf" \
        "$@"
fi

# Create a CSI index for random access.
# -f: overwrite an existing index.
# -c: use CSI rather than TBI.
bcftools index \
    -f \
    -c \
    "$output_vcf"
