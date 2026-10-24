#!/usr/bin/env bash
set -euo pipefail

BASE=/mnt/tank/scratch/plem/data/Bogdan/Minimal-SNP-Panels
DIR="$BASE/data/full_wgs_20261022/Ursus_maritimus/SRR15170742"

URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/042/SRR15170742/SRR15170742_2.fastq.gz"

FINAL="$DIR/SRR15170742_2.fastq.gz"
PART="$FINAL.part"

EXPECTED_SIZE=17532490796
EXPECTED_MD5=8762c901ad35ba5da30ec145c97b8b74

[[ ! -e "$FINAL" ]] || {
    echo "ERROR: final FASTQ already exists"
    exit 1
}

echo "=== DOWNLOADING R2 ==="

curl --http1.1 \
    --fail --location \
    --retry 12 \
    --retry-all-errors \
    --retry-delay 15 \
    --connect-timeout 30 \
    --continue-at - \
    --output "$PART" \
    "$URL"

echo "=== CHECKING SIZE ==="

ACTUAL_SIZE=$(stat -c %s "$PART")
echo "Expected: $EXPECTED_SIZE"
echo "Actual:   $ACTUAL_SIZE"

if [[ "$ACTUAL_SIZE" != "$EXPECTED_SIZE" ]]; then
    echo "ERROR: size mismatch"
    exit 1
fi

echo "=== CHECKING MD5 ==="

ACTUAL_MD5=$(md5sum "$PART" | awk '{print $1}')
echo "Expected: $EXPECTED_MD5"
echo "Actual:   $ACTUAL_MD5"

if [[ "$ACTUAL_MD5" != "$EXPECTED_MD5" ]]; then
    echo "ERROR: MD5 mismatch; preserving .part"
    exit 1
fi

mv -- "$PART" "$FINAL"

echo "SUCCESS: FASTQ downloaded and verified"
echo "$FINAL"
