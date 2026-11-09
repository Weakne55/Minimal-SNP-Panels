#!/usr/bin/env bash
set -euo pipefail

out=$1; threads=$2; shift 2
mkdir -p "$(dirname "$out")"
[[ ! -e "$out" ]] || { echo "Output exists: $out" >&2; exit 2; }

if (( $# == 1 )); then
    ln -s "$(realpath "$1")" "$out"
else
    tmp="${out}.incomplete.bam"
    trap 'rm -f -- "$tmp"' EXIT
    samtools merge -@ "$((threads > 1 ? threads - 1 : 0))" -o "$tmp" "$@"
    samtools quickcheck -v "$tmp"
    mv -- "$tmp" "$out"
fi
