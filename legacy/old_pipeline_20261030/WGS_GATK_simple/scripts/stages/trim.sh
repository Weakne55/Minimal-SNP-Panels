#!/usr/bin/env bash
set -euo pipefail
[[ $# == 6 ]] || { echo 'Usage: trim.sh R1 R2 OUT_R1 OUT_R2 JSON HTML' >&2; exit 2; }
mkdir -p "$(dirname "$3")" "$(dirname "$4")" "$(dirname "$5")" "$(dirname "$6")"
fastp -i "$1" -I "$2" -o "$3" -O "$4" -j "$5" -h "$6" -w 4 --detect_adapter_for_pe
