#!/usr/bin/env bash
set -euo pipefail
[[ $# == 2 ]] || { echo 'Usage: multiqc.sh QC_ROOT OUTDIR' >&2; exit 2; }
mkdir -p "$2"
multiqc --force --outdir "$2" --filename multiqc_report.html "$1"
