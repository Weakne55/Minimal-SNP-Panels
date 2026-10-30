#!/usr/bin/env bash
# Назначение: объединить BAM разных runs/библиотек одного образца.
# Вход: один или несколько coordinate-sorted BAM с корректными read groups.
# Выход: один BAM (для одного входа — символическая ссылка).
# Следующий шаг: единый markdup для объединённого BAM.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: merge_bam.sh OUT.bam THREADS INPUT_1.bam [INPUT_2.bam ...]

Объединить выравнивания одного образца, сохранив read groups.
Если вход один, создаётся ссылка без копирования большого BAM.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then usage; exit 0; fi
if (( $# < 3 )); then usage >&2; exit 2; fi

output_bam=$1
threads=$2
shift 2
input_bams=("$@")

mkdir -p "$(dirname "$output_bam")"

if (( ${#input_bams[@]} == 1 )); then
  ln -sf "$(realpath "${input_bams[0]}")" "$output_bam"
else
  samtools merge -@ "$threads" -f "${output_bam}.tmp" "${input_bams[@]}"
  mv -- "${output_bam}.tmp" "$output_bam"
fi

samtools quickcheck -v "$output_bam"
