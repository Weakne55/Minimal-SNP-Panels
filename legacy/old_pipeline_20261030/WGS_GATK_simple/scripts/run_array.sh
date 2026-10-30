#!/usr/bin/env bash
# SLURM array worker: map SLURM_ARRAY_TASK_ID to one planned task.
set -euo pipefail

if [[ $# != 2 ]]; then
  echo "Usage: run_array.sh PLAN_JSON TASK_IDS_FILE" >&2
  exit 2
fi

: "${SLURM_ARRAY_TASK_ID:?Must be executed inside a SLURM array}"
plan=$1
task_list=$2
line_number=$((SLURM_ARRAY_TASK_ID + 1))
task_id=$(sed -n "${line_number}p" "$task_list")

if [[ -z "$task_id" ]]; then
  echo "No task for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec "${PIPELINE_PYTHON:-python3}" "$root/scripts/pipeline.py" \
  run-task --plan "$plan" --task-id "$task_id"
