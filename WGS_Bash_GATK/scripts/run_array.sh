#!/usr/bin/env bash
set -euo pipefail
[[ $# == 2 ]] || { echo 'Usage: run_array.sh PLAN_JSON IDS_FILE' >&2; exit 2; }
: "${SLURM_ARRAY_TASK_ID:?Needs SLURM_ARRAY_TASK_ID}"
task_id="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$2")"
[[ -n "$task_id" ]] || { echo 'Array index missing in IDs list' >&2; exit 1; }
py="${PIPELINE_PYTHON:-python3}"
pipeline_py="$("$py" -c 'import json, pathlib, sys; print(pathlib.Path(json.load(open(sys.argv[1]))["root"]) / "scripts" / "pipeline.py")' "$1")"
[[ -f "$pipeline_py" ]] || { echo "ERROR: pipeline.py not found: $pipeline_py" >&2; exit 2; }
exec "$py" "$pipeline_py" run-task --plan "$1" --task-id "$task_id"
