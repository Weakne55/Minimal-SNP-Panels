#!/usr/bin/env bash
# One entry point for users. The Python program handles the workflow.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec python3 "$ROOT/scripts/pipeline.py" "$@"
