"""Test SLURM array dispatch without SLURM or genomics commands."""

import os
import subprocess
from pathlib import Path

RUN = Path(__file__).resolve().parents[1] / "scripts" / "run_array.sh"


def test_array_worker_dispatch(tmp_path):
    tasks = tmp_path / "tasks.txt"
    tasks.write_text("align.A\nalign.B\n")
    env = {**os.environ, "SLURM_ARRAY_TASK_ID": "1", "PIPELINE_PYTHON": "/bin/echo"}
    result = subprocess.run(["bash", str(RUN), "plan.json", str(tasks)],
                            env=env, capture_output=True, text=True, check=True)
    assert result.stdout.strip() == (
        f"{RUN.parents[1]}/scripts/pipeline.py run-task --plan plan.json --task-id align.B"
    )


def test_array_worker_rejects_invalid_job(tmp_path):
    tasks = tmp_path / "tasks.txt"
    tasks.write_text("align.A\n")
    env = {**os.environ, "SLURM_ARRAY_TASK_ID": "2", "PIPELINE_PYTHON": "/bin/echo"}
    result = subprocess.run(["bash", str(RUN), "plan.json", str(tasks)],
                            env=env, capture_output=True, text=True)
    assert result.returncode == 1
    assert "No task for SLURM_ARRAY_TASK_ID=2" in result.stderr
