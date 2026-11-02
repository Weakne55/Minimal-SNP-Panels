"""Регрессионные тесты SLURM-зависимостей и безопасного повторного запуска."""
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest
from minsnp import slurm_wgs


def test_slurm_array_one_for_all_species(tmp_path, monkeypatch):
    plan = {'samples': [{'species': 'beluga', 'sample': 'B1'},
                        {'species': 'polar_bear', 'sample': 'P1'}]}
    (tmp_path / 'slurm_plan.json').write_text(json.dumps(plan))
    cfg = {'slurm': {'partition': 'main', 'max_parallel_samples': 2,
                     'sample': {'cpus': 16, 'memory': '64G', 'time': '2-00:00:00'},
                     'joint': {'cpus': 8, 'memory': '48G', 'time': '12:00:00'}}}
    requests = []
    def fake(cmd, **kw):
        requests.append(cmd)
        return SimpleNamespace(stdout=f'{700+len(requests)}\n', returncode=0)
    monkeypatch.setattr(slurm_wgs.subprocess, 'run', fake)
    jobs = slurm_wgs.submit_jobs(tmp_path, cfg, sys.executable)
    assert '--array=0-1%2' in requests[0]
    assert len(requests) == 3
    assert all('--dependency=afterok:701' in cmd for cmd in requests[1:])
    assert 'joint_beluga' in jobs and 'joint_polar_bear' in jobs
    with pytest.raises(ValueError):
        slurm_wgs.submit_jobs(tmp_path, cfg, sys.executable)


def test_partial_submit_recorded(tmp_path, monkeypatch):
    (tmp_path / 'slurm_plan.json').write_text(json.dumps({'samples': [{'species': 'beluga', 'sample': 'B1'}]}))
    cfg = {'slurm': {'max_parallel_samples': 2,
                     'sample': {'cpus': 16, 'memory': '64G', 'time': '2-00:00:00'},
                     'joint': {'cpus': 8, 'memory': '48G', 'time': '12:00:00'}}}
    calls = []
    def fake(cmd, **kw):
        calls.append(cmd)
        if len(calls) > 1:
            raise RuntimeError('sbatch failed')
        return SimpleNamespace(stdout='501\n', returncode=0)
    monkeypatch.setattr(slurm_wgs.subprocess, 'run', fake)
    with pytest.raises(RuntimeError):
        slurm_wgs.submit_jobs(tmp_path, cfg, sys.executable)
    assert json.loads((tmp_path / 'slurm_jobs.partial.json').read_text()) == {'samples': '501'}
