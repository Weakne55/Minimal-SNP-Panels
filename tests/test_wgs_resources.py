"""Проверки расчёта CPU/RAM и формирования SLURM job без реальных FASTQ."""
from __future__ import annotations

from types import SimpleNamespace
import subprocess

import pytest

from minsnp.resources import configure_resources, memory_mib
from minsnp import cli


@pytest.mark.parametrize("threads", [1, 2, 8, 64, 128])
def test_any_positive_cpu_budget(monkeypatch, threads):
    """Вычисления принимают любой положительный CPU-бюджет."""
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    monkeypatch.delenv("SLURM_MEM_PER_NODE", raising=False)
    monkeypatch.delenv("SLURM_MEM_PER_CPU", raising=False)
    args = SimpleNamespace(threads=threads, memory="64G")
    cfg = configure_resources(args, {"threads": 8, "slurm": {"memory": "64G"}})
    assert cfg["threads"] == threads
    assert cfg["_memory_mb"] == 65536
    assert cfg["_requested_threads"] == threads


def test_cpu_and_memory_are_capped_by_slurm(monkeypatch):
    """Реальное ограничение SLURM приоритетнее пользовательских настроек."""
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "1")
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "16384")
    cfg = configure_resources(SimpleNamespace(threads=16, memory="64G"),
                              {"threads": 8, "slurm": {"memory": "64G"}})
    assert cfg["threads"] == 1
    assert cfg["_memory_mb"] == 16384


def test_cpu_defaults_to_slurm_allocation(monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "64")
    monkeypatch.delenv("SLURM_MEM_PER_NODE", raising=False)
    monkeypatch.delenv("SLURM_MEM_PER_CPU", raising=False)
    cfg = configure_resources(SimpleNamespace(threads=None, memory="64G"),
                              {"threads": 8, "slurm": {"memory": "64G"}})
    assert cfg["threads"] == 64


@pytest.mark.parametrize("bad", [0, -1, "four", None])
def test_bad_threads_rejected(monkeypatch, bad):
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    with pytest.raises(ValueError):
        configure_resources(SimpleNamespace(threads=bad, memory="64G"),
                            {"threads": 0, "slurm": {"memory": "64G"}})


def test_memory_units():
    assert memory_mib("1G") == 1024
    assert memory_mib("64G") == 65536
    assert memory_mib("65536M") == 65536
    with pytest.raises(ValueError):
        memory_mib("64GB")


def test_submission_requests_cpu_and_ram(monkeypatch, tmp_path):
    """Проверка sbatch без запуска SLURM и биологических программ."""
    commands = []

    def fake_run(command, **kwargs):
        commands.append(command)
        return subprocess.CompletedProcess(command, 0, stdout="12345\n", stderr="")

    monkeypatch.setattr(cli, "require_commands", lambda *names: None)
    monkeypatch.setattr(cli.subprocess, "run", fake_run)
    monkeypatch.setattr(cli.sys, "argv", ["cli.py", "wgs-vcf", "--species", "beluga", "--backend", "slurm"])
    args = SimpleNamespace(species="beluga", threads=None, memory=None)
    cfg = {"slurm": {"time": "1-00:00:00", "partition": ""},
           "_requested_threads": 1, "_requested_memory": "32G"}
    cli.submit_slurm(args, cfg, tmp_path)
    sbatch = next(cmd for cmd in commands if cmd[0] == "sbatch")
    assert sbatch[sbatch.index("--cpus-per-task") + 1] == "1"
    assert sbatch[sbatch.index("--mem") + 1] == "32G"
    script = (tmp_path / "slurm_job.sh").read_text()
    assert "--threads 1" in script
    assert "--memory 32G" in script
    assert "--resume" in script
