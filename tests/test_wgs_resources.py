"""Проверка единого бюджета WGS-образца и joint-задания."""
import pytest
from minsnp.resources import memory_mib, validate_resources, stage_effective


def config():
    return {'slurm': {'max_parallel_samples': 3,
                      'sample': {'cpus': 16, 'memory': '64G', 'time': '2-00:00:00'},
                      'joint': {'cpus': 8, 'memory': '48G', 'time': '12:00:00'}}}


def test_memory_units():
    assert memory_mib('64G') == 65536
    assert memory_mib('49152M') == 49152
    with pytest.raises(ValueError):
        memory_mib('64GB')


def test_stage_uses_group_allocation(monkeypatch):
    for name in ('SLURM_CPUS_PER_TASK', 'SLURM_MEM_PER_NODE', 'SLURM_MEM_PER_CPU'):
        monkeypatch.delenv(name, raising=False)
    cfg = config()
    assert stage_effective(cfg, 'align.sh') == (16, 65536)
    assert stage_effective(cfg, 'haplotypecaller.sh') == (16, 65536)
    assert stage_effective(cfg, 'genomicsdb.sh') == (8, 49152)


def test_slurm_allocation_too_small_fails(monkeypatch):
    monkeypatch.setenv('SLURM_CPUS_PER_TASK', '4')
    monkeypatch.setenv('SLURM_MEM_PER_NODE', '8192')
    with pytest.raises(RuntimeError):
        stage_effective(config(), 'align.sh')


def test_bad_config_rejected():
    cfg = config()
    cfg['slurm']['joint']['memory'] = '-4G'
    with pytest.raises(ValueError):
        validate_resources(cfg)
