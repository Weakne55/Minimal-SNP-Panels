"""Проверка и предоставление ресурсов для WGS и joint genotyping."""
from __future__ import annotations

import os
import re

WGS_STAGES = frozenset({'trim', 'align', 'merge_bam', 'markdup', 'bam_qc', 'haplotypecaller'})
JOINT_STAGES = frozenset({'genomicsdb', 'genotype', 'concat', 'export_vcf'})


def memory_mib(value: str) -> int:
    """Перевести значение памяти из YAML в MiB (используется двоичная шкала)."""
    match = re.fullmatch(r'([1-9][0-9]*)([MGT]?)', str(value).strip().upper())
    if not match:
        raise ValueError(f'Некорректный объём памяти: {value!r}; пример: 64G')
    return int(match[1]) * {'': 1, 'M': 1, 'G': 1024, 'T': 1024 * 1024}[match[2]]


def validate_resources(config: dict) -> None:
    """Проверить границы конфигурации один раз перед созданием заданий."""
    slurm = config.get('slurm')
    if not isinstance(slurm, dict):
        raise ValueError('В config/pipeline.yaml требуется секция slurm')
    max_parallel = slurm.get('max_parallel_samples')
    if type(max_parallel) is not int or max_parallel < 1:
        raise ValueError('slurm.max_parallel_samples должен быть целым числом >= 1')
    for group in ('sample', 'joint'):
        resource = slurm.get(group)
        if not isinstance(resource, dict):
            raise ValueError(f'Отсутствует slurm.{group}')
        cpu = resource.get('cpus')
        if type(cpu) is not int or cpu < 1:
            raise ValueError(f'slurm.{group}.cpus должен быть целым числом >= 1')
        if memory_mib(resource.get('memory', '')) < 1024:
            raise ValueError(f'slurm.{group}.memory должен быть не менее 1G')
        value = resource.get('time')
        match = (re.fullmatch(r'(?:(\d+)-)?(\d{1,2}):(\d{2}):(\d{2})', value)
                 if isinstance(value, str) else None)
        if match is None:
            raise ValueError(f'slurm.{group}.time: ожидается D-HH:MM:SS или HH:MM:SS')
        days, hours, minutes, seconds = match.groups()
        # При наличии дней часы ограничены 0–23; без дней допускаются 00–99 часов.
        if (int(minutes) >= 60 or int(seconds) >= 60 or
                (days is not None and int(hours) >= 24) or
                not (int(days or 0) or int(hours) or int(minutes) or int(seconds))):
            raise ValueError(f'slurm.{group}.time: некорректные часы, минуты или секунды: {value!r}')


def get_parallelism(config: dict) -> int:
    validate_resources(config)
    return config['slurm']['max_parallel_samples']


def get_stage_resources(config: dict, stage: str) -> tuple[int, int]:
    """Получить единый бюджет задания для любой его внутренней Bash-стадии."""
    name = stage.removesuffix('.sh')
    if name in WGS_STAGES:
        group = 'sample'
    elif name in JOINT_STAGES:
        group = 'joint'
    else:
        raise ValueError(f'Неизвестная стадия: {stage}')
    resource = config['slurm'][group]
    return resource['cpus'], memory_mib(resource['memory'])


def stage_allocation(config: dict, script: str) -> tuple[int, int]:
    """Совместимый интерфейс для существующих модулей WGS и core."""
    validate_resources(config)
    return get_stage_resources(config, script)


def job_allocation(config: dict, group: str) -> tuple[int, str]:
    validate_resources(config)
    if group not in {'sample', 'joint'}:
        raise ValueError(f'Неизвестная группа ресурсов: {group}')
    resource = config['slurm'][group]
    return resource['cpus'], str(resource['memory']).upper()


def stage_effective(config: dict, script: str) -> tuple[int, int]:
    """Запретить тихое превышение фактической аллокации SLURM."""
    threads, memory = stage_allocation(config, script)
    cpu_limit = int(os.environ.get('SLURM_CPUS_PER_TASK') or 0)
    memory_limit = int(os.environ.get('SLURM_MEM_PER_NODE') or 0)
    if not memory_limit and os.environ.get('SLURM_MEM_PER_CPU') and cpu_limit:
        memory_limit = int(os.environ['SLURM_MEM_PER_CPU']) * cpu_limit
    if cpu_limit and threads > cpu_limit:
        raise RuntimeError(f'{script}: запрошено {threads} CPU, SLURM выделил {cpu_limit}')
    if memory_limit and memory > memory_limit:
        raise RuntimeError(f'{script}: запрошено {memory} MiB, SLURM выделил {memory_limit} MiB')
    return threads, memory
