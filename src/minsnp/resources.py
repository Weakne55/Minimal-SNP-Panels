"""Поэтапные SLURM allocations и контроль фактических лимитов."""
from __future__ import annotations

import os
import re

LANE_STAGES = ('trim', 'align')
SAMPLE_STAGES = ('merge_bam', 'markdup', 'bam_qc', 'haplotypecaller')
JOINT_STAGES = ('genomicsdb', 'genotype', 'export_vcf')
STAGE_ORDER = LANE_STAGES + SAMPLE_STAGES + JOINT_STAGES
STAGE_NAMES = frozenset(STAGE_ORDER)


def memory_mib(value: str) -> int:
    match = re.fullmatch(r'([1-9][0-9]*)([MGT]?)', str(value).strip().upper())
    if not match:
        raise ValueError(f'Некорректный объём памяти: {value!r}; пример: 64G')
    return int(match[1]) * {'': 1, 'M': 1, 'G': 1024, 'T': 1024 * 1024}[match[2]]


def _validate_entry(stage: str, resource: dict) -> None:
    if not isinstance(resource, dict):
        raise ValueError(f'Отсутствует slurm.stages.{stage}')
    cpu = resource.get('cpus')
    if type(cpu) is not int or cpu < 1:
        raise ValueError(f'slurm.stages.{stage}.cpus должен быть целым >= 1')
    if memory_mib(resource.get('memory', '')) < 1024:
        raise ValueError(f'slurm.stages.{stage}.memory должен быть не менее 1G')
    value = resource.get('time')
    match = (re.fullmatch(r'(?:(\d+)-)?(\d{1,2}):(\d{2}):(\d{2})', value)
             if isinstance(value, str) else None)
    if match is None:
        raise ValueError(f'slurm.stages.{stage}.time: ожидается D-HH:MM:SS или HH:MM:SS')
    days, hours, minutes, seconds = match.groups()
    if (int(minutes) >= 60 or int(seconds) >= 60 or
            (days is not None and int(hours) >= 24) or
            not (int(days or 0) or int(hours) or int(minutes) or int(seconds))):
        raise ValueError(f'slurm.stages.{stage}.time: некорректное значение {value!r}')


def validate_resources(config: dict) -> None:
    slurm = config.get('slurm')
    if not isinstance(slurm, dict):
        raise ValueError('В pipeline.yaml требуется секция slurm')
    maximum = slurm.get('max_concurrent_jobs')
    if type(maximum) is not int or maximum < 1:
        raise ValueError('slurm.max_concurrent_jobs должен быть целым >= 1')
    stages = slurm.get('stages')
    if not isinstance(stages, dict):
        raise ValueError('В pipeline.yaml требуется slurm.stages')
    for name in STAGE_ORDER:
        _validate_entry(name, stages.get(name))
    level = config.get('logging', {}).get('level', 'INFO')
    if level not in {'DEBUG', 'INFO', 'WARNING', 'ERROR'}:
        raise ValueError('logging.level: DEBUG, INFO, WARNING или ERROR')


def get_parallelism(config: dict) -> int:
    validate_resources(config)
    return config['slurm']['max_concurrent_jobs']


def job_allocation(config: dict, stage: str) -> tuple[int, str]:
    validate_resources(config)
    if stage not in STAGE_NAMES or stage not in config['slurm']['stages']:
        raise ValueError(f'Нет ресурсов для стадии {stage}')
    resource = config['slurm']['stages'][stage]
    return resource['cpus'], str(resource['memory']).upper()


def stage_allocation(config: dict, script: str) -> tuple[int, int]:
    name = script.removesuffix('.sh')
    cpus, memory = job_allocation(config, name)
    return cpus, memory_mib(memory)


def stage_effective(config: dict, script: str) -> tuple[int, int]:
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
