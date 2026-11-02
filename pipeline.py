#!/usr/bin/env python3
"""Minimal SNP Panels: единый пользовательский интерфейс для FASTQ → VCF.

Варианты использования:
    python pipeline.py --fastq inputs/beluga_two.fastq.tsv --run-id test --backend slurm
    python pipeline.py status --run-id test
    python pipeline.py --vcf cohort.raw.bcf
    python pipeline.py --vcf cohort.raw.bcf --labels labels.tsv

Последние два сценария — честные заглушки до разработки научного анализа.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

# При запуске из исходников не требуется предварительное pip install -e .
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'src'))

from minsnp import analysis_stub, slurm_wgs, wgs
from minsnp.core import capture_versions, file_exists, identifier
from minsnp.resources import get_parallelism, job_allocation, memory_mib, validate_resources


def build_parser() -> argparse.ArgumentParser:
    """Определить единственный публичный CLI и служебные команды SLURM."""
    parser = argparse.ArgumentParser(
        description='Minimal SNP Panels: FASTQ → multi-sample VCF/BCF; VCF-анализ пока заглушка',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('command', nargs='?', choices=['run', 'status', '_worker'], default='run',
                        help='run — обычный запуск, status — проверка SLURM')
    sources = parser.add_mutually_exclusive_group()
    sources.add_argument('--fastq', nargs='+', metavar='TSV',
                         help='Одна или несколько таблиц с FASTQ, виды берутся из species')
    sources.add_argument('--vcf', metavar='VCF', help='Входной индексированный multi-sample BCF/VCF.gz')
    parser.add_argument('--labels', help='Таблица sample_id/population для будущего отбора SNP')
    parser.add_argument('--run-id', help='Уникальный идентификатор вычислительного запуска')
    parser.add_argument('--outdir', default='results', help='Корень выходных каталогов')
    parser.add_argument('--settings', default=str(ROOT / 'config/pipeline.yaml'),
                        help='YAML с ресурсами WGS/joint и настройками')
    parser.add_argument('--references', default=str(ROOT / 'config/references.yaml'),
                        help='YAML с путями к референсам каждого вида')
    parser.add_argument('--backend', choices=['local', 'slurm'], default='local')
    parser.add_argument('--submit', action='store_true',
                        help='Отправить SLURM jobs; без флага показывается только план')
    parser.add_argument('--resume', action='store_true',
                        help='Продолжить только существующий запуск с неизменными входами')
    parser.add_argument('--check-tools', action='store_true', help='Показать версии установленных программ')
    # Эти параметры использует только автоматически созданный SLURM worker.
    parser.add_argument('--plan', help=argparse.SUPPRESS)
    parser.add_argument('--group', help=argparse.SUPPRESS)
    parser.add_argument('--worker-species', help=argparse.SUPPRESS)
    parser.add_argument('--index', type=int, default=0, help=argparse.SUPPRESS)
    return parser


def resolve_mode(args: argparse.Namespace) -> str:
    """Режим однозначно определяется видом входных файлов."""
    if args.fastq:
        return 'fastq_to_vcf'
    if args.vcf and args.labels:
        return 'snp_panel'
    if args.vcf:
        return 'population_structure'
    raise ValueError('Укажите --fastq или --vcf; для справки используйте --help')


def load_settings(filename: str) -> dict:
    """Считать YAML и проверить ресурсы до создания каталогов результатов."""
    import yaml
    cfg = yaml.safe_load(file_exists(filename).read_text())
    if not isinstance(cfg, dict):
        raise ValueError(f'Настройки должны быть YAML-словарём: {filename}')
    validate_resources(cfg)
    return cfg


def print_allocation(cfg: dict, cohort: dict) -> None:
    """Показать максимальные заявки, не навязывая глобальный лимит CPU/RAM."""
    sample_n = sum(len(items) for items in cohort.values())
    active = min(sample_n, get_parallelism(cfg))
    scpu, smem = job_allocation(cfg, 'sample')
    jcpu, jmem = job_allocation(cfg, 'joint')
    print(f'Образцы: {sample_n}; виды: {len(cohort)}; одновременно WGS: {active}')
    print(f'WGS: на животное {scpu} CPU / {smem}; '
          f'максимум {active * scpu} CPU / {active * memory_mib(smem)} MiB')
    print(f'Joint: по одному заданию на вид, {jcpu} CPU / {jmem}; '
          f'максимум {len(cohort) * jcpu} CPU / {len(cohort) * memory_mib(jmem)} MiB')
    print('Хромосомы внутри WGS и joint выполняются последовательно.')


def run_fastq(args: argparse.Namespace) -> int:
    """Подготовить план и запустить локальный WGS либо общий SLURM array."""
    if not args.run_id:
        raise ValueError('FASTQ → VCF требует --run-id')
    identifier(args.run_id)
    if args.labels:
        raise ValueError('--labels используется только вместе с --vcf')
    if args.submit and args.backend != 'slurm':
        raise ValueError('--submit допускается только при --backend slurm')
    cfg = load_settings(args.settings)
    cohort = wgs.load_cohort(args.fastq)
    print_allocation(cfg, cohort)
    run_dir = Path(args.outdir).expanduser().resolve() / args.run_id
    if args.backend == 'slurm' and not args.submit:
        print('DRY RUN: задания не отправлены и каталоги результатов не изменены.')
        return 0
    plan = slurm_wgs.prepare_plan(run_dir, args.fastq, args.references,
                                  args.settings, cfg, resume=args.resume)
    versions_file = run_dir / 'versions.json'
    if not versions_file.exists():
        versions = capture_versions(['bwa-mem3', 'samtools', 'bcftools', 'gatk', 'fastp', 'mosdepth'])
        versions_file.write_text(json.dumps(versions, ensure_ascii=False, indent=2) + '\n')
    if args.backend == 'slurm':
        if args.resume and (run_dir / 'slurm_jobs.json').exists():
            raise ValueError('Задания уже отправлены; используйте status, не повторяйте --submit')
        from minsnp.core import require_commands
        require_commands('sbatch')
        jobs = slurm_wgs.submit_jobs(run_dir, cfg, sys.executable)
        print('Отправлены SLURM jobs:', json.dumps(jobs, ensure_ascii=False))
        return 0
    # Один и тот же worker для локального и SLURM запуска.
    plan_obj = json.loads(plan.read_text())
    for i in range(len(plan_obj['samples'])):
        slurm_wgs.worker(str(plan), 'sample', index=i)
    for species in sorted(cohort):
        slurm_wgs.worker(str(plan), 'joint', species=species)
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.check_tools:
        print(json.dumps(capture_versions(['bwa-mem3', 'samtools', 'bcftools', 'gatk', 'fastp', 'mosdepth']),
                         indent=2, ensure_ascii=False))
        return 0
    if args.command == '_worker':
        if not args.plan or not args.group:
            raise ValueError('SLURM worker требует --plan и --group')
        slurm_wgs.worker(args.plan, args.group, args.index, args.worker_species)
        return 0
    if args.command == 'status':
        if not args.run_id:
            raise ValueError('Для status требуется --run-id')
        slurm_wgs.status(Path(args.outdir).resolve() / identifier(args.run_id))
        return 0
    mode = resolve_mode(args)
    if mode == 'fastq_to_vcf':
        return run_fastq(args)
    if args.submit or args.resume:
        raise ValueError('--submit и --resume предназначены для FASTQ-режима')
    analysis_stub.pending_analysis(args.vcf, args.labels)
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except NotImplementedError as exc:
        print(f'НЕ РЕАЛИЗОВАНО: {exc}', file=sys.stderr)
        sys.exit(2)
    except (OSError, RuntimeError, ValueError) as exc:
        print(f'ОШИБКА: {exc}', file=sys.stderr)
        sys.exit(1)
