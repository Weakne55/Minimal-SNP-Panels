#!/usr/bin/env python3
"""Minimal SNP Panels — FASTQ→VCF, VCF→structure (stub), VCF+metadata→panel.

Примеры:
  python pipeline.py --fastq inputs/fastq.tsv --run-id beluga1 --backend slurm --submit
  python pipeline.py status --run-id beluga1
  python pipeline.py --vcf cohort.raw.bcf                       # structure stub
  python pipeline.py --vcf cohort.raw.bcf --metadata labels.tsv --species beluga --run-id panel1
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'src'))

from minsnp import slurm_wgs, wgs
from minsnp.core import capture_versions, config, file_exists, identifier
from minsnp.resources import STAGE_ORDER, get_parallelism, job_allocation, memory_mib, validate_resources

LOG = logging.getLogger('minsnp')
TOOLS = ['bwa-mem3', 'samtools', 'bcftools', 'gatk', 'fastp', 'mosdepth']


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='FASTQ→VCF, VCF→structure (stub), VCF+metadata→SNP panel (stub)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('command', nargs='?', choices=['run', 'status', '_worker'], default='run')
    source = parser.add_mutually_exclusive_group()
    source.add_argument('--fastq', nargs='+', metavar='TSV', help='FASTQ manifests')
    source.add_argument('--vcf', metavar='VCF', help='Индексированный multi-sample .bcf/.vcf.gz')
    parser.add_argument('--metadata', '--labels', dest='metadata',
                        help='TSV с sample_id/population (режим пока заглушка)')
    parser.add_argument('--species', help='Зарезервировано для будущего режима SNP panel')
    parser.add_argument('--run-id', help='Идентификатор отдельного запуска')
    parser.add_argument('--outdir', default='results')
    parser.add_argument('--settings', default=str(ROOT / 'config/pipeline.yaml'))
    parser.add_argument('--references', default=str(ROOT / 'config/references.yaml'))
    parser.add_argument('--backend', choices=['local', 'slurm'], default='local')
    parser.add_argument('--submit', action='store_true', help='Для FASTQ: отправить SLURM arrays')
    parser.add_argument('--resume', action='store_true', help='Возобновить существующий local run')
    parser.add_argument('--check-tools', action='store_true')
    parser.add_argument('--plan', help=argparse.SUPPRESS)
    parser.add_argument('--stage', help=argparse.SUPPRESS)
    parser.add_argument('--index', type=int, default=0, help=argparse.SUPPRESS)
    return parser


def resolve_mode(args: argparse.Namespace) -> str:
    if args.fastq:
        return 'fastq_to_vcf'
    if args.vcf and args.metadata:
        return 'snp_panel'
    if args.vcf:
        return 'population_structure'
    raise ValueError('Укажите --fastq или --vcf')


def load_settings(filename: str) -> dict:
    import yaml
    cfg = yaml.safe_load(file_exists(filename).read_text())
    if not isinstance(cfg, dict):
        raise ValueError(f'Настройки должны быть YAML-словарём: {filename}')
    validate_resources(cfg)
    return cfg


def print_allocation(cfg: dict, cohort: dict) -> None:
    sample_n = sum(len(items) for items in cohort.values())
    lane_n = sum(len(lanes) for samples in cohort.values() for lanes in samples.values())
    print(f'Образцы: {sample_n}; lanes: {lane_n}; виды: {len(cohort)}')
    cap = get_parallelism(cfg)
    print(f'Максимум одновременно работающих SLURM задач в этом run-id: {cap}')
    for stage in STAGE_ORDER:
        if stage == 'trim' and not cfg.get('trim_fastq', False):
            continue
        cpus, ram = job_allocation(cfg, stage)
        print(f'  {stage}: {cpus} CPU, {ram} per task; максимум {cap * cpus} CPU, '
              f'{cap * memory_mib(ram)} MiB для полного массива')
    print('HaplotypeCaller: один gVCF на sample (все выбранные аутосомы).')
    print('GenomicsDB/GenotypeGVCFs: по одному заданию на вид, все аутосомы; барьеры afterok.')


def run_fastq(args: argparse.Namespace) -> int:
    if not args.run_id:
        raise ValueError('Для FASTQ → VCF нужен --run-id')
    identifier(args.run_id)
    if args.metadata or args.species:
        raise ValueError('--metadata и --species не используются с --fastq')
    if args.submit and args.backend != 'slurm':
        raise ValueError('--submit допускается только с --backend slurm')
    cfg = load_settings(args.settings)
    cohort = wgs.load_cohort(args.fastq)
    # В сухом запуске тоже проверяем ссылки на FASTA и список аутосом.
    for species in cohort:
        config(Path(args.references), Path(args.settings), species)
    print_allocation(cfg, cohort)
    run = Path(args.outdir).expanduser().resolve() / args.run_id
    if args.backend == 'slurm' and not args.submit:
        print('DRY RUN: никаких файлов и заданий не создано.')
        return 0
    plan_path = slurm_wgs.prepare_plan(run, args.fastq, args.references,
                                       args.settings, cfg, resume=args.resume)
    slurm_wgs.setup_logging(cfg.get('logging', {}).get('level', 'INFO'), run / 'logs' / 'pipeline.log')
    LOG.info('FASTQ mode run=%s backend=%s submit=%s samples=%d',
             run, args.backend, args.submit, sum(map(len, cohort.values())))
    versions_file = run / 'versions.json'
    if not versions_file.exists():
        versions_file.write_text(json.dumps(capture_versions(TOOLS), ensure_ascii=False, indent=2) + '\n')
    if args.backend == 'slurm':
        from minsnp.core import require_commands
        require_commands('sbatch')
        jobs = slurm_wgs.submit_jobs(run, cfg, sys.executable)
        print('Отправлены SLURM arrays:', json.dumps(jobs, ensure_ascii=False))
        return 0
    plan = json.loads(plan_path.read_text())
    for stage_name in plan['stages']:
        for index in range(len(plan['tasks'][stage_name])):
            slurm_wgs.worker(str(plan_path), stage_name, index=index)
    LOG.info('FASTQ→VCF pipeline completed run=%s', run)
    return 0


def run_panel(args: argparse.Namespace) -> int:
    """Режим SNP panel пока не реализован."""
    from minsnp import panels
    panels.run(Path(args.vcf), Path(args.metadata))
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.check_tools:
        print(json.dumps(capture_versions(TOOLS), indent=2, ensure_ascii=False))
        return 0
    if args.command == '_worker':
        if not args.plan or not args.stage:
            raise ValueError('SLURM worker требует --plan и --stage')
        slurm_wgs.worker(args.plan, args.stage, args.index)
        return 0
    if args.command == 'status':
        if not args.run_id:
            raise ValueError('Для status требуется --run-id')
        slurm_wgs.status(Path(args.outdir).resolve() / identifier(args.run_id))
        return 0
    mode = resolve_mode(args)
    if mode == 'fastq_to_vcf':
        return run_fastq(args)
    if mode == 'snp_panel':
        return run_panel(args)
    if args.submit or args.resume or args.backend != 'local':
        raise ValueError('VCF→structure сейчас заглушка; --submit/--resume/SLURM не поддерживаются')
    from minsnp import structure
    structure.run(Path(args.vcf))
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except NotImplementedError as exc:
        print(f'НЕ РЕАЛИЗОВАНО: {exc}', file=sys.stderr)
        sys.exit(2)
    except Exception as exc:
        LOG.exception("FAILED pipeline: %s", exc)
        print(f'ОШИБКА: {exc}', file=sys.stderr)
        sys.exit(1)
