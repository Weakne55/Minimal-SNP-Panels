"""Контрактные тесты нового запуска, без настоящего sbatch или обработки FASTQ."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest
from minsnp import resources, slurm_wgs, wgs

ROOT = Path(__file__).resolve().parents[1]


def load_entry():
    spec = importlib.util.spec_from_file_location('minsnp_entrypoint', ROOT / 'pipeline.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sample_config():
    return {'slurm': {'partition': 'main', 'nodelist': 'meduza-1',
                      'max_parallel_samples': 3,
                      'sample': {'cpus': 16, 'memory': '64G', 'time': '2-00:00:00'},
                      'joint': {'cpus': 8, 'memory': '48G', 'time': '12:00:00'}}}


def test_cli_single_entrypoint_modes():
    entry = load_entry()
    args = entry.build_parser().parse_args(['--fastq', 'one.tsv', 'two.tsv', '--run-id', 'r1'])
    assert args.fastq == ['one.tsv', 'two.tsv']
    assert entry.resolve_mode(args) == 'fastq_to_vcf'
    assert entry.resolve_mode(entry.build_parser().parse_args(['--vcf', 'a.bcf'])) == 'population_structure'
    assert entry.resolve_mode(entry.build_parser().parse_args(['--vcf', 'a.bcf', '--labels', 'l.tsv'])) == 'snp_panel'
    with pytest.raises(SystemExit):
        entry.build_parser().parse_args(['--fastq', 'x', '--vcf', 'y'])


def test_memory_and_per_group_resources():
    cfg = sample_config()
    resources.validate_resources(cfg)
    assert resources.stage_allocation(cfg, 'align.sh') == (16, 65536)
    assert resources.stage_allocation(cfg, 'haplotypecaller.sh') == (16, 65536)
    assert resources.stage_allocation(cfg, 'genotype.sh') == (8, 49152)
    assert resources.get_parallelism(cfg) == 3
    assert resources.job_allocation(cfg, 'sample') == (16, '64G')
    assert resources.job_allocation(cfg, 'joint') == (8, '48G')


def test_reject_bad_resources():
    cfg = sample_config()
    cfg['slurm']['sample']['cpus'] = 0
    with pytest.raises(ValueError):
        resources.validate_resources(cfg)
    cfg = sample_config()
    cfg['slurm']['max_parallel_samples'] = 0
    with pytest.raises(ValueError):
        resources.validate_resources(cfg)


def test_joint_array_not_created_and_global_limit(tmp_path, monkeypatch):
    calls = []
    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        return SimpleNamespace(returncode=0, stdout=str(101 + len(calls) - 1) + '\n', stderr='')
    monkeypatch.setattr(slurm_wgs.subprocess, 'run', fake_run)
    plan = {'samples': [{'species': 'beluga', 'sample': 'B1'}, {'species': 'beluga', 'sample': 'B2'},
                        {'species': 'polar_bear', 'sample': 'P1'}, {'species': 'polar_bear', 'sample': 'P2'},
                        {'species': 'polar_bear', 'sample': 'P3'}]}
    (tmp_path / 'slurm_plan.json').write_text(json.dumps(plan))
    jobs = slurm_wgs.submit_jobs(tmp_path, sample_config(), sys.executable)
    assert len(calls) == 3
    assert '--array=0-4%3' in calls[0]
    assert '--cpus-per-task=16' in calls[0]
    assert '--mem=64G' in calls[0]
    for cmd in calls[1:]:
        assert not any(x.startswith('--array=') for x in cmd)
        assert '--dependency=afterok:101' in cmd
        assert '--cpus-per-task=8' in cmd
    assert set(jobs) == {'samples', 'joint_beluga', 'joint_polar_bear'}


def test_cohort_two_species_combined(tmp_path):
    def sheet(path, ident, species):
        r1 = tmp_path / f'{ident}_1.fastq.gz'; r1.write_bytes(b'x')
        r2 = tmp_path / f'{ident}_2.fastq.gz'; r2.write_bytes(b'x')
        path.write_text('sample_id\tspecies\trun_id\tlibrary\tfastq_r1\tfastq_r2\n'
                        f'{ident}\t{species}\tR{ident}\tL1\t{r1.name}\t{r2.name}\n')
    p1 = tmp_path / 'a.tsv'; p2 = tmp_path / 'b.tsv'
    sheet(p1, 'B1', 'beluga'); sheet(p2, 'P1', 'polar_bear')
    result = wgs.load_cohort([p1, p2])
    assert set(result) == {'beluga', 'polar_bear'}
    assert set(result['beluga']) == {'B1'}
    assert set(result['polar_bear']) == {'P1'}


def test_status_does_not_claim_success_on_empty_queue(tmp_path, monkeypatch, capsys):
    (tmp_path / 'slurm_jobs.json').write_text('{"samples":"123"}')
    (tmp_path / 'slurm_plan.json').write_text('{"samples":[]}')
    monkeypatch.setattr(slurm_wgs.subprocess, 'run', lambda *a, **kw: SimpleNamespace(returncode=0,stdout='',stderr=''))
    slurm_wgs.status(tmp_path)
    assert 'не подтверждён' in capsys.readouterr().out


def test_plan_snapshot_and_resume_no_duplicates(tmp_path, monkeypatch):
    refs = tmp_path / 'references.yaml'; refs.write_text('temporary references')
    settings = tmp_path / 'pipeline.yaml'; settings.write_text('temporary settings')
    fastq = tmp_path / 'read.fastq.gz'; fastq.write_bytes(b'123')
    mate = tmp_path / 'mate.fastq.gz'; mate.write_bytes(b'456')
    sheet = tmp_path / 'cohort.tsv'
    sheet.write_text('sample_id\tspecies\trun_id\tlibrary\tfastq_r1\tfastq_r2\n'
                     'B1\tbeluga\tR1\tL1\tread.fastq.gz\tmate.fastq.gz\n')
    monkeypatch.setattr(slurm_wgs, 'config', lambda *a: ({}, {'contigs': ['chr1']}))
    run = tmp_path / 'results' / 'r1'
    first = slurm_wgs.prepare_plan(run, [str(sheet)], str(refs), str(settings), sample_config())
    assert len(json.loads(first.read_text())['samples']) == 1
    assert first == slurm_wgs.prepare_plan(run, [str(sheet)], str(refs), str(settings),
                                          sample_config(), resume=True)
    with pytest.raises(ValueError, match='уже существует'):
        slurm_wgs.prepare_plan(run, [str(sheet)], str(refs), str(settings), sample_config())
    sheet.write_text(sheet.read_text() + '\n')
    with pytest.raises(ValueError, match='Изменился план'):
        slurm_wgs.prepare_plan(run, [str(sheet)], str(refs), str(settings),
                               sample_config(), resume=True)


def test_slurm_dry_run_keeps_output_untouched(tmp_path, monkeypatch, capsys):
    entry = load_entry()
    fastq = tmp_path / 'example.tsv'
    fastq.write_text('dummy')
    monkeypatch.setattr(entry, 'load_settings', lambda *a: sample_config())
    monkeypatch.setattr(entry.wgs, 'load_cohort', lambda *a: {'beluga': {'B1': []}})
    args = entry.build_parser().parse_args(['--fastq', str(fastq), '--backend', 'slurm',
                                            '--run-id', 'r1', '--outdir', str(tmp_path / 'results')])
    assert entry.run_fastq(args) == 0
    assert 'DRY RUN' in capsys.readouterr().out
    assert not (tmp_path / 'results').exists()


def test_failed_stage_exits_and_records_log(tmp_path, monkeypatch):
    from minsnp import core
    monkeypatch.setattr(core.subprocess, 'run', lambda *a, **kw: SimpleNamespace(returncode=9))
    with pytest.raises(RuntimeError, match='код 9'):
        core.call(['bash', 'fake.sh'], tmp_path / 'logs', 'failed')
    assert (tmp_path / 'logs/failed.log').is_file()
