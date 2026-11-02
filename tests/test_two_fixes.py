"""Регрессии для двух исправлений: SLURM time и прерванный план.

Тесты работают только во временных каталогах; sbatch и GATK не вызываются.
"""
from __future__ import annotations

from pathlib import Path

import pytest
from minsnp import resources, slurm_wgs


def config(time: str = '02:00:00') -> dict:
    return {'slurm': {
        'max_parallel_samples': 2,
        'sample': {'cpus': 8, 'memory': '16G', 'time': time},
        'joint': {'cpus': 4, 'memory': '8G', 'time': '12:00:00'},
    }}


@pytest.mark.parametrize('value', [
    '99:99:99', '12:60:00', '12:00:60', '2-24:00:00',
    '2-04:99:00', '2-04:00:99', '00:00:00', '-01:00:00',
])
def test_reject_invalid_slurm_time(value: str) -> None:
    with pytest.raises(ValueError, match='time'):
        resources.validate_resources(config(value))


@pytest.mark.parametrize('value', [
    '12:00:00', '99:00:00', '2-04:30:59', '0-23:59:59', '00:00:01',
])
def test_accept_valid_slurm_time(value: str) -> None:
    resources.validate_resources(config(value))


def inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    (tmp_path / 'r1.fastq.gz').write_bytes(b'FASTQ1')
    (tmp_path / 'r2.fastq.gz').write_bytes(b'FASTQ2')
    sheet = tmp_path / 'samples.tsv'
    sheet.write_text('sample_id\tspecies\trun_id\tlibrary\tfastq_r1\tfastq_r2\n'
                     'B1\tbeluga\tR1\tL1\tr1.fastq.gz\tr2.fastq.gz\n')
    refs = tmp_path / 'references.yaml'; refs.write_text('species: {}\n')
    settings = tmp_path / 'settings.yaml'; settings.write_text('trim_fastq: false\n')
    return sheet, refs, settings


def test_resume_after_partial_manifest(tmp_path: Path, monkeypatch) -> None:
    sheet, refs, settings = inputs(tmp_path)
    run = tmp_path / 'results' / 'test'
    monkeypatch.setattr(slurm_wgs, 'config', lambda *_: ({}, {'contigs': ['chr1']}))
    real_write = slurm_wgs.write_tsv

    def interrupted_write(path: Path, *args):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('НЕПОЛНАЯ ЗАПИСЬ\n')
        raise OSError('Смоделировано прерывание во время записи TSV')

    monkeypatch.setattr(slurm_wgs, 'write_tsv', interrupted_write)
    with pytest.raises(OSError):
        slurm_wgs.prepare_plan(run, [str(sheet)], str(refs), str(settings), config())

    monkeypatch.setattr(slurm_wgs, 'write_tsv', real_write)
    saved = slurm_wgs.prepare_plan(run, [str(sheet)], str(refs), str(settings),
                                    config(), resume=True)
    assert saved.is_file()
    manifest = run / 'manifests/beluga/B1.tsv'
    assert manifest.read_text().startswith('sample_id\tspecies\trun_id')
    assert not list(run.rglob('*.tmp'))
    # Повторное восстановление не должно менять уже готовый план.
    before = saved.read_bytes()
    assert slurm_wgs.prepare_plan(run, [str(sheet)], str(refs), str(settings),
                                  config(), resume=True).read_bytes() == before


def test_partial_resume_does_not_overwrite_unknown_files(tmp_path: Path, monkeypatch) -> None:
    sheet, refs, settings = inputs(tmp_path)
    run = tmp_path / 'results/test'
    run.mkdir(parents=True)
    other = run / 'other.txt'; other.write_text('Важные данные')
    monkeypatch.setattr(slurm_wgs, 'config', lambda *_: ({}, {'contigs': ['chr1']}))
    with pytest.raises(ValueError, match='посторонн|план|без плана'):
        slurm_wgs.prepare_plan(run, [str(sheet)], str(refs), str(settings),
                               config(), resume=True)
    assert other.read_text() == 'Важные данные'
    assert not (run / 'slurm_plan.json').exists()
