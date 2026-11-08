"""Контрактные тесты границ данных, безопасности запуска и аналитики."""
from __future__ import annotations

import json
from pathlib import Path
import sys

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from minsnp import panels, slurm_wgs, wgs


def _fastq_manifest(tmp_path: Path, entries: list[tuple]) -> Path:
    manifest = tmp_path / 'manifest.tsv'
    manifest.write_text('sample_id\tspecies\trun_id\tlibrary\tfastq_r1\tfastq_r2\n' +
                        ''.join('\t'.join(map(str, row)) + '\n' for row in entries))
    return manifest


def test_fastq_validation_shared_between_cohort_and_sample(tmp_path):
    fastq = [tmp_path / f'r{i}.fq.gz' for i in range(4)]
    for read in fastq:
        read.write_bytes(b'not a real FASTQ; validation checks metadata only')
    manifest = _fastq_manifest(tmp_path, [
        ('B1', 'beluga', 'R1', 'LIB1', fastq[0].name, fastq[1].name),
        ('B1', 'beluga', 'R2', 'LIB1', fastq[2].name, fastq[3].name),
    ])
    cohort = wgs.load_cohort([manifest])
    assert cohort['beluga'] == wgs.validate_fastq(str(manifest), 'beluga')
    assert len(cohort['beluga']['B1']) == 2
    assert cohort['beluga']['B1'][0]['fastq_r1'] == str(fastq[0].resolve())
    with pytest.raises(ValueError, match='Неверный вид'):
        wgs.validate_fastq(str(manifest), 'polar_bear')
    manifest.write_text(manifest.read_text() +
                        f'B1\tbeluga\tR3\tLIB1\t{fastq[0]}\t{fastq[1]}\n')
    with pytest.raises(ValueError, match='FASTQ используется'):
        wgs.load_cohort([manifest])


def test_concurrent_submission_lock_blocks_duplicate_sbatch(tmp_path, monkeypatch):
    run = tmp_path / 'run'
    run.mkdir()
    (run / 'slurm_plan.json').write_text(json.dumps({'stages': ['align'],
                                                     'tasks': {'align': [{'sample_index': 0}]}}))
    (run / '.submission.lock').mkdir()
    cfg = yaml.safe_load((ROOT / 'config/pipeline.yaml').read_text())
    monkeypatch.setattr(slurm_wgs, 'verify_plan', lambda _: None)
    def no_sbatch(*args, **kwargs):
        raise AssertionError('sbatch must not be called when run is locked')
    monkeypatch.setattr(slurm_wgs.subprocess, 'run', no_sbatch)
    with pytest.raises(ValueError, match='уже выполняется'):
        slurm_wgs.submit_jobs(run, cfg, sys.executable)


def test_panel_is_explicit_stub(tmp_path):
    """Режим панели не выполняет вычисления и не создаёт результатов."""
    import pipeline

    vcf = tmp_path / "cohort.bcf"
    metadata = tmp_path / "labels.tsv"

    with pytest.raises(NotImplementedError, match="не реализован"):
        panels.run(vcf, metadata)

    with pytest.raises(NotImplementedError, match="не реализован"):
        pipeline.main([
            "--vcf", str(vcf),
            "--metadata", str(metadata),
        ])

    assert not list(tmp_path.iterdir())


def test_slurm_dry_run_validates_reference_without_creating_outputs(tmp_path, capsys):
    """Dry run не требует GATK, но должен отвергать неверный список аутосом."""
    sys.path.insert(0, str(ROOT))
    import pipeline

    fasta = tmp_path / 'reference.fa'
    fasta.write_text('>chr1\nACGT\n')
    (tmp_path / 'reference.fa.fai').write_text('chr1\t4\t6\t4\t5\n')
    refs_file = tmp_path / 'references.yaml'
    refs_file.write_text(yaml.safe_dump({'species': {'beluga': {
        'assembly': 'test', 'fasta': str(fasta), 'contigs': ['chr1'],
    }}}))
    reads = [tmp_path / name for name in ('r1.fq.gz', 'r2.fq.gz')]
    for fastq in reads:
        fastq.write_bytes(b'dummy')
    manifest = _fastq_manifest(tmp_path, [
        ('B1', 'beluga', 'R1', 'LIB1', reads[0], reads[1]),
    ])
    args = [
        '--fastq', str(manifest), '--backend', 'slurm', '--run-id', 'dry',
        '--references', str(refs_file), '--settings', str(ROOT / 'config/pipeline.yaml'),
        '--outdir', str(tmp_path / 'results'),
    ]
    assert pipeline.main(args) == 0
    assert 'DRY RUN' in capsys.readouterr().out
    assert not (tmp_path / 'results').exists()

    spec = yaml.safe_load(refs_file.read_text())
    spec['species']['beluga']['contigs'] = ['missing_chr']
    refs_file.write_text(yaml.safe_dump(spec))
    with pytest.raises(ValueError, match='Contigs отсутствуют'):
        pipeline.main(args)
