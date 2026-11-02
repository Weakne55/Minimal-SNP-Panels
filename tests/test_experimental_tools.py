"""Контрактные тесты экспериментальной замены mosdepth и samtools markdup.

Fake executables проверяют вызовы Bash и выходные файлы, но не биологию.
"""
from __future__ import annotations

import gzip
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _fake_toolbox(path: Path) -> Path:
    """Создать имитации команд для тестирования shell orchestration."""
    path.mkdir(parents=True)
    fake = path / 'tool.py'
    fake.write_text('''#!/usr/bin/env python3
import gzip, os, pathlib, sys
name = pathlib.Path(sys.argv[0]).name
args = sys.argv[1:]
with open(os.environ['TRACE'], 'a') as log:
    log.write(name + ' ' + ' '.join(args) + '\\n')
if name == 'samtools':
    sub = args[0]
    if sub == 'quickcheck':
        sys.exit(0 if pathlib.Path(args[-1]).is_file() else 1)
    if sub == 'collate':
        sys.stdout.buffer.write(pathlib.Path(args[-2]).read_bytes())
    if sub == 'sort':
        pathlib.Path(args[args.index('-o') + 1]).write_bytes(sys.stdin.buffer.read())
    if sub == 'index':
        bam = args[-1]
        pathlib.Path(bam + '.bai').write_bytes(b'INDEX')
    if sub in {'flagstat', 'stats'}:
        print('mock statistics')
elif name == 'dupblaster':
    assert '--metrics-prefix' in args
    assert '--sequencing-duplicate-detection' in args
    assert args[args.index('--sequencing-duplicate-detection')+1] == 'off'
    prefix = pathlib.Path(args[args.index('--metrics-prefix') + 1])
    pathlib.Path(str(prefix)+'.duplicate-metrics.tsv').write_text('library\\tduplicate_pairs\\nL1\\t1\\n')
    sys.stdout.buffer.write(sys.stdin.buffer.read())
elif name == 'pandepth':
    assert '-i' in args and '-o' in args and '-t' in args
    prefix = args[args.index('-o') + 1]
    with gzip.open(prefix+'.chr.stat.gz', 'wt') as out:
        out.write('#Chr\\tLength\\tCoveredSite\\tTotalDepth\\tCoverage(%)\\tMeanDepth\\n')
        out.write('chr1\\t100\\t90\\t400\\t90.0\\t4.0\\n')
''')
    fake.chmod(0o755)
    for tool in ('samtools', 'dupblaster', 'pandepth'):
        (path / tool).symlink_to(fake)
    return path


def _env(tmp_path: Path) -> dict[str, str]:
    bindir = _fake_toolbox(tmp_path / 'bin')
    return dict(os.environ, PATH=f'{bindir}:{os.environ["PATH"]}',
                TRACE=str(tmp_path / 'trace.txt'), MIN_SNP_MEMORY_MB='16384',
                TMPDIR=str(tmp_path))


def test_dupblaster_after_collate_before_sort(tmp_path: Path) -> None:
    source = tmp_path / 'merged.bam'
    source.write_bytes(b'MOCK_BAM')
    result = tmp_path / 'output.bam'
    metrics = tmp_path / 'metrics' / 'sample.dupblaster'
    env = _env(tmp_path)
    subprocess.run(['bash', str(ROOT / 'stages/markdup.sh'),
                    str(source), str(result), str(metrics), '8'], env=env, check=True)
    assert result.read_bytes() == b'MOCK_BAM'
    assert (tmp_path / 'output.bam.bai').exists()
    assert Path(str(metrics) + '.duplicate-metrics.tsv').is_file()
    trace = (tmp_path / 'trace.txt').read_text().splitlines()
    collate = next(i for i, s in enumerate(trace) if s.startswith('samtools collate '))
    dupblaster = next(i for i, s in enumerate(trace) if s.startswith('dupblaster '))
    sort = next(i for i, s in enumerate(trace) if s.startswith('samtools sort '))
    # В pipe процессы стартуют одновременно: порядок логов не определён.
    # Проверяем состав команд и порядок потоков по Bash-скрипту отдельно.
    assert len({collate, dupblaster, sort}) == 3
    script = (ROOT / 'stages/markdup.sh').read_text()
    assert script.index('samtools collate -O') < script.index('| dupblaster') < script.index('| samtools sort')
    assert 'samtools markdup' not in script


def test_pandepth_outputs_match_resume_contract(tmp_path: Path) -> None:
    source = tmp_path / 'final.bam'
    source.write_bytes(b'MOCK_BAM')
    prefix = tmp_path / 'qc' / 'sample'
    env = _env(tmp_path)
    subprocess.run(['bash', str(ROOT / 'stages/bam_qc.sh'),
                    str(source), str(prefix), '8'], env=env, check=True)
    for suffix in ('.flagstat.txt', '.stats.txt', '.chr.stat.gz'):
        assert Path(str(prefix) + suffix).stat().st_size > 0
    with gzip.open(str(prefix) + '.chr.stat.gz', 'rt') as handle:
        assert 'MeanDepth' in handle.readline()
    assert 'pandepth ' in (tmp_path / 'trace.txt').read_text()
    assert 'mosdepth' not in (ROOT / 'stages/bam_qc.sh').read_text().lower()


def test_python_pipeline_expects_pandepth_and_dupblaster() -> None:
    code = (ROOT / 'src/minsnp/wgs.py').read_text()
    assert 'Path(f"{prefix}.chr.stat.gz")' in code
    assert 'f"{sample}.dupblaster"' in code
    assert 'mosdepth.summary.txt' not in code
    cli = (ROOT / 'pipeline.py').read_text()
    assert "'pandepth', 'dupblaster'" in cli


def test_test_environment_has_pinned_tools() -> None:
    import yaml
    spec = yaml.safe_load((ROOT / 'envs/wgs-gatk-mem3-test.yaml').read_text())
    assert spec['name'] == 'wgs-gatk-mem3-test'
    pan = yaml.safe_load((ROOT / 'envs/pandepth-test.yaml').read_text())
    assert pan['name'] == 'minsnp-pandepth'
    assert 'pandepth=2.26' in pan['dependencies']
    assert 'dupblaster=0.1.1' in spec['dependencies']
    assert not any(str(s).startswith('mosdepth=') for s in spec['dependencies'])


def test_pandepth_version_from_help(monkeypatch) -> None:
    from minsnp import core
    class Completed:
        stdout = "Usage: pandepth ...\n  -h  show this help [v2.26]\n"
        stderr = ""
    monkeypatch.setattr(core.shutil, 'which', lambda name: '/bin/' + name)
    monkeypatch.setattr(core.subprocess, 'run', lambda *args, **kwargs: Completed())
    assert core.capture_versions(['pandepth']) == {'pandepth': 'PanDepth 2.26'}
