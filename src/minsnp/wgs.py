"""Обработка FASTQ по sequencing units и получение одного gVCF на образец.

FASTQ-манифест содержит идентификаторы sample/species/run/library и парные
FASTQ. Несколько sequencing units могут относиться к одному образцу.
"""
from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Callable

from .core import check_index, file_exists, identifier, read_tsv, stage, write_tsv
from .resources import stage_effective

FASTQ_COLUMNS = (
    'sample_id', 'species', 'run_id', 'library', 'fastq_r1', 'fastq_r2',
)
BWA_INDEX_SUFFIXES = ('.amb', '.ann', '.pac', '.bwt.2bit.64')


def load_cohort(sheets: list[str | Path]) -> dict[str, dict[str, list[dict[str, str]]]]:
    """Проверить входные TSV и объединить sequencing units по виду/образцу.

    Пути к FASTQ интерпретируются относительно соответствующего TSV. Один
    FASTQ и одна тройка (sample, run, library) могут встретиться только раз.
    """
    cohort: dict[str, dict[str, list[dict[str, str]]]] = {}
    sample_species: dict[str, str] = {}
    used_units: set[tuple[str, str, str]] = set()
    used_fastq: set[Path] = set()

    for sheet in sheets:
        manifest = Path(sheet).expanduser().resolve()
        for row in read_tsv(manifest, FASTQ_COLUMNS):
            sample, species, run_id, library = (
                identifier(row[field]) for field in FASTQ_COLUMNS[:4]
            )
            if sample in sample_species and sample_species[sample] != species:
                raise ValueError(f'Образец {sample} встречается у нескольких видов')
            sample_species[sample] = species
            unit = (sample, run_id, library)
            if unit in used_units:
                raise ValueError(f'Повторная строка sample/run/library: {unit}')
            used_units.add(unit)

            resolved = dict(row)
            for field in ('fastq_r1', 'fastq_r2'):
                fastq = Path(row[field]).expanduser()
                if not fastq.is_absolute():
                    fastq = manifest.parent / fastq
                fastq = file_exists(fastq)
                if fastq in used_fastq:
                    raise ValueError(f'FASTQ используется более одного раза: {fastq}')
                used_fastq.add(fastq)
                resolved[field] = str(fastq)
            cohort.setdefault(species, {}).setdefault(sample, []).append(resolved)

    if not cohort:
        raise ValueError('Нет животных в FASTQ-манифестах')
    return cohort


def validate_fastq(sheet: str, species: str) -> dict[str, list[dict[str, str]]]:
    """Прочитать манифест одного вида, используя общую проверку FASTQ."""
    cohort = load_cohort([sheet])
    if set(cohort) != {species}:
        raise ValueError(f'Неверный вид в FASTQ TSV: ожидался {species}, получены {sorted(cohort)}')
    return cohort[species]


def quarantine_incomplete(path: Path) -> None:
    """Сохранить повреждённый выход и его индексы с суффиксом .invalid_N."""
    for target in (path, *(Path(str(path) + suffix) for suffix in ('.bai', '.csi', '.tbi'))):
        if target.exists() or target.is_symlink():
            version = 1
            while Path(f'{target}.invalid_{version}').exists():
                version += 1
            target.rename(Path(f'{target}.invalid_{version}'))


def _valid_bam(path: Path, index: bool = False) -> bool:
    """Проверить структуру BAM (не полный аудит прочтений)."""
    if not path.is_file() or path.stat().st_size == 0:
        return False
    if index and not any(Path(f'{path}{suffix}').is_file() for suffix in ('.bai', '.csi')):
        return False
    return subprocess.run(['samtools', 'quickcheck', '-v', str(path)], capture_output=True).returncode == 0


def _valid_gvcf(path: Path) -> bool:
    """Проверить наличие индекса и читаемость заголовка/индекса BCF/VCF."""
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        check_index(path)
    except FileNotFoundError:
        return False
    return subprocess.run(['bcftools', 'index', '-n', str(path)], capture_output=True).returncode == 0


def _require_bwa_mem3_index(fasta: str) -> None:
    missing = [f'{fasta}{suffix}' for suffix in BWA_INDEX_SUFFIXES
               if not Path(f'{fasta}{suffix}').is_file() or Path(f'{fasta}{suffix}').stat().st_size == 0]
    if missing:
        raise ValueError('Отсутствует индекс BWA-MEM3: ' + ', '.join(missing)
                         + f'. Выполните bwa-mem3 index {fasta!r}')


def _lane_dir(out: Path, sample: str, row: dict[str, str]) -> Path:
    return out / 'lanes' / sample / row['run_id'] / row['library']


def _lane_reads(out: Path, sample: str, row: dict[str, str], trim: bool) -> tuple[str, str]:
    if not trim:
        return row['fastq_r1'], row['fastq_r2']
    lane = _lane_dir(out, sample, row)
    r1, r2, report = lane / 'trim_R1.fastq.gz', lane / 'trim_R2.fastq.gz', lane / 'trim_R1.json'
    if not all(p.is_file() and p.stat().st_size > 0 for p in (r1, r2, report)):
        raise FileNotFoundError(f'Отсутствует завершённый trim: {lane}')
    return str(r1), str(r2)


def _run_checked(output: Path, valid: Callable[[Path], bool], script: str,
                 arguments: list, logs: Path, label: str, cfg: dict) -> Path:
    """Переиспользовать корректный выход; иначе сохранить старый и пересоздать."""
    if not valid(output):
        quarantine_incomplete(output)
        stage(script, arguments, logs, label, cfg=cfg)
    if not valid(output):
        raise RuntimeError(f'Стадия {label} создала некорректный результат: {output}')
    return output


def _run_lane(name: str, out: Path, logs: Path, sample: str, rows: list[dict[str, str]],
              lane_id: tuple[str, str] | None, fasta: str, cpus: int, cfg: dict) -> Path:
    if lane_id is None:
        raise ValueError(f'{name}: отсутствует lane_id')
    matching = [r for r in rows if (r['run_id'], r['library']) == lane_id]
    if len(matching) != 1:
        raise ValueError(f'{name}: не найдена sequencing unit {lane_id}')
    row = matching[0]
    lane = _lane_dir(out, sample, row)
    lane.mkdir(parents=True, exist_ok=True)
    label = f"{sample}_{row['run_id']}_{row['library']}"

    if name == 'trim':
        r1, r2, report = lane / 'trim_R1.fastq.gz', lane / 'trim_R2.fastq.gz', lane / 'trim_R1.json'
        if not all(p.is_file() and p.stat().st_size for p in (r1, r2, report)):
            for path in (r1, r2, report):
                quarantine_incomplete(path)
            stage('trim.sh', [row['fastq_r1'], row['fastq_r2'], r1, r2, cpus],
                  logs, f'trim_{label}', cfg=cfg)
        _lane_reads(out, sample, row, True)
        return r1

    _require_bwa_mem3_index(fasta)
    reads = _lane_reads(out, sample, row, bool(cfg.get('trim_fastq', False)))
    output = lane / 'fixed.bam'
    return _run_checked(output, _valid_bam, 'align.sh',
                        [fasta, *reads, output, sample, row['library'], row['run_id'], cpus],
                        logs, f'align_{label}', cfg)


def _run_sample(name: str, out: Path, logs: Path, sample: str,
                rows: list[dict[str, str]], spec: dict, cpus: int, cfg: dict) -> Path:
    bam_dir = out / 'bam'
    bam_dir.mkdir(parents=True, exist_ok=True)
    merged = bam_dir / f'{sample}.merged.bam'
    final = bam_dir / f'{sample}.bam'

    if name == 'merge_bam':
        lane_bams = [_lane_dir(out, sample, row) / 'fixed.bam' for row in rows]
        for lane_bam in lane_bams:
            if not _valid_bam(lane_bam):
                raise RuntimeError(f'Нет корректного lane BAM: {lane_bam}')
        return _run_checked(merged, _valid_bam, 'merge_bam.sh',
                            [merged, cpus, *lane_bams], logs, f'merge_{sample}', cfg)

    if name == 'markdup':
        if not _valid_bam(merged):
            raise RuntimeError(f'Отсутствует merge BAM: {merged}')
        return _run_checked(final, lambda p: _valid_bam(p, index=True),
                            'markdup.sh', [merged, final, out / 'qc' / f'{sample}.dup.txt', cpus],
                            logs, f'markdup_{sample}', cfg)

    if not _valid_bam(final, index=True):
        raise RuntimeError(f'Отсутствует финальный BAM: {final}')
    if name == 'bam_qc':
        prefix = out / 'qc' / sample
        reports = [Path(f'{prefix}.{suffix}') for suffix in
                   ('flagstat.txt', 'stats.txt', 'mosdepth.summary.txt')]
        if not all(p.is_file() and p.stat().st_size > 0 for p in reports):
            stage('bam_qc.sh', [final, prefix, cpus], logs, f'qc_{sample}', cfg=cfg)
        if not all(p.is_file() and p.stat().st_size > 0 for p in reports):
            raise RuntimeError(f'BAM QC не создал полные отчёты: {prefix}')
        return reports[0]

    if name == 'haplotypecaller':
        gvcf = out / 'gvcf' / f'{sample}.g.vcf.gz'
        gvcf.parent.mkdir(parents=True, exist_ok=True)
        result = _run_checked(gvcf, _valid_gvcf, 'haplotypecaller.sh',
                              [spec['fasta'], final, gvcf, cpus, *spec['contigs']],
                              logs, f'hc_{sample}', cfg)
        manifest = out / 'gvcf_manifest.tsv'
        tmp = manifest.with_name(manifest.name + '.tmp')
        write_tsv(tmp, ['sample_id', 'gvcf_path'], [(sample, str(result.resolve()))])
        tmp.replace(manifest)
        return result

    raise ValueError(f'Неизвестная WGS-стадия: {name}')


def run_stage(name: str, sheet: str, species: str, spec: dict, cfg: dict,
              out: Path, logs: Path, lane_id: tuple[str, str] | None = None) -> Path:
    """Выполнить одну WGS-стадию для sequencing unit либо образца."""
    grouped = validate_fastq(sheet, species)
    if len(grouped) != 1:
        raise ValueError(f'Ожидается манифест одного образца: {sheet}')
    sample, rows = next(iter(grouped.items()))
    rows.sort(key=lambda row: (row['run_id'], row['library']))
    cpus, _ = stage_effective(cfg, name)
    out.mkdir(parents=True, exist_ok=True)
    if name in ('trim', 'align'):
        return _run_lane(name, out, logs, sample, rows, lane_id, spec['fasta'], cpus, cfg)
    return _run_sample(name, out, logs, sample, rows, spec, cpus, cfg)
