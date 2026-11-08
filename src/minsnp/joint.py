"""Совместное генотипирование всех аутосом: одна GenomicsDB и один callset на вид."""
from __future__ import annotations

from pathlib import Path
import shutil

from .core import check_index, file_exists, samples_in_vcf, stage
from .wgs import _valid_gvcf, quarantine_incomplete


def valid_variant(path: Path) -> bool:
    return _valid_gvcf(path)


def sample_map_from_plan(plan: dict, species: str, run: Path) -> dict[str, str]:
    """Проверять входные gVCF и состав вида перед импортом."""
    mapping = {}
    for item in plan['samples']:
        if item['species'] != species:
            continue
        sample = item['sample']
        gvcf = (run / 'samples' / species / sample / 'wgs' / 'gvcf' /
                f'{sample}.g.vcf.gz')
        file_exists(gvcf)
        check_index(gvcf)
        if not _valid_gvcf(gvcf):
            raise ValueError(f"Повреждён gVCF: {gvcf}")
        actual = samples_in_vcf(gvcf)
        if actual != [sample]:
            raise ValueError(f"gVCF {gvcf}: sample columns {actual}, ожидался {sample}")
        if sample in mapping:
            raise ValueError(f"Повторный образец: {sample}")
        mapping[sample] = str(gvcf.resolve())
    if not mapping:
        raise ValueError(f"Не найдены образцы вида {species}")
    return mapping


def _setup_sample_map(plan: dict, species: str, run: Path, out: Path) -> Path:
    mapping = sample_map_from_plan(plan, species, run)
    sample_map = out / 'sample_map.tsv'
    sample_map.parent.mkdir(parents=True, exist_ok=True)
    expected = ''.join(f'{sample}\t{mapping[sample]}\n' for sample in sorted(mapping))
    if sample_map.exists():
        if sample_map.read_text() != expected:
            raise ValueError(f"Изменился состав sample-map: {sample_map}")
    else:
        tmp = sample_map.with_name(sample_map.name + '.tmp')
        tmp.write_text(expected)
        tmp.replace(sample_map)
    return sample_map


def _db_valid(db: Path) -> bool:
    return db.is_dir() and (db / '.completed').is_file()


def _check_samples(plan: dict, species: str, path: Path) -> None:
    expected = {s['sample'] for s in plan['samples'] if s['species'] == species}
    actual = samples_in_vcf(path)
    if len(actual) != len(expected) or set(actual) != expected:
        raise ValueError(f"Не совпадает состав образцов в {path}: {actual}")


def run_stage(name: str, plan: dict, species: str, spec: dict,
              run: Path, cfg: dict, out: Path, logs: Path) -> Path:
    """Каждая стадия выполняется один раз на вид, без per-contig shards."""
    out.mkdir(parents=True, exist_ok=True)
    contigs = spec['contigs']
    if not contigs or len(set(contigs)) != len(contigs):
        raise ValueError(f"Некорректный список аутосом вида {species}")
    db = out / 'genomicsdb'
    bcf = out / 'cohort.raw.bcf'
    vcf = out / 'cohort.raw.vcf.gz'
    if name == 'genomicsdb':
        samplemap = _setup_sample_map(plan, species, run, out)
        if not _db_valid(db):
            if db.exists() or db.is_symlink():
                if db.is_symlink() or db.parent.resolve() != out.resolve():
                    raise ValueError(f"Небезопасный путь GenomicsDB: {db}")
                shutil.rmtree(db)
            stage('genomicsdb.sh', [spec['fasta'], samplemap, db, *contigs],
                  logs, 'import_autosomes', cfg=cfg)
        if not _db_valid(db):
            raise RuntimeError(f"GenomicsDB не завершился: {db}")
        return db
    if name == 'genotype':
        if not _db_valid(db):
            raise RuntimeError(f"Отсутствует готовый GenomicsDB: {db}")
        if not valid_variant(bcf):
            quarantine_incomplete(bcf)
            stage('genotype.sh', [spec['fasta'], db, bcf, *contigs],
                  logs, 'genotype_autosomes', cfg=cfg)
        if not valid_variant(bcf):
            raise RuntimeError(f"Некорректный общий BCF: {bcf}")
        _check_samples(plan, species, bcf)
        return bcf
    if name == 'export_vcf':
        if not valid_variant(bcf):
            raise RuntimeError(f"Отсутствует итоговый BCF: {bcf}")
        if not valid_variant(vcf):
            quarantine_incomplete(vcf)
            stage('export_vcf.sh', [bcf, vcf], logs, 'export_vcf', cfg=cfg)
        if not valid_variant(vcf):
            raise RuntimeError(f"Некорректный итоговый VCF: {vcf}")
        _check_samples(plan, species, vcf)
        return vcf
    raise ValueError(f"Неизвестная joint стадия: {name}")
