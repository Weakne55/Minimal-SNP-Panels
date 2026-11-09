"""Место для будущего VCF → популяционная структура; научных результатов нет."""
from __future__ import annotations

from pathlib import Path

from .core import check_index, file_exists


def run(bcf: Path, contigs: list[str] | None = None, cfg: dict | None = None,
        out: Path | None = None, logs: Path | None = None,
        resume: bool = False) -> Path:
    """Проверить вход и явно сообщить, что режим ещё не разработан."""
    path = file_exists(bcf)
    if not (str(path).endswith('.bcf') or str(path).endswith('.vcf.gz')):
        raise ValueError('Ожидается индексированный .bcf или .vcf.gz')
    check_index(path)
    raise NotImplementedError(
        'VCF/BCF → QC → population structure → population_labels.tsv пока заглушка; '
        'PCA/кластеризация не выполнялись и labels.tsv не создан.'
    )
