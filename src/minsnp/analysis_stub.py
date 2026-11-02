"""Явные контракты будущих VCF workflows без фиктивных научных результатов.

Вызов stub *не* создаёт labels.tsv, SNP-панели или пустые файлы результата.
Входы только проверяются, после чего выполнение завершается с кодом 2.
"""
from __future__ import annotations

from pathlib import Path

from .core import check_index, file_exists, read_tsv


def pending_analysis(vcf_path: str, labels_path: str | None = None) -> None:
    """Проверить существование входов и сообщить, какой модуль предстоит реализовать."""
    vcf = file_exists(vcf_path)
    if not (str(vcf).endswith(".bcf") or str(vcf).endswith(".vcf.gz")):
        raise ValueError("Ожидается индексированный .bcf или .vcf.gz")
    check_index(vcf)
    if labels_path is None:
        raise NotImplementedError(
            "VCF/BCF → QC → population structure → population_labels.tsv. "
            "Это заглушка; labels.tsv не создан. Сейчас поддерживается полноценный FASTQ → VCF."
        )
    labels = file_exists(labels_path)
    rows = read_tsv(labels, ("sample_id", "population"))
    ids = [row["sample_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError(f"Повторяющиеся sample_id в {labels}")
    raise NotImplementedError(
        "VCF/BCF + population_labels.tsv → QC → SNP ranking/LD → CV → minimal SNP panel. "
        "Критерий минимальности и метрика будут определены отдельно. "
        "Это заглушка; SNP-панель и метрики не созданы."
    )
