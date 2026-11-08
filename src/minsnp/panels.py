"""Заглушка VCF/BCF + metadata → SNP panel."""
from pathlib import Path


def run(vcf: Path, metadata: Path) -> None:
    """Не выполняет анализ и не создаёт результатов."""
    raise NotImplementedError(
        "VCF/BCF + metadata → SNP panel пока не реализован."
    )
