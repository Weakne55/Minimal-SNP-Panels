"""Проверки FASTQ-манифеста без генерации геномных данных."""
# MIN_SNP_FASTQ_VCF_PATCH_V1
from pathlib import Path

import pytest

from minsnp import wgs
from minsnp.cli import parser


def _sheet(path: Path, rows: list[str]) -> None:
    path.write_text(
        "sample_id\tspecies\trun_id\tlibrary\tfastq_r1\tfastq_r2\n"
        + "\n".join(rows) + "\n"
    )


def test_fastq_paths_are_resolved_relative_to_manifest(tmp_path, monkeypatch):
    # Проверяем контракт путей, без создания FASTQ/генотипов.
    monkeypatch.setattr(wgs, "file_exists", lambda p: Path(p).resolve())
    sheet = tmp_path / "samples.tsv"
    _sheet(sheet, ["S1\tbeluga\tRUN1\tLIB1\treads/R1.fastq.gz\treads/R2.fastq.gz"])
    rows = wgs.validate_fastq(str(sheet), "beluga")
    assert rows["S1"][0]["fastq_r1"] == str(tmp_path / "reads/R1.fastq.gz")
    assert rows["S1"][0]["fastq_r2"] == str(tmp_path / "reads/R2.fastq.gz")


def test_fastq_reuse_is_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr(wgs, "file_exists", lambda p: Path(p).resolve())
    sheet = tmp_path / "samples.tsv"
    _sheet(sheet, [
        "S1\tbeluga\tRUN1\tLIB1\ta.fastq.gz\tb.fastq.gz",
        "S2\tbeluga\tRUN2\tLIB1\ta.fastq.gz\tc.fastq.gz",
    ])
    with pytest.raises(ValueError, match="FASTQ используется более одного раза"):
        wgs.validate_fastq(str(sheet), "beluga")


def test_fastq_duplicate_unit_is_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr(wgs, "file_exists", lambda p: Path(p).resolve())
    sheet = tmp_path / "samples.tsv"
    _sheet(sheet, [
        "S1\tbeluga\tRUN1\tLIB1\ta.fastq.gz\tb.fastq.gz",
        "S1\tbeluga\tRUN1\tLIB1\tc.fastq.gz\td.fastq.gz",
    ])
    with pytest.raises(ValueError, match="Повторная строка"):
        wgs.validate_fastq(str(sheet), "beluga")


def test_cli_has_fastq_to_joint_vcf_command():
    args = parser().parse_args(["wgs-vcf", "--fastq", "samples.tsv", "--species", "beluga"])
    assert args.command == "wgs-vcf"
