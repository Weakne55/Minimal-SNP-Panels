"""Контрактные тесты без создания искусственных геномных данных."""
from pathlib import Path
import pytest
from minsnp.core import identifier, read_tsv
from pipeline import build_parser


def test_identifier_security():
    for invalid in ("../bad", "", "A B", "/tmp/path", ".", ".."):
        with pytest.raises(ValueError):
            identifier(invalid)
    assert identifier("run.2026-10-30") == "run.2026-10-30"


def test_manifest_requires_columns(tmp_path):
    path = tmp_path / "fastq.tsv"
    path.write_text("sample_id\tspecies\n", encoding="utf-8")
    with pytest.raises(ValueError, match="столбцы"):
        read_tsv(path, ("sample_id", "species", "fastq_r1"))


def test_empty_manifest_rejected(tmp_path):
    path = tmp_path / "fastq.tsv"
    path.write_text("sample_id\tspecies\n", encoding="utf-8")
    with pytest.raises(ValueError, match="пустая"):
        read_tsv(path, ("sample_id", "species"))


def test_cli_has_two_user_inputs():
    assert build_parser().parse_args(["--fastq", "input.tsv"]).fastq == ["input.tsv"]
    assert build_parser().parse_args(["--vcf", "cohort.bcf"]).vcf == "cohort.bcf"
    with pytest.raises(SystemExit):
        build_parser().parse_args(["--fastq", "a", "--vcf", "b"])


def test_shell_stages_parse():
    import subprocess
    from minsnp.core import ROOT
    scripts = sorted((ROOT / "stages").glob("*.sh"))
    assert len(scripts) >= 8
    subprocess.run(["bash", "-n", *(str(p) for p in scripts)], check=True)


def test_reference_lists_from_old_pipeline():
    import yaml
    from minsnp.core import ROOT
    refs = yaml.safe_load((ROOT / "config/references.yaml").read_text())
    assert len(refs["species"]["beluga"]["contigs"]) == 21
    assert len(refs["species"]["polar_bear"]["contigs"]) == 36
    assert "CM197971.1" not in refs["species"]["polar_bear"]["contigs"]
