"""BAM/CRAM header validation without external samtools or actual WGS."""

from types import SimpleNamespace
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import validate_alignment


@pytest.fixture
def reference(tmp_path, monkeypatch):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\nACGT\n>chr2\nTTA\n")
    (tmp_path / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\nchr2\t3\t17\t3\t4\n")
    (tmp_path / "ref.dict").write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:4\tM5:abc\n@SQ\tSN:chr2\tLN:3\n"
    )
    def fake_run(args, **kwargs):
        assert args == ["samtools", "view", "-H", "-T", str(fasta), "reads.bam"]
        assert kwargs == {"check": True, "capture_output": True, "text": True}
        return SimpleNamespace(stdout=fake_run.header)

    fake_run.header = (
        "@HD\tVN:1.6\tSO:coordinate\n"
        "@SQ\tSN:chr1\tLN:4\tM5:abc\n"
        "@SQ\tSN:chr2\tLN:3\n"
        "@RG\tID:lane1\tSM:S1\n"
    )
    monkeypatch.setattr(validate_alignment.subprocess, "run", fake_run)
    return fasta, fake_run


def test_valid_header(reference):
    fasta, stub = reference
    assert validate_alignment.validate("reads.bam", "S1", str(fasta)) is True
    stub.header += "@RG\tID:lane2\tSM:S1\n"
    assert validate_alignment.validate("reads.bam", "S1", str(fasta)) is True


@pytest.mark.parametrize("old,new,error", [
    ("SO:coordinate", "SO:queryname", "SO:coordinate"),
    ("SM:S1", "SM:OTHER", "@RG SM"),
    ("@RG\tID:lane1\tSM:S1\n", "", "@RG SM"),
    ("SN:chr2", "SN:unknown", "contig set"),
    ("LN:4", "LN:7", "length mismatch"),
    ("M5:abc", "M5:def", "MD5 differs"),
    ("@SQ\tSN:chr2\tLN:3\n", "@SQ\tSN:chr2\tLN:3\n@SQ\tSN:chr2\tLN:3\n", "Duplicate @SQ"),
])
def test_invalid_headers(reference, old, new, error):
    fasta, stub = reference
    stub.header = stub.header.replace(old, new)
    with pytest.raises(ValueError, match=error):
        validate_alignment.validate("reads.bam", "S1", str(fasta))


def test_missing_optional_md5_is_accepted(reference):
    fasta, stub = reference
    stub.header = stub.header.replace("\tM5:abc", "")
    assert validate_alignment.validate("reads.bam", "S1", str(fasta)) is True
