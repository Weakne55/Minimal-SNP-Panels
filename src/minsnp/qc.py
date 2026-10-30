"""Техническая фильтрация, без использования меток популяций."""
from __future__ import annotations
from pathlib import Path
import subprocess
from .core import check_index, file_exists, samples_in_vcf, stage, write_tsv


def run(vcf: Path, fasta: str, contigs: list[str], cfg: dict, out: Path, logs: Path, resume: bool) -> Path:
    file_exists(vcf)
    check_index(vcf)
    if not samples_in_vcf(vcf):
        raise ValueError("Во входном VCF/BCF отсутствуют образцы")
    out.mkdir(parents=True, exist_ok=True)
    normalized = out / "normalized.bcf"
    filtered = out / "filtered.bcf"
    if not resume or not filtered.is_file():
        q = cfg["qc"]
        indexed = subprocess.run(["bcftools", "index", "-s", str(vcf)], check=True, capture_output=True, text=True)
        present = {line.split("\t")[0] for line in indexed.stdout.splitlines() if line.strip()}
        selected = [c for c in contigs if c in present]
        if not selected:
            raise ValueError("В VCF нет аутосом, указанных для выбранной сборки")
        stage("qc.sh", [fasta, vcf, normalized, filtered, q["min_dp"], q["min_gq"], q["max_missing"], q["min_maf"], ",".join(selected)], logs, "variant_qc")
    if not resume or not (out / "bcftools.stats.txt").is_file():
        with (out / "bcftools.stats.txt").open("w") as fh:
            subprocess.run(["bcftools", "stats", "-s", "-", str(filtered)], stdout=fh, check=True)
    # PSC: per-sample статистика, полученная из настоящего VCF.
    sample_rows = []
    for line in (out / "bcftools.stats.txt").read_text().splitlines():
        if line.startswith("PSC\t"):
            parts = line.split("\t")
            if len(parts) >= 7:
                sample_rows.append((parts[2], parts[3], parts[4], parts[5], parts[6]))
    if sample_rows:
        write_tsv(out / "sample_qc.tsv", ["sample_id", "nRefHom", "nNonRefHom", "nHets", "nTransitions"], sample_rows)
    result = subprocess.run(["bcftools", "index", "-n", str(filtered)], check=True, capture_output=True, text=True)
    if int(result.stdout.strip()) == 0:
        raise ValueError("После технического QC не осталось SNP")
    write_tsv(out / "qc_summary.tsv", ["metric", "value"], [("variants", result.stdout.strip()), ("samples", len(samples_in_vcf(filtered)))])
    check_index(filtered)
    return filtered
