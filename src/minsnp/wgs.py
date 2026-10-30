"""FASTQ (несколько runs на sample) → BAM → хромосомные gVCF.

Usage: вызывается из CLI командой ``run_pipeline.sh wgs-vcf --fastq ...``.
Вход: TSV sample_id, species, run_id, library, fastq_r1, fastq_r2.
Выход: BAM каждого sample, QC и gvcf_manifest.tsv.
"""
# MIN_SNP_FASTQ_VCF_PATCH_V1
from __future__ import annotations

import subprocess
from collections import defaultdict
from pathlib import Path

from .core import check_index, file_exists, identifier, read_tsv, stage, write_tsv


def validate_fastq(sheet: str, species: str) -> dict[str, list[dict[str, str]]]:
    """Проверить таблицу и раскрыть относительные FASTQ-пути от директории TSV."""
    path = Path(sheet).expanduser().resolve()
    rows = read_tsv(path, ("sample_id", "species", "run_id", "library", "fastq_r1", "fastq_r2"))
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    units: set[tuple[str, str, str]] = set()
    reads: set[Path] = set()

    for row in rows:
        sample, run_id, library = (identifier(row[key]) for key in ("sample_id", "run_id", "library"))
        if row["species"] != species:
            raise ValueError(f"Неверный вид в FASTQ TSV: {row['species']!r}, ожидался {species!r}")
        unit = (sample, run_id, library)
        if unit in units:
            raise ValueError(f"Повторная строка sample/run/library: {unit}")
        units.add(unit)

        resolved = dict(row)
        for key in ("fastq_r1", "fastq_r2"):
            filename = Path(row[key]).expanduser()
            if not filename.is_absolute():
                filename = path.parent / filename
            filename = file_exists(filename)
            if filename in reads:
                raise ValueError(f"FASTQ используется более одного раза: {filename}")
            reads.add(filename)
            resolved[key] = str(filename)
        grouped[sample].append(resolved)
    return dict(grouped)


def _valid_bam(path: Path, index: bool = False) -> bool:
    """Проверка целостности BAM перед пропуском этапа при --resume."""
    if not path.is_file() or path.stat().st_size == 0:
        return False
    if index and not any(Path(str(path) + suffix).is_file() for suffix in (".bai", ".csi")):
        return False
    return subprocess.run(["samtools", "quickcheck", "-v", str(path)], capture_output=True).returncode == 0


def _valid_gvcf(path: Path) -> bool:
    """Проверить индекс и доступность хромосомного gVCF (не full-content QC)."""
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        check_index(path)
    except FileNotFoundError:
        return False
    return subprocess.run(["bcftools", "index", "-n", str(path)], capture_output=True).returncode == 0


# MIN_SNP_BWA_MEM3_PREFLIGHT_V1
def _require_bwa_mem3_index(reference_fasta: str) -> None:
    """Проверить индекс BWA-MEM3 до запуска трудоёмкого картирования."""
    missing = [
        f"{reference_fasta}{suffix}"
        for suffix in (".amb", ".ann", ".pac", ".bwt.2bit.64")
        if not Path(f"{reference_fasta}{suffix}").is_file()
        or Path(f"{reference_fasta}{suffix}").stat().st_size == 0
    ]
    if missing:
        raise ValueError(
            "Референс не проиндексирован для BWA-MEM3: "
            + ", ".join(missing)
            + f". Постройте индекс: bwa-mem3 index {reference_fasta!r}"
        )


def run(sheet: str, species: str, spec: dict, cfg: dict, out: Path, logs: Path, resume: bool) -> Path:
    """Подготовить BAM и per-contig gVCF; вернуть manifest для joint calling."""
    grouped = validate_fastq(sheet, species)
    threads = int(cfg.get("threads", 8))
    if threads < 2:
        raise ValueError("Для WGS требуется threads >= 2")
    fasta = spec["fasta"]
    _require_bwa_mem3_index(fasta)
    manifest = []

    for sample, rows in sorted(grouped.items()):
        lane_bams: list[Path] = []
        for row in rows:
            # Отдельные каталоги run/library исключают коллизии при похожих ID.
            lane = out / "lanes" / sample / row["run_id"] / row["library"]
            lane.mkdir(parents=True, exist_ok=True)
            r1, r2 = row["fastq_r1"], row["fastq_r2"]

            if cfg.get("trim_fastq", False):
                t1, t2 = lane / "trim_R1.fastq.gz", lane / "trim_R2.fastq.gz"
                report = lane / "trim_R1.json"
                if not resume or not all(p.is_file() and p.stat().st_size for p in (t1, t2, report)):
                    stage("trim.sh", [r1, r2, t1, t2, threads], logs, f"trim_{sample}_{row['run_id']}_{row['library']}")
                r1, r2 = str(t1), str(t2)

            bam = lane / "fixed.bam"
            if not resume or not _valid_bam(bam):
                stage(
                    "align.sh",
                    [fasta, r1, r2, bam, sample, row["library"], row["run_id"], threads],
                    logs, f"align_{sample}_{row['run_id']}_{row['library']}",
                )
            if not _valid_bam(bam):
                raise RuntimeError(f"Не удалось получить корректный BAM: {bam}")
            lane_bams.append(bam)

        merged = out / "bam" / f"{sample}.merged.bam"
        merged.parent.mkdir(parents=True, exist_ok=True)
        if not resume or not _valid_bam(merged):
            stage("merge_bam.sh", [merged, threads, *lane_bams], logs, f"merge_{sample}")
        if not _valid_bam(merged):
            raise RuntimeError(f"Некорректный объединённый BAM: {merged}")

        final = out / "bam" / f"{sample}.bam"
        if not resume or not _valid_bam(final, index=True):
            stage("markdup.sh", [merged, final, out / "qc" / f"{sample}.dup.txt", threads], logs, f"markdup_{sample}")
        if not _valid_bam(final, index=True):
            raise RuntimeError(f"Некорректный markdup BAM: {final}")

        prefix = out / "qc" / sample
        qc_outputs = [
            Path(f"{prefix}.flagstat.txt"),
            Path(f"{prefix}.stats.txt"),
            Path(f"{prefix}.mosdepth.summary.txt"),
        ]
        if not resume or not all(p.is_file() and p.stat().st_size for p in qc_outputs):
            stage("bam_qc.sh", [final, prefix, threads], logs, f"qc_{sample}")

        for contig in spec["contigs"]:
            gvcf = out / "gvcf" / sample / f"{contig}.g.vcf.gz"
            gvcf.parent.mkdir(parents=True, exist_ok=True)
            if not resume or not _valid_gvcf(gvcf):
                stage("haplotypecaller.sh", [fasta, final, contig, gvcf, threads], logs, f"hc_{sample}_{contig}")
            if not _valid_gvcf(gvcf):
                raise RuntimeError(f"Некорректный gVCF или индекс: {gvcf}")
            manifest.append((sample, contig, str(gvcf.resolve())))

    output = out / "gvcf_manifest.tsv"
    write_tsv(output, ["sample_id", "contig", "gvcf_path"], manifest)
    return output
