"""GenomicsDBImport и GenotypeGVCFs по выбранным хромосомам."""
from __future__ import annotations
from pathlib import Path
from collections import defaultdict
import subprocess
import json
from .core import check_index, file_exists, identifier, read_tsv, samples_in_vcf, stage


def validate_gvcf(sheet: Path, contigs: list[str], fasta: str) -> dict[str, list[tuple[str, str]]]:
    rows = read_tsv(sheet, ("sample_id", "contig", "gvcf_path"))
    fai_lengths = {}
    for line in Path(fasta + ".fai").read_text().splitlines():
        chrom, size, *_ = line.split("\t")
        fai_lengths[chrom] = int(size)
    by_contig = defaultdict(list)
    observed = set()
    all_samples = set()
    for row in rows:
        sample = identifier(row["sample_id"])
        contig = row["contig"]
        if contig not in contigs or (sample, contig) in observed:
            raise ValueError(f"Неожиданная или повторная пара sample × contig: {sample}, {contig}")
        observed.add((sample, contig))
        all_samples.add(sample)
        path = file_exists(row["gvcf_path"])
        check_index(path)
        indexed = subprocess.run(["bcftools", "index", "-s", str(path)], check=True, capture_output=True, text=True)
        sizes = {parts[0]: parts[1] for line in indexed.stdout.splitlines() if (parts := line.split("\t")) and len(parts) >= 2}
        if contig not in sizes:
            raise ValueError(f"gVCF не содержит ожидаемый contig: {path}, {contig}")
        if sizes[contig] != "." and int(sizes[contig]) != fai_lengths[contig]:
            raise ValueError(f"Длина contig в gVCF не совпадает с FASTA: {contig}, {path}")
        actual = samples_in_vcf(path)
        if actual != [sample]:
            raise ValueError(f"Неверный sample ID в {path}: {actual}, ожидался {sample}")
        by_contig[contig].append((sample, str(path)))
    if set(by_contig) != set(contigs):
        raise ValueError("Для некоторых contigs gVCF отсутствуют")
    for contig in contigs:
        if {s for s, _ in by_contig[contig]} != all_samples:
            raise ValueError(f"Неполная матрица sample × contig: {contig}")
    return dict(by_contig)


def run(sheet: Path, spec: dict, out: Path, logs: Path, resume: bool) -> Path:
    mapping = validate_gvcf(sheet, spec["contigs"], spec["fasta"])
    out.mkdir(parents=True, exist_ok=True)
    signature = {"assembly": spec["assembly"], "fasta": spec["fasta"], "entries": []}
    for contig in spec["contigs"]:
        for sample, filename in sorted(mapping[contig]):
            st = Path(filename).stat()
            signature["entries"].append([sample, contig, filename, st.st_size, st.st_mtime_ns])
    lock = out / "cohort.lock.json"
    if lock.exists():
        if json.loads(lock.read_text()) != signature:
            raise ValueError(f"Входы существующей когорты изменились: {lock}")
    else:
        lock.write_text(json.dumps(signature, indent=2) + "\n")
    contig_vcfs = []
    for contig in spec["contigs"]:
        samplemap = out / "sample_maps" / f"{contig}.tsv"
        # Для GATK sample map не содержит заголовка.
        samplemap.parent.mkdir(parents=True, exist_ok=True)
        text = "".join(f"{s}\t{p}\n" for s, p in sorted(mapping[contig]))
        if samplemap.exists() and samplemap.read_text() != text:
            raise ValueError(f"Запрещено менять состав существующей когорты: {samplemap}")
        samplemap.write_text(text)
        db = out / "genomicsdb" / contig
        vcf = out / "per_contig" / f"{contig}.vcf.gz"
        vcf.parent.mkdir(parents=True, exist_ok=True)
        if vcf.is_file() and not (Path(str(vcf) + ".tbi").is_file() or Path(str(vcf) + ".csi").is_file()):
            raise ValueError(f"VCF существует, но индекс отсутствует: {vcf}")
        if not resume or not vcf.is_file():
            if db.exists():
                if not resume or not (db / ".completed").is_file():
                    raise ValueError(f"Незавершённый/чужой GenomicsDB: {db}. Используйте новый run-id")
            else:
                stage("genomicsdb.sh", [spec["fasta"], contig, samplemap, db], logs, f"import_{contig}")
            stage("genotype.sh", [spec["fasta"], db, vcf], logs, f"genotype_{contig}")
        contig_vcfs.append(vcf)
    result = out / "cohort.raw.bcf"
    if result.exists() and not Path(str(result) + ".csi").is_file():
        raise ValueError(f"BCF существует, но индекс отсутствует: {result}")
    if not resume or not result.exists():
        stage("concat.sh", [result, *contig_vcfs], logs, "concat")
    check_index(result)
    expected = {sample for items in mapping.values() for sample, _ in items}
    observed = samples_in_vcf(result)
    if len(observed) != len(expected) or set(observed) != expected:
        raise ValueError(f"Состав образцов в итоговом VCF не совпадает с manifest: {observed}")
    return result
