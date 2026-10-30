"""LD pruning, PCA и предварительные генетические кластеры."""
from __future__ import annotations
from pathlib import Path
import csv
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from .core import call, samples_in_vcf, write_tsv


def run(bcf: Path, contigs: list[str], cfg: dict, out: Path, logs: Path, resume: bool) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    c = cfg["structure"]
    prefix = out / "genotypes"
    prune = out / "prune"
    numeric_bcf = out / "plink_input.bcf"
    chromosome_map = out / "chromosome_map.tsv"
    # PLINK не рассматривает CM... как аутосомы. Перенумеровываем только временную копию.
    from .core import write_tsv
    write_tsv(chromosome_map, ["original_contig", "plink_chromosome"], [(contig, i) for i, contig in enumerate(contigs, 1)])
    rename = out / "rename_chrs.tsv"
    rename.write_text("".join(f"{contig}\t{i}\n" for i, contig in enumerate(contigs, 1)))
    call(["bcftools", "annotate", "--rename-chrs", str(rename), "-Ob", "-o", str(numeric_bcf), str(bcf)], logs, "rename_contigs")
    call(["bcftools", "index", "-f", str(numeric_bcf)], logs, "index_plink_input")
    chrom_set = ["--chr-set", str(len(contigs))]
    pca = out / "pca"
    labels_file = out / "population_labels.tsv"
    if resume and labels_file.is_file():
        from .reporting import plot_structure
    plot_structure(out, rows, pcs, best[1])
    return labels_file
    call(["plink2", "--bcf", str(numeric_bcf), "--double-id", *chrom_set, "--set-all-var-ids", "@:#:$r:$a", "--make-bed", "--out", str(prefix)], logs, "plink_import")
    # При kb-окне PLINK2 допускает шаг только в один SNP.
    if int(c["ld_step"]) != 1:
        raise ValueError("PLINK2 с kb-window требует ld_step=1")
    call(["plink2", "--bfile", str(prefix), *chrom_set, "--indep-pairwise", f"{int(c['ld_window_kb'])}kb", str(c["ld_r2"]), "--out", str(prune)], logs, "ld_prune")
    # PCA не использует SNP, исключённые LD pruning.
    n_pcs = min(10, len(samples_in_vcf(bcf)) - 1)
    if n_pcs < 2:
        raise ValueError("Для PCA нужно минимум три образца")
    call(["plink2", "--bfile", str(prefix), *chrom_set, "--extract", str(prune) + ".prune.in", "--pca", str(n_pcs), "--out", str(pca)], logs, "pca")
    with Path(str(pca) + ".eigenvec").open() as fh:
        raw = [line.split() for line in fh if line.strip()]
    header = raw[0]
    rows = [dict(zip(header, values)) for values in raw[1:]]
    if len(rows) < 4:
        raise ValueError("Для кластеризации нужно минимум четыре образца")
    pcs = sorted((col for col in rows[0] if col.startswith("PC")), key=lambda x: int(x[2:]))
    if len(pcs) < 2:
        raise ValueError("Недостаточно компонент PCA")
    x = np.array([[float(row[key]) for key in pcs[:min(10, len(pcs))]] for row in rows])
    best = None
    metrics = []
    for k in range(2, min(int(c["max_clusters"]), len(rows) - 1) + 1):
        groups = KMeans(n_clusters=k, random_state=int(c["random_state"]), n_init=20).fit_predict(x)
        if len(set(groups)) < 2 or len(set(groups)) >= len(rows):
            continue
        score = silhouette_score(x, groups)
        metrics.append((k, round(float(score), 6)))
        if best is None or score > best[0]:
            best = (score, groups)
    if best is None:
        raise ValueError("Не удалось определить устойчивый набор кластеров")
    write_tsv(out / "cluster_diagnostics.tsv", ["k", "silhouette"], metrics)
    write_tsv(out / "pca.tsv", ["sample_id", *pcs], [(row["IID"], *(row[col] for col in pcs)) for row in rows])
    write_tsv(labels_file, ["sample_id", "population"], [(row["IID"], f"POP_{int(group) + 1}") for row, group in zip(rows, best[1])])
    (out / "IMPORTANT.txt").write_text("Кластеры предварительные и определены на тех же SNP. Подтвердите родство, batch effects и географическую интерпретацию перед внешним assignment.\n", encoding="utf-8")
    from .reporting import plot_structure
    plot_structure(out, rows, pcs, best[1])
    return labels_file
