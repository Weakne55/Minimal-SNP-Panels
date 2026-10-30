"""Отбор информативных SNP внутри folds и оценка классификации."""
from __future__ import annotations
import csv
import json
from collections import Counter
from pathlib import Path
import subprocess
import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold
from sklearn.metrics import balanced_accuracy_score, f1_score, accuracy_score, confusion_matrix
import joblib
from .core import read_tsv, samples_in_vcf, write_tsv


def extract_genotypes(vcf: Path, samples: list[str], out: Path) -> tuple[np.memmap, list[tuple[str, int, str, str]]]:
    """Потоковая запись SNP × sample в бинарный файл, без загрузки всего VCF в RAM."""
    out.mkdir(parents=True, exist_ok=True)
    raw = out / "genotypes.i8"
    sites = []
    fmt = r"%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n"
    cmd = ["bcftools", "query", "-f", fmt, str(vcf)]
    with raw.open("wb") as stream:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        assert proc.stdout is not None
        for line in proc.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 4 + len(samples):
                proc.kill()
                raise ValueError("Несовпадение числа sample columns в bcftools query")
            chrom, pos, ref, alt = fields[:4]
            # Только диплоидные двухаллельные SNP; остальные игнорируются.
            if len(ref) != 1 or len(alt) != 1 or "," in alt:
                continue
            dosages = []
            valid = True
            for gt in fields[4:]:
                alleles = gt.replace("|", "/").split("/")
                if len(alleles) != 2:
                    valid = False
                    break
                if "." in alleles:
                    dosages.append(-1)
                elif all(a in {"0", "1"} for a in alleles):
                    dosages.append(int(alleles[0]) + int(alleles[1]))
                else:
                    valid = False
                    break
            if not valid:
                continue
            stream.write(np.asarray(dosages, dtype=np.int8).tobytes())
            sites.append((chrom, int(pos), ref, alt))
        error = proc.stderr.read() if proc.stderr is not None else ""
        if proc.wait():
            raise RuntimeError("bcftools query: " + error[-1000:])
    if not sites:
        raise ValueError("После отбора диплоидных biallelic SNP нет данных")
    write_tsv(out / "variants.tsv", ["contig", "pos", "ref", "alt"], sites)
    mat = np.memmap(raw, mode="r", dtype="i1", shape=(len(sites), len(samples)))
    return mat, sites


def ranking(mat: np.ndarray, y: np.ndarray, train: np.ndarray, batch: int = 4096) -> tuple[np.ndarray, np.ndarray]:
    """FST-подобная мера (Ht-Hs)/Ht, вычисленная только по train."""
    classes = np.unique(y[train])
    scores = np.zeros(mat.shape[0], dtype=np.float32)
    for start in range(0, mat.shape[0], batch):
        stop = min(start + batch, mat.shape[0])
        x = np.asarray(mat[start:stop, train], dtype=np.float32)
        x[x < 0] = np.nan
        with np.errstate(invalid="ignore", divide="ignore"):
            p = np.nanmean(x, axis=1) / 2
            ht = 2 * p * (1 - p)
            hs = np.zeros(stop - start, dtype=np.float32)
            den = np.zeros(stop - start, dtype=np.float32)
            for c in classes:
                group = x[:, y[train] == c]
                cnt = np.isfinite(group).sum(axis=1)
                pg = np.nanmean(group, axis=1) / 2
                hg = 2 * pg * (1 - pg)
                hs += np.nan_to_num(hg) * cnt
                den += cnt
            hs /= np.maximum(den, 1)
            scores[start:stop] = np.maximum(0, np.nan_to_num((ht - hs) / ht, nan=0, posinf=0, neginf=0))
    return np.argsort(-scores, kind="stable"), scores


def ld_select(mat: np.ndarray, ranked: np.ndarray, train: np.ndarray, size: int, r2: float, pool: int) -> list[int]:
    """LD контроль по train; отсутствующие GT заменяются train средним только для LD."""
    chosen, vectors = [], []
    for idx in ranked[:pool]:
        vals = np.asarray(mat[int(idx), train], dtype=float)
        vals[vals < 0] = np.nan
        if np.isfinite(vals).sum() < 2:
            continue
        vals[np.isnan(vals)] = np.nanmean(vals)
        vals -= vals.mean()
        norm = np.linalg.norm(vals)
        if norm < 1e-8:
            continue
        vals /= norm
        if any(float(np.dot(vals, other)) ** 2 > r2 for other in vectors):
            continue
        chosen.append(int(idx))
        vectors.append(vals)
        if len(chosen) == size:
            break
    return chosen


def train_predict(mat: np.ndarray, y: np.ndarray, train: np.ndarray, test: np.ndarray, selected: list[int]) -> tuple[np.ndarray, object]:
    if not selected:
        raise ValueError("Не удалось выбрать SNP: проверьте MAF, sample size и LD")
    x_train = np.asarray(mat[selected, :][:, train].T, dtype=float)
    x_test = np.asarray(mat[selected, :][:, test].T, dtype=float)
    x_train[x_train < 0] = np.nan
    x_test[x_test < 0] = np.nan
    model = make_pipeline(SimpleImputer(strategy="most_frequent", keep_empty_features=True), LogisticRegression(max_iter=1000, class_weight="balanced"))
    model.fit(x_train, y[train])
    return model.predict(x_test), model


def cv_splits(y: np.ndarray, groups: np.ndarray | None, folds: int, seed: int):
    labels, counts = np.unique(y, return_counts=True)
    k = min(folds, int(counts.min()))
    if k < 2:
        raise ValueError("Минимум 2 образца в каждой популяции для cross-validation")
    if groups is None:
        return list(StratifiedKFold(n_splits=k, shuffle=True, random_state=seed).split(np.zeros(len(y)), y))
    if len(np.unique(groups)) < k:
        raise ValueError("Недостаточно независимых групп для CV")
    pairs = list(StratifiedGroupKFold(n_splits=k, shuffle=True, random_state=seed).split(np.zeros(len(y)), y, groups))
    for tr, te in pairs:
        if set(np.unique(y[tr])) != set(labels) or set(np.unique(y[te])) != set(labels):
            raise ValueError("Group CV не обеспечивает присутствие популяций в folds; укрупните группы или уменьшите folds")
    return pairs


def run(vcf: Path, labels: Path, cfg: dict, out: Path, logs: Path, assembly: str) -> Path:
    samples = samples_in_vcf(vcf)
    metadata = read_tsv(labels, ("sample_id", "population"))
    indexed = {row["sample_id"]: row for row in metadata}
    if len(indexed) != len(metadata) or set(indexed) != set(samples):
        raise ValueError("sample_id в labels должны точно совпадать с заголовком VCF, без дубликатов")
    y = np.asarray([indexed[s]["population"] for s in samples])
    group_values = [indexed[s].get("group", "") for s in samples]
    if any(group_values) and not all(group_values):
        raise ValueError("Столбец group должен быть заполнен у всех образцов или ни у кого")
    groups = np.asarray(group_values) if all(group_values) else None
    conf = cfg["panels"]
    sizes = sorted(set(int(n) for n in conf["sizes"]))
    if not sizes or any(n < 1 for n in sizes):
        raise ValueError("Размеры панелей должны быть положительными")
    folds = cv_splits(y, groups, int(conf["folds"]), int(conf["seed"]))
    mat, sites = extract_genotypes(vcf, samples, out / "features")
    all_metrics, predictions, selections = [], [], []
    scores_by_size = {n: [] for n in sizes}
    for fold, (tr, te) in enumerate(folds, start=1):
        order, scores = ranking(mat, y, tr)
        for size in sizes:
            selected = ld_select(mat, order, tr, size, float(conf["ld_r2"]), int(conf["selection_pool"]))
            pred, _ = train_predict(mat, y, tr, te, selected)
            bacc = balanced_accuracy_score(y[te], pred)
            f1 = f1_score(y[te], pred, average="macro", zero_division=0)
            all_metrics.append((fold, size, len(selected), round(float(bacc), 6), round(float(f1), 6)))
            scores_by_size[size].append(float(bacc))
            selections.extend((fold, size, *sites[idx], float(scores[idx])) for idx in selected)
            predictions.extend((fold, size, samples[i], y[i], pred[j]) for j, i in enumerate(te))
    write_tsv(out / "cv_metrics.tsv", ["fold", "target_size", "actual_size", "balanced_accuracy", "macro_f1"], all_metrics)
    write_tsv(out / "cv_predictions.tsv", ["fold", "target_size", "sample_id", "true", "pred"], predictions)
    write_tsv(out / "fold_selections.tsv", ["fold", "target_size", "contig", "pos", "ref", "alt", "score"], selections)
    write_tsv(out / "panel_size_summary.tsv", ["size", "mean_balanced_accuracy"], [(n, np.mean(v)) for n, v in scores_by_size.items()])
    # Выбор по той же CV может быть оптимистичен: фиксируем это в отчёте.
    # Nested CV: подбор размера панели внутри outer-train, оценка только на outer-test.
    nested = []
    nested_error = None
    for fold, (outer_train, outer_test) in enumerate(folds, start=1):
        try:
            inner_groups = groups[outer_train] if groups is not None else None
            inner = cv_splits(y[outer_train], inner_groups, min(3, int(conf["folds"])), int(conf["seed"]) + fold)
        except ValueError as exc:
            nested_error = str(exc)
            break
        inner_scores = {size: [] for size in sizes}
        for inner_train, inner_test in inner:
            train = outer_train[inner_train]
            test = outer_train[inner_test]
            order, _ = ranking(mat, y, train)
            for size in sizes:
                chosen = ld_select(mat, order, train, size, float(conf["ld_r2"]), int(conf["selection_pool"]))
                pred, _ = train_predict(mat, y, train, test, chosen)
                inner_scores[size].append(float(balanced_accuracy_score(y[test], pred)))
        inner_best = max(np.mean(v) for v in inner_scores.values())
        picked_size = min(size for size, values in inner_scores.items() if np.mean(values) >= inner_best - 0.02)
        ranked, _ = ranking(mat, y, outer_train)
        chosen = ld_select(mat, ranked, outer_train, picked_size, float(conf["ld_r2"]), int(conf["selection_pool"]))
        pred, _ = train_predict(mat, y, outer_train, outer_test, chosen)
        nested.append((fold, picked_size, len(chosen), round(float(balanced_accuracy_score(y[outer_test], pred)), 6), round(float(f1_score(y[outer_test], pred, average="macro", zero_division=0)), 6)))
    if nested_error is None:
        write_tsv(out / "nested_cv_metrics.tsv", ["outer_fold", "inner_selected_size", "actual_size", "balanced_accuracy", "macro_f1"], nested)
    else:
        (out / "nested_cv_unavailable.txt").write_text("Недостаточно данных для nested CV: " + nested_error + "\n")
    best = max(np.mean(v) for v in scores_by_size.values())
    size = min(n for n in sizes if np.mean(scores_by_size[n]) >= best - 0.02)
    final_order, final_scores = ranking(mat, y, np.arange(len(y)))
    final_ids = ld_select(mat, final_order, np.arange(len(y)), size, float(conf["ld_r2"]), int(conf["selection_pool"]))
    _, model = train_predict(mat, y, np.arange(len(y)), np.arange(len(y)), final_ids)
    write_tsv(out / "minimal_panel.tsv", ["assembly_id", "contig", "pos", "ref", "alt", "score"], [(assembly, *sites[i], round(float(final_scores[i]), 8)) for i in final_ids])
    joblib.dump({"model": model, "assembly_id": assembly, "sites": [sites[i] for i in final_ids],
                 "encoding": "ALT allele dosage 0/1/2; -1 missing; most_frequent imputation from train",
                 "samples": samples, "populations": sorted(set(y))}, out / "model.joblib")
    (out / "panel_selection.json").write_text(json.dumps({"assembly_id": assembly, "selection_method": "Nei-style (Ht-Hs)/Ht + train-fold LD r2", "cv_balanced_accuracy_by_size": {str(n): round(float(np.mean(scores)), 6) for n, scores in scores_by_size.items()}, "tolerance": 0.02, "chosen_target_size": size, "actual_size": len(final_ids)}, ensure_ascii=False, indent=2) + "\n")
    (out / "LIMITATIONS.txt").write_text("CV по фиксированным размерам выполняется с отбором SNP внутри folds. При наличии nested_cv_metrics.tsv это независимая outer-fold оценка процедуры выбора размера; без него качество подобранного размера требует внешней проверки. Финальный размер выбран по CV на полной когорте и не является внешне валидированным. Если labels получены кластеризацией этих же WGS, оценивается воспроизведение кластеров, а не независимая география.\n", encoding="utf-8")
    from .reporting import plot_panels
    plot_panels(out, all_metrics, predictions)
    return out / "minimal_panel.tsv"
