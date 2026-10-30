"""Графики на основе реальных результатов анализов."""
from collections import defaultdict
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import ConfusionMatrixDisplay, confusion_matrix


def plot_structure(out: Path, rows: list[dict], pcs: list[str], labels: np.ndarray) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))
    for group in sorted(set(labels)):
        idx = [i for i, value in enumerate(labels) if value == group]
        ax.scatter([float(rows[i][pcs[0]]) for i in idx], [float(rows[i][pcs[1]]) for i in idx], label=f"POP_{group+1}")
    ax.set(xlabel=pcs[0], ylabel=pcs[1], title="PCA: генетические кластеры")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "pca.png", dpi=180)
    plt.close(fig)


def plot_panels(out: Path, metrics: list[tuple], predictions: list[tuple]) -> None:
    scores = defaultdict(list)
    for fold, size, actual, bacc, macro in metrics:
        scores[size].append(bacc)
    fig, ax = plt.subplots(figsize=(7, 5))
    sizes = sorted(scores)
    ax.errorbar(sizes, [np.mean(scores[n]) for n in sizes], yerr=[np.std(scores[n]) for n in sizes], marker="o")
    ax.set(xlabel="Число SNP", ylabel="Balanced accuracy", ylim=(0, 1.02), title="Фиксированные размеры панелей: outer CV")
    fig.tight_layout()
    fig.savefig(out / "accuracy_vs_panel_size.png", dpi=180)
    plt.close(fig)
    for size in sizes:
        subset = [(real, pred) for fold, n, sample, real, pred in predictions if n == size]
        labels = sorted({a for a, _ in subset} | {b for _, b in subset})
        cm = confusion_matrix([a for a, _ in subset], [b for _, b in subset], labels=labels)
        fig, ax = plt.subplots(figsize=(6, 5))
        ConfusionMatrixDisplay(cm, display_labels=labels).plot(ax=ax, cmap="Blues", colorbar=False)
        ax.set_title(f"SNP panel: {size} SNP")
        fig.tight_layout()
        fig.savefig(out / f"confusion_matrix_{size}.png", dpi=180)
        plt.close(fig)
