from pathlib import Path

import matplotlib
import numpy as np
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402 - select the noninteractive backend before pyplot


def qq_values(pvalues):
    values = np.asarray(pvalues, dtype=float)
    values = np.sort(values[np.isfinite(values)])
    if ((values < 0) | (values > 1)).any():
        raise ValueError("Invalid p-value for QQ plot")
    expected = -np.log10((np.arange(1, len(values) + 1) - 0.5) / len(values)) if len(values) else np.array([])
    return expected, -np.log10(np.clip(values, np.finfo(float).tiny, 1))


def make_plots(frame, directory, alpha, planted_id):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white"})
    data = frame[np.isfinite(frame["logistic_p"])].copy()
    fig, ax = plt.subplots(figsize=(11, 4.8), constrained_layout=True)
    offset, ticks, labels = 0, [], []
    def order(chrom):
        return int(chrom) if str(chrom).isdigit() else {"X": 23, "Y": 24, "MT": 25}.get(str(chrom), 26)
    for index, chrom in enumerate(sorted(data["chrom"].astype(str).unique(), key=order)):
        group = data[data["chrom"].astype(str) == chrom].sort_values("pos")
        positions = group["pos"].to_numpy() + offset
        ax.scatter(positions, -np.log10(np.clip(group["logistic_p"], np.finfo(float).tiny, 1)), s=15, color=["#27637f", "#5a9f9a"][index % 2])
        ticks.append(float(positions.mean()))
        labels.append(chrom)
        planted = group[group["variant_id"] == planted_id]
        if not planted.empty:
            x = float(planted.iloc[0]["pos"] + offset)
            y = -np.log10(max(float(planted.iloc[0]["logistic_p"]), np.finfo(float).tiny))
            ax.scatter([x], [y], color="#bf3f55", s=45, zorder=4)
            ax.annotate(planted_id + " (planted)", (x, y), xytext=(7, -12), textcoords="offset points", fontsize=9)
        offset = float(positions.max()) + 1e6
    ax.axhline(-np.log10(alpha / len(frame)), linestyle="--", color="#bf3f55", linewidth=1, label=f"Bonferroni alpha={alpha}, m={len(frame)}")
    ax.set(xticks=ticks, xticklabels=labels, xlabel="Chromosome (sparse synthetic markers)", ylabel="-log10 logistic Wald p-value", title="Synthetic association scan")
    ax.legend(loc="upper right", frameon=False, fontsize=8)
    fig.savefig(directory / "manhattan.png", dpi=150)
    plt.close(fig)
    expected, observed = qq_values(data["logistic_p"])
    fig, ax = plt.subplots(figsize=(5.5, 5), constrained_layout=True)
    ax.scatter(expected, observed, s=15, color="#27637f")
    limit = max(float(expected.max()) if len(expected) else 1, 1)
    ax.plot([0, limit], [0, limit], "--", color="#8998a8", linewidth=1)
    ax.set(xlabel="Expected -log10 p-value", ylabel="Observed -log10 p-value", title="QQ plot: synthetic simulation")
    fig.savefig(directory / "qq.png", dpi=150)
    plt.close(fig)
