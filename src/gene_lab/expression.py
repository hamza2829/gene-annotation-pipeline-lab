import csv
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402 - select the noninteractive backend before pyplot


def summarize_expression(path, symbols, store):
    with open(path, encoding="utf-8") as stream:
        headers = next(csv.reader(stream), [])
    if len(headers) < 2 or len(set(headers)) != len(headers):
        raise ValueError(f"{path}: expression matrix needs unique cell columns")
    frame = pd.read_csv(path, index_col=0)
    if frame.index.duplicated().any() or frame.empty:
        raise ValueError(f"{path}: expression matrix is empty or has duplicate genes")
    try:
        frame = frame.apply(pd.to_numeric, errors="raise")
    except ValueError as error:
        raise ValueError(f"{path}: expression counts must be numeric") from error
    if not np.isfinite(frame.to_numpy()).all() or (frame.to_numpy() < 0).any():
        raise ValueError(f"{path}: expression counts must be finite and nonnegative")
    rows = []
    for symbol in symbols:
        if symbol not in frame.index:
            rows.append({"gene": symbol, "cells": frame.shape[1], "mean": np.nan, "median": np.nan, "detected_fraction": np.nan, "status": "missing"})
        else:
            values = frame.loc[symbol]
            rows.append({"gene": symbol, "cells": len(values), "mean": float(values.mean()), "median": float(values.median()),
                "detected_fraction": float((values > 0).mean()), "status": "ok"})
    result = store.write("expression", rows)
    valid = result[result["status"] == "ok"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    axes[0].bar(valid["gene"], valid["mean"], color="#27637f")
    axes[1].bar(valid["gene"], valid["detected_fraction"] * 100, color="#5a9f9a")
    for ax in axes:
        ax.tick_params(axis="x", rotation=55)
    axes[0].set(title="Toy expression counts", ylabel="Mean synthetic counts per cell")
    axes[1].set(title="Toy detection fraction", ylabel="Cells with nonzero counts (%)", ylim=(0, 100))
    fig.savefig(Path(store.directory) / "expression.png", dpi=150)
    plt.close(fig)
    return rows
