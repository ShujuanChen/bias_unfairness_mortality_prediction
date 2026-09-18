#!/usr/bin/env python3
"""The three figures of the weighting phase, for the weight set this run reports.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "common"))
import reporting as rp
from helpers import WEIGHT_METHOD, WEIGHT_METHOD_LABEL
from paths import path_results, path_temp, setting
from tabular import read_table

PHASE = "phase3_weighting"
WINSORISE_PROBS = tuple(setting("phase3", "winsorise_percentiles"))

BEFORE, AFTER = "#B22222", "#2C5D8A"
COHORT, REFERENCE = "#640404", "#2C5D8A"


def _table(name):
    path = path_results(PHASE, name, create=False)
    if not path.exists():
        raise SystemExit(f"No {name} in this run's results. Run the phase 3 "
                         "diagnostics stage first.")
    return read_table(path)


def overlap_figure():
    path = path_temp("weight_probability_density.csv", create=False)
    positivity = path_results(PHASE, "weight_positivity.xlsx", create=False)
    d = read_table(path) if path.exists() else None
    if d is not None:
        d = d[d["weight_source"] == WEIGHT_METHOD]
    bounds = read_table(positivity) if positivity.exists() else None
    if bounds is not None:
        bounds = bounds[bounds["weight_source"] == WEIGHT_METHOD]
    if d is None or not len(d) or bounds is None or not len(bounds):
        print("no fitted probabilities for this weight set, overlap not drawn")
        return
    bounds = bounds.iloc[0]

    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    for sample, colour, label in (("cohort", COHORT, "UKB"),
                                  ("reference", REFERENCE, "Reference sample")):
        side = d[d["sample"] == sample]
        ax.fill_between(side["probability"], side["density"], color=colour,
                        alpha=0.35, linewidth=0)
        ax.plot(side["probability"], side["density"], color=colour,
                linewidth=1.4, label=label)
    for bound in (bounds["overlap_lower"], bounds["overlap_upper"]):
        ax.axvline(float(bound), color="#444444", linestyle="--", linewidth=0.9)
    ax.set_xlim(0, 1)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Fitted probability of cohort membership", fontsize=10)
    ax.set_ylabel("Density", fontsize=10)
    ax.set_title(WEIGHT_METHOD_LABEL, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    fig.tight_layout()
    rp.save_fig(fig, path_results(PHASE), "weight_overlap")


def _density(x, grid, n=None):
    iqr = float(np.subtract(*np.percentile(x, [75, 25])))
    sd = float(x.std(ddof=1))
    scale = min(sd, iqr / 1.34) if iqr > 0 else sd
    bw = 0.9 * scale * len(x) ** (-0.2)
    step = grid[1] - grid[0]
    pad = min(int(np.ceil(4 * bw / step)), (len(grid) - 1) // 2)
    edges = np.concatenate([grid - step / 2, [grid[-1] + step / 2]])
    counts, _ = np.histogram(x, bins=edges)
    n = len(x) if n is None else n
    k = np.arange(-pad, pad + 1) * step
    kernel = np.exp(-0.5 * (k / bw) ** 2)
    kernel /= kernel.sum() * step
    return np.convolve(counts / (n * step), kernel * step, mode="same")


def distribution_figure():
    path = path_temp("weights", WEIGHT_METHOD, "ukb_weights.csv", create=False)
    if not path.exists():
        raise SystemExit(f"No weights at {path}. Run phase 3 first.")
    w = read_table(path)["w"].to_numpy(dtype=float)
    w = w[np.isfinite(w) & (w > 0)]

    hi = 10 ** np.ceil(np.log10(w.max()) + 0.3)
    floor = hi / 1e6
    lo = 10 ** np.floor(np.log10(max(w.min(), floor)) - 0.3)
    inside = (w >= lo) & (w <= hi)
    outside = int((~inside).sum())
    if outside:
        print(f"  {outside} of {len(w)} weights lie outside the drawn range "
              f"{lo:.3g} to {hi:.3g}, smallest {w.min():.3g}, "
              f"largest {w.max():.3g}")
    grid = np.linspace(np.log10(lo), np.log10(hi), 1200)
    d = _density(np.log10(w[inside]), grid, n=len(w))
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    ax.fill_between(10 ** grid, d, color=COHORT, alpha=0.30, linewidth=0)
    ax.plot(10 ** grid, d, color=COHORT, linewidth=1.8)
    cuts = np.quantile(w, WINSORISE_PROBS)
    for bound, prob, upper in zip(cuts, WINSORISE_PROBS, (False, True)):
        if not (lo <= bound <= hi):
            continue
        ax.axvline(bound, color="#444444", linestyle="--", linewidth=0.9)
        ax.annotate(f"{100 * prob:g} pct", xy=(bound, 1.0),
                    xycoords=("data", "axes fraction"),
                    xytext=(3 if upper else -3, -4),
                    textcoords="offset points", fontsize=8, color="#444444",
                    rotation=90, ha="left" if upper else "right", va="top")
    ax.set_xscale("log")
    ax.set_xlim(lo, hi)
    ticks = [10.0 ** e for e in range(int(np.floor(np.log10(lo))),
                                      int(np.ceil(np.log10(hi))) + 1)]
    ticks = [t for t in ticks if lo <= t <= hi]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.set_ylim(0, float(d.max()) * 1.08)
    label = "Normalised participation weights"
    if outside:
        label += (f"\n{outside:,} of {len(w):,} lie outside the axis, "
                  f"smallest {w.min():.2g}")
    ax.set_xlabel(label, fontsize=10)
    ax.set_ylabel("Density per decade", fontsize=10)
    ax.set_title(WEIGHT_METHOD_LABEL, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    rp.save_fig(fig, path_results(PHASE), "weight_distribution")


def main():
    summary = _table("weight_summary.xlsx").set_index("weight_source")
    if WEIGHT_METHOD not in summary.index:
        raise SystemExit(f"weight_summary.xlsx has no row for {WEIGHT_METHOD}.")
    row = summary.loc[WEIGHT_METHOD]
    print(f"{WEIGHT_METHOD_LABEL}: effective sample size "
          f"{float(row['effective_sample_size']):,.0f} of "
          f"{int(row['n_cohort']):,}, "
          f"{100 * float(row['effective_fraction']):.1f} per cent")

    overlap_figure()
    distribution_figure()
    print("phase 3 results written")


if __name__ == "__main__":
    main()
