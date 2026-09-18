#!/usr/bin/env python3
"""Result of the correction phase.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.transforms as mtransforms
from matplotlib.gridspec import GridSpec

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "common"))
import consequences as cq
import reporting as rp
from helpers import (ALL_CAUSE, CAUSES, CAUSE_LABELS, WUKB_SOURCE,
                     compute_rates, load_cause, require_weighted_column)
from paths import path_results, path_temp
from survival import HORIZONS_YEARS, risk_column
from tabular import read_table, write_table

BIAS_LABEL = "Prediction bias (deaths per 1,000,000)"
CORRECTION_LABEL = "Bias correction (%)"
UNFLAGGED_LABEL = "Additional unflagged death"

FIG_WIDTH = 8.0
WIDTH_RATIOS = (3.0, 1.2)
WSPACE = 0.08
PAD = 0.4

INTERVAL_TABLES = {"cause": "bootstrap_intervals_by_cause.xlsx",
                   "subgroup": "bootstrap_intervals_by_subgroup.xlsx"}


def load_intervals(shape):
    path = path_results(PHASE5, INTERVAL_TABLES[shape], create=False)
    if not path.exists():
        return None
    frame = read_table(path)
    lo, hi = f"{WUKB_SOURCE}_correction_lo", f"{WUKB_SOURCE}_correction_hi"
    if lo not in frame.columns:
        return None
    keys = (["cause", "horizon_years"] if shape == "cause"
            else ["horizon_years", "variable", "level"])
    out = {}
    for _, r in frame.iterrows():
        key = tuple(str(r[k]) for k in keys)
        out[key] = (float(r[lo]), float(r[hi]))
    return out

ARMS = [("ukb", "unweighted"), (WUKB_SOURCE, "reweighted")]

PHASE = "phase4_bias_correction"
PHASE5 = "phase5_uncertainty"


def correction_pct(ukb_bias, wukb_bias):
    a = abs(ukb_bias)
    return (a - abs(wukb_bias)) / a * 100 if a > 0 else np.nan


def _rates(df, cause, hz):
    pred = {k: risk_column(cause, s, hz) for k, s in
            (("pmr", "pmr"), ("ukb", "ukb"), ("wukb", WUKB_SOURCE))}
    require_weighted_column(df, cause, pred["wukb"])
    r = compute_rates(df, f"event_{hz}y", pred)
    ukb_bias = (r["ukb_pred_rate"] - r["observed_rate"]) * 1e6
    wukb_bias = (r["wukb_pred_rate"] - r["observed_rate"]) * 1e6
    return r, ukb_bias, wukb_bias


# ── 1. Correction by cause ───────────────────────────────────────────────────

def correction_by_cause():
    rows = []
    for cause in CAUSES:
        df = load_cause(cause)
        if df is None:
            print(f"[SKIP] {cause}")
            continue
        for hz in HORIZONS_YEARS:
            r, ukb_bias, wukb_bias = _rates(df, cause, hz)
            rows.append({
                "cause": cause, "cause_label": CAUSE_LABELS[cause],
                "horizon_years": hz, "weight_source": WUKB_SOURCE,
                "n": r["n"], "n_events": r["n_events"],
                "observed_per_mil": r["observed_rate"] * 1e6,
                "pmr_predicted_per_mil": r["pmr_pred_rate"] * 1e6,
                "ukb_predicted_per_mil": r["ukb_pred_rate"] * 1e6,
                "wukb_predicted_per_mil": r["wukb_pred_rate"] * 1e6,
                "pmr_bias_per_mil": (r["pmr_pred_rate"] - r["observed_rate"]) * 1e6,
                "ukb_bias_per_mil": ukb_bias,
                "wukb_bias_per_mil": wukb_bias,
                "ukb_relative_bias": r["ukb_relative_bias"],
                "wukb_relative_bias": r["wukb_relative_bias"],
                "bias_correction_pct": correction_pct(ukb_bias, wukb_bias),
            })
        del df
    out = pd.DataFrame(rows)
    write_table(out, path_temp("bias_correction_by_cause.csv"))

    outdir = path_results(PHASE)
    intervals = load_intervals("cause")
    for hz in HORIZONS_YEARS:
        sub = out[out["horizon_years"] == hz].copy()
        sub["cause"] = pd.Categorical(sub["cause"], categories=CAUSES, ordered=True)
        sub = sub.sort_values("cause").reset_index(drop=True)
        bounds = [intervals.get((str(c), str(hz)), (np.nan, np.nan))
                  if intervals else (np.nan, np.nan) for c in sub["cause"]]
        sub["correction_lo"] = [b[0] for b in bounds]
        sub["correction_hi"] = [b[1] for b in bounds]
        rows = pd.DataFrame({
            "kind": "level",
            "label": [CAUSE_LABELS[c] for c in sub["cause"]],
            "pmr_bias_per_mil": sub["pmr_bias_per_mil"].to_numpy(),
            "ukb_bias_per_mil": sub["ukb_bias_per_mil"].to_numpy(),
            "wukb_bias_per_mil": sub["wukb_bias_per_mil"].to_numpy(),
            "bias_correction_pct": sub["bias_correction_pct"].to_numpy(),
            "correction_lo": sub["correction_lo"].to_numpy(),
            "correction_hi": sub["correction_hi"].to_numpy()})
        fig = draw_pair(rows, show_pmr=True, marker_size=50, bar_height=0.28,
                        cap=0.07, label_fontsize=11, axis_fontsize=12,
                        fig_height=max(3.0, len(rows) * 0.42 + 1.2),
                        y_labels_on_axis=True, xtick_bins=5)
        rp.save_fig(fig, outdir, f"bias_correction_by_cause_{hz}y")


# ── The two-panel figure ─────────────────────────────────────────────────────

def draw_pair(rows, show_pmr, marker_size, bar_height, cap, label_fontsize,
              axis_fontsize, fig_height, y_labels_on_axis, xtick_bins="auto"):
    n = len(rows)
    fig, (ax, ax_bc) = plt.subplots(
        1, 2, figsize=(FIG_WIDTH, fig_height),
        gridspec_kw={"width_ratios": list(WIDTH_RATIOS), "wspace": WSPACE},
        sharey=True)

    ax_bc.set_xlim(0, 100)
    ax_bc.set_xticks([0, 25, 50, 75, 100])

    frame = rows.reset_index(drop=True)
    for i, row in frame.iterrows():
        if row["kind"] != "level":
            continue
        ukb, wukb = row["ukb_bias_per_mil"], row["wukb_bias_per_mil"]
        if show_pmr and np.isfinite(row.get("pmr_bias_per_mil", np.nan)):
            ax.scatter(row["pmr_bias_per_mil"], i, c=rp.C_PMR, marker=rp.M_PMR,
                       s=marker_size, zorder=5, edgecolors="none")
        if np.isfinite(ukb) and np.isfinite(wukb):
            ax.annotate("", xy=(wukb, i), xytext=(ukb, i),
                        arrowprops=dict(arrowstyle="->", color="#888888", lw=0.8))
        if np.isfinite(ukb):
            ax.scatter(ukb, i, c=rp.C_UKB, marker=rp.M_UKB, s=marker_size,
                       zorder=6, edgecolors="none")
        if np.isfinite(wukb):
            ax.scatter(wukb, i, c=rp.C_WEIGHTED, marker=rp.M_WEIGHTED,
                       s=marker_size, zorder=7, edgecolors="none")
        rp.bar_with_interval(ax_bc, i, row["bias_correction_pct"],
                             row.get("correction_lo", np.nan),
                             row.get("correction_hi", np.nan),
                             rp.C_WEIGHTED, height=bar_height, cap=cap)

    y_bot, y_top = n - 0.5 + PAD, -0.5 - PAD
    for a in (ax, ax_bc):
        a.plot([0, 0], [y_bot, y_top], color="grey", linewidth=1.6, alpha=0.7,
               zorder=3)
    ax.invert_yaxis()
    ax.set_ylim(y_bot, y_top)

    levels = [i for i, r in frame.iterrows() if r["kind"] == "level"]
    ax.set_yticks(levels)
    if y_labels_on_axis:
        ax.set_yticklabels([frame.loc[i, "label"] for i in levels],
                           fontsize=label_fontsize)
        ax.tick_params(axis="y", length=4, width=1.2)
    else:
        ax.set_yticklabels([""] * len(levels))
        ax.tick_params(axis="y", length=3.5, width=0.8)
        trans = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
        for i, row in frame.iterrows():
            bold = row["kind"] == "group"
            ax.text(-0.02, i, row["label"] if bold else "  " + row["label"],
                    transform=trans, ha="right", va="center",
                    fontweight="bold" if bold else "normal",
                    fontsize=label_fontsize, clip_on=False)

    ax.xaxis.set_major_locator(
        mticker.MaxNLocator(nbins=xtick_bins, steps=[1, 2, 5, 10]))
    ax.xaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.set_xlabel(BIAS_LABEL, fontsize=axis_fontsize)
    ax_bc.set_xlabel(CORRECTION_LABEL, fontsize=axis_fontsize)
    ax_bc.tick_params(axis="y", left=False, labelleft=False)
    for a in (ax, ax_bc):
        rp.box_spine(a, linewidth=1.0)
        a.tick_params(axis="x", labelsize=axis_fontsize - 2)
    fig.tight_layout()
    return fig


# ── 2. Correction by stratum ─────────────────────────────────────────────────

def correction_by_stratum(df):
    rows = rp.stratum_rows(df)
    outdir = path_results(PHASE)
    intervals = load_intervals("subgroup")
    table = []
    for hz in HORIZONS_YEARS:
        pred = {k: risk_column(ALL_CAUSE, s, hz) for k, s in
                (("pmr", "pmr"), ("ukb", "ukb"), ("wukb", WUKB_SOURCE))}
        require_weighted_column(df, ALL_CAUSE, pred["wukb"])
        drawn = rows.copy()
        for col in ("pmr_bias_per_mil", "ukb_bias_per_mil", "wukb_bias_per_mil",
                    "bias_correction_pct", "correction_lo", "correction_hi"):
            drawn[col] = np.nan
        for i, row in drawn.iterrows():
            if row["kind"] != "level":
                continue
            mask = rp.as_level(df[row["var"]]) == row["level"]
            r = compute_rates(df[mask], f"event_{hz}y", pred)
            ukb_bias = (r["ukb_pred_rate"] - r["observed_rate"]) * 1e6
            wukb_bias = (r["wukb_pred_rate"] - r["observed_rate"]) * 1e6
            drawn.at[i, "pmr_bias_per_mil"] = (r["pmr_pred_rate"] - r["observed_rate"]) * 1e6
            drawn.at[i, "ukb_bias_per_mil"] = ukb_bias
            drawn.at[i, "wukb_bias_per_mil"] = wukb_bias
            drawn.at[i, "bias_correction_pct"] = correction_pct(ukb_bias, wukb_bias)
            if intervals:
                lo, hi = intervals.get((str(hz), row["var"], row["level"]),
                                       (np.nan, np.nan))
                drawn.at[i, "correction_lo"] = lo
                drawn.at[i, "correction_hi"] = hi
            table.append({
                "cause": ALL_CAUSE, "horizon_years": hz,
                "weight_source": WUKB_SOURCE,
                "variable": row["var"], "level": row["level"],
                "n": r["n"], "n_events": r["n_events"],
                "observed_per_mil": r["observed_rate"] * 1e6,
                "pmr_predicted_per_mil": r["pmr_pred_rate"] * 1e6,
                "ukb_predicted_per_mil": r["ukb_pred_rate"] * 1e6,
                "wukb_predicted_per_mil": r["wukb_pred_rate"] * 1e6,
                "ukb_bias_per_mil": ukb_bias, "wukb_bias_per_mil": wukb_bias,
                "bias_correction_pct": drawn.at[i, "bias_correction_pct"]})

        fig = draw_pair(drawn, show_pmr=True, marker_size=45, bar_height=0.45,
                        cap=0.11, label_fontsize=8.5, axis_fontsize=9,
                        fig_height=max(3.0, len(drawn) * 0.22 + 0.5),
                        y_labels_on_axis=False)
        rp.save_fig(fig, outdir,
                    f"bias_correction_by_strata_{hz}y")

    write_table(pd.DataFrame(table),
                path_temp("bias_correction_by_subgroup.csv"))


# ── 3. The unflagged death sweep, before and after ───────────────────────────

def draw_pairs(part, years, outdir, name):
    strata = list(dict.fromkeys(part["stratum"]))
    npair = min(2, len(strata))
    nrow = -(-len(strata) // npair)
    ncol = npair * 2
    W, H = 2.6 * ncol, 3.0 * nrow
    fig = plt.figure(figsize=(W, H))
    outer = GridSpec(nrow, npair, figure=fig, left=1.05 / W, right=1 - 0.14 / W,
                     bottom=0.60 / H, top=1 - 0.52 / H, wspace=0.12, hspace=0.42)

    for i, stratum in enumerate(strata):
        r, c = divmod(i, npair)
        inner = outer[r, c].subgridspec(1, 2, wspace=0.07)
        left = fig.add_subplot(inner[0])
        right = fig.add_subplot(inner[1], sharey=left)
        sub = part[part["stratum"] == stratum]
        for ax, arm, style in ((left, "unweighted", "-"),
                               (right, "reweighted", "-.")):
            for lv in dict.fromkeys(sub["level"]):
                g = sub[sub["level"] == lv].sort_values("threshold")
                ax.plot(g["threshold"] * 100,
                        g[f"additional_unflagged_death_{arm}"], style,
                        linewidth=1.2, label=g["level_label"].iloc[0])
        rp.sweep_panel(left, f"{stratum}\nunweighted", years, legend=False,
                       ylabel=UNFLAGGED_LABEL if c == 0 else None)
        rp.sweep_panel(right, f"{stratum}\nreweighted", years, ynumbers=False)
    rp.save_fig(fig, outdir, f"{name}_{years}y")


def consequences(df):
    sweep = pd.concat([cq.sweep(df, ALL_CAUSE, hz, ARMS)
                       for hz in HORIZONS_YEARS], ignore_index=True)
    outdir = path_results(PHASE)
    for hz in HORIZONS_YEARS:
        draw_pairs(sweep[sweep["horizon_years"] == hz], hz, outdir,
                   "additional_unflagged_death")


def main():
    correction_by_cause()
    df = load_cause(ALL_CAUSE)
    if df is None:
        raise SystemExit("No assembled all-cause predictions. Run phase 2 first.")
    correction_by_stratum(df)
    consequences(df)
    print("phase 4 results written")


if __name__ == "__main__":
    main()
