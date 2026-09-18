#!/usr/bin/env python3
"""Every reported result of the prediction-bias phase.
"""

from __future__ import annotations

import gc
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "common"))
import consequences as cq
import reporting as rp
from helpers import (ALL_CAUSE, CAUSES, CAUSE_LABELS, compute_rates,
                     load_cause, performance_rows)
from paths import path_results
from survival import HORIZONS_YEARS, risk_column

BIAS_LABEL = "Prediction bias (deaths per 1,000,000)"
UNFLAGGED_LABEL = "Additional unflagged death"
FALSE_NEGATIVE_LABEL = "Difference in false negative rate"

ARMS = [("ukb", "unweighted")]

PHASE = "phase2_prediction_bias"


# ── 1. Prediction bias by cause ──────────────────────────────────────────────

def bias_by_cause():
    rows = []
    for cause in CAUSES:
        df = load_cause(cause)
        if df is None:
            print(f"[SKIP] {cause}")
            continue
        for hz in HORIZONS_YEARS:
            pred = {s: risk_column(cause, s, hz) for s in ("pmr", "ukb")}
            r = compute_rates(df, f"event_{hz}y", pred)
            rows.append({
                "cause": cause, "horizon_years": hz,
                "pmr_bias_per_mil": (r["pmr_pred_rate"] - r["observed_rate"]) * 1e6,
                "ukb_bias_per_mil": (r["ukb_pred_rate"] - r["observed_rate"]) * 1e6,
            })
        del df
        gc.collect()
    out = pd.DataFrame(rows)

    outdir = path_results(PHASE)
    for hz in HORIZONS_YEARS:
        sub = out[out["horizon_years"] == hz].copy()
        sub["cause"] = pd.Categorical(sub["cause"], categories=CAUSES, ordered=True)
        sub = sub.sort_values("cause")
        rows_f = pd.DataFrame({"kind": "level",
                               "label": [CAUSE_LABELS[c] for c in sub["cause"]],
                               "pmr_bias_per_mil": sub["pmr_bias_per_mil"].to_numpy(),
                               "ukb_bias_per_mil": sub["ukb_bias_per_mil"].to_numpy()})
        fig, ax = plt.subplots(figsize=(5.6, 2.8))
        rp.lollipop(ax, rows_f,
                    [("pmr_bias_per_mil", rp.C_PMR, rp.M_PMR),
                     ("ukb_bias_per_mil", rp.C_UKB, rp.M_UKB)], BIAS_LABEL)
        rp.pad_axis(ax, list(rows_f["pmr_bias_per_mil"]) + list(rows_f["ukb_bias_per_mil"]))
        ax.legend(handles=_legend_handles(), fontsize=9, frameon=False, ncol=2,
                  loc="lower center", bbox_to_anchor=(0.5, 1.0))
        fig.tight_layout()
        rp.save_fig(fig, outdir, f"prediction_bias_by_cause_{hz}y")


def _legend_handles():
    from matplotlib.lines import Line2D
    return [Line2D([], [], color=rp.C_PMR, marker=rp.M_PMR, linestyle="none",
                   markersize=6, label="PMR"),
            Line2D([], [], color=rp.C_UKB, marker=rp.M_UKB, linestyle="none",
                   markersize=6, label="UKB")]


# ── 2. Prediction bias by stratum ────────────────────────────────────────────

def bias_by_stratum(df):
    rows = rp.stratum_rows(df)
    outdir = path_results(PHASE)
    for hz in HORIZONS_YEARS:
        pred = {s: risk_column(ALL_CAUSE, s, hz) for s in ("pmr", "ukb")}
        drawn = rows.copy()
        for col in ("pmr_bias_per_mil", "ukb_bias_per_mil"):
            drawn[col] = np.nan
        for i, row in drawn.iterrows():
            if row["kind"] != "level":
                continue
            mask = rp.as_level(df[row["var"]]) == row["level"]
            r = compute_rates(df[mask], f"event_{hz}y", pred)
            drawn.at[i, "pmr_bias_per_mil"] = (r["pmr_pred_rate"] - r["observed_rate"]) * 1e6
            drawn.at[i, "ukb_bias_per_mil"] = (r["ukb_pred_rate"] - r["observed_rate"]) * 1e6

        fig, ax = plt.subplots(figsize=(6.4, max(4.0, len(drawn) * 0.26 + 1.0)))
        rp.lollipop(ax, drawn,
                    [("pmr_bias_per_mil", rp.C_PMR, rp.M_PMR),
                     ("ukb_bias_per_mil", rp.C_UKB, rp.M_UKB)], BIAS_LABEL)
        rp.pad_axis(ax, list(drawn["pmr_bias_per_mil"]) + list(drawn["ukb_bias_per_mil"]))
        ax.legend(handles=_legend_handles(), fontsize=9, frameon=False, ncol=2,
                  loc="lower center", bbox_to_anchor=(0.5, 1.0))
        fig.tight_layout()
        rp.save_fig(fig, outdir, f"prediction_bias_by_strata_{hz}y")



# ── 3 and 4. The two consequence sweeps ──────────────────────────────────────

def draw_sweep(part, years, outdir, value, ylabel, name, select=None):
    strata = list(dict.fromkeys(part["stratum"]))
    ncol = min(3, len(strata))
    nrow = -(-len(strata) // ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 3.0 * nrow),
                             squeeze=False)
    for i, (ax, stratum) in enumerate(zip(axes.ravel(), strata)):
        sub = part[part["stratum"] == stratum]
        levels = list(dict.fromkeys(sub["level"]))
        if select is not None:
            var = sub["variable"].iloc[0]
            want = [str(x) for x in select.get(var, levels)]
            levels = [lv for lv in levels if lv in want]
        for lv in levels:
            g = sub[sub["level"] == lv].sort_values("threshold")
            ax.plot(g["threshold"] * 100, g[value], linewidth=1.2,
                    label=g["level_label"].iloc[0])
        rp.sweep_panel(ax, stratum, years,
                       ylabel if i % ncol == 0 else None)
    for ax in axes.ravel()[len(strata):]:
        ax.axis("off")
    fig.tight_layout()
    rp.save_fig(fig, outdir, f"{name}_{years}y")


def consequences(df):
    frames = [cq.sweep(df, ALL_CAUSE, hz, ARMS) for hz in HORIZONS_YEARS]
    sweep = pd.concat(frames, ignore_index=True)

    outdir = path_results(PHASE)
    for hz in HORIZONS_YEARS:
        part = sweep[sweep["horizon_years"] == hz]
        draw_sweep(part, hz, outdir, "additional_unflagged_death_unweighted",
                   UNFLAGGED_LABEL, "additional_unflagged_death_unweighted")
        draw_sweep(part, hz, outdir, "false_negative_rate_difference_unweighted",
                   FALSE_NEGATIVE_LABEL,
                   "difference_in_false_negative_rate_unweighted",
                   select=rp.FALSE_NEGATIVE_RATE_LEVELS)


# ── 5. The parity workbook ───────────────────────────────────────────────────

CAPACITY = cq.CAPACITY


def _fmt(x, spec="{:.2f}"):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return ""
    return spec.format(x)


def parity(df):
    cals, pars = [], []
    for hz in HORIZONS_YEARS:
        c, p = cq.parity(df, ALL_CAUSE, hz, ARMS)
        cals.append(c)
        pars.append(p)
    cal = pd.concat(cals, ignore_index=True)
    par = pd.concat(pars, ignore_index=True)

    cal = cal[cal["variable"].isin(rp.PARITY_STRATA)]
    par = par[par["variable"].isin(rp.PARITY_STRATA)
              & (par["capacity"].round(4) == round(CAPACITY, 4))]

    cal_r = cq.level_range(cal, "relative_bias_pc", ["variable", "stratum", "horizon_years"])
    par_r = cq.level_range(par, "ppv_pc", ["variable", "stratum", "horizon_years"])
    rows = []
    order = {v: i for i, v in enumerate(rp.PARITY_STRATA)}
    for hz in HORIZONS_YEARS:
        for var in rp.PARITY_STRATA:
            def ratio(frame):
                sel = frame[(frame["variable"] == var) & (frame["horizon_years"] == hz)]
                by = sel.set_index("source")["range"]
                if "pmr" not in by or "ukb" not in by or by["pmr"] == 0:
                    return np.nan
                return float(by["ukb"] / by["pmr"])
            rows.append({
                "Horizon": f"{int(hz)}-year",
                "Stratum": rp.STRATA_LABELS.get(var, var),
                "Calibration parity, relative bias range ratio": _fmt(ratio(cal_r)),
                f"Predictive parity, predictive value range ratio at "
                f"{CAPACITY:.0%} capacity": _fmt(ratio(par_r)),
                "_h": int(hz), "_o": order[var]})
    table = (pd.DataFrame(rows).sort_values(["_h", "_o"], kind="mergesort")
             .drop(columns=["_h", "_o"]).reset_index(drop=True))
    table.loc[table["Horizon"].duplicated(), "Horizon"] = ""
    rp.write_workbook(
        path_results(PHASE, "calibration_predictive_parity.xlsx"),
        {"Parity ratios": table})


# ── 6. The performance workbook ──────────────────────────────────────────────

MODELS = [("PMR-trained", "pmr"), ("UKB-trained", "ukb")]


def performance():
    rows = []
    for cause in CAUSES:
        df = load_cause(cause)
        if df is None:
            print(f"[SKIP] {cause}")
            continue
        rows.extend(performance_rows(df, cause, MODELS))
        del df
        gc.collect()
    rp.write_workbook(path_results(PHASE, "model_performance.xlsx"),
                      {"Model performance": pd.DataFrame(rows)})


def main():
    bias_by_cause()

    df = load_cause(ALL_CAUSE)
    if df is None:
        raise SystemExit("No assembled all-cause predictions. Run stages 1 to 4.")
    bias_by_stratum(df)
    consequences(df)
    parity(df)
    del df
    gc.collect()
    performance()
    print("phase 2 results written")


if __name__ == "__main__":
    main()
