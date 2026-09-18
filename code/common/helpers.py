"""Constants and shared readers for the reporting stages."""

from pathlib import Path
import numpy as np
import pandas as pd

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parent))
from tabular import read_table
from paths import path_temp
import survival as sv

CAUSE_LABELS = {
    "all_cause_mortality": "All-cause",
    "cancer_mortality": "Cancer",
    "cardiovascular_mortality": "Cardiovascular",
    "respiratory_mortality": "Respiratory",
    "digestive_mortality": "Digestive",
}
CAUSES = list(CAUSE_LABELS)
_reportable = {k for k, v in sv.setting("phase2", "outcomes").items()
               if v["type"] != "residual"}
if set(CAUSES) != _reportable:
    raise RuntimeError(
        "This module reports " + ", ".join(sorted(CAUSES))
        + " and framework_config.json models " + ", ".join(sorted(_reportable))
        + " as non-residual causes. They must be the same set.")

ALL_CAUSE = "all_cause_mortality"

WEIGHT_METHOD = str(sv.setting("phase3", "weight_method"))
WEIGHT_METHOD_OPTIONS = list(sv.setting("phase3", "weight_method_options"))
if WEIGHT_METHOD not in WEIGHT_METHOD_OPTIONS:
    raise RuntimeError(
        f"phase3.weight_method is {WEIGHT_METHOD!r}, which is not one of "
        + ", ".join(WEIGHT_METHOD_OPTIONS) + ".")
if WEIGHT_METHOD not in sv.setting("phase2", "weight_sources"):
    raise RuntimeError(
        f"phase3.weight_method is {WEIGHT_METHOD!r} and phase2.weight_sources "
        "does not declare it.")
WUKB_SOURCE = "ukbw_" + WEIGHT_METHOD
WEIGHT_METHOD_LABEL = sv.setting("phase2", "weight_sources")[WEIGHT_METHOD]["label"]


def find_individual_predictions(cause):
    return path_temp("predictions", f"{cause}_individual_predictions.csv", create=False)


def load_cause(cause):
    path = find_individual_predictions(cause)
    if not path.exists():
        return None
    return read_table(path)


def require_weighted_column(df, cause, column):
    if column not in df.columns:
        raise SystemExit(
            f"{cause}: no {column} column. The weighted series is "
            f"{WUKB_SOURCE}, which framework_config.json names under "
            "phase3.weight_method. Run phase 3 for that weight set.")


def weighted_mean(vals, w):
    ok = np.isfinite(vals) & np.isfinite(w) & (w > 0)
    if not np.any(ok):
        return np.nan
    return float(np.sum(vals[ok] * w[ok]) / np.sum(w[ok]))


def compute_rates(df, event_col, pred_cols, w_col="w"):
    obs = df[event_col].to_numpy(dtype=float)
    w = df[w_col].to_numpy(dtype=float) if w_col in df.columns else np.ones(len(df))
    obs_rate = weighted_mean(obs, w)
    result = {"observed_rate": obs_rate, "n": len(df), "n_events": int(np.nansum(obs))}
    for name, col in pred_cols.items():
        if col in df.columns:
            pred = df[col].to_numpy(dtype=float)
            result[f"{name}_pred_rate"] = weighted_mean(pred, w)
            result[f"{name}_relative_bias"] = (result[f"{name}_pred_rate"] - obs_rate) / obs_rate if obs_rate > 0 else np.nan
        else:
            result[f"{name}_pred_rate"] = np.nan
            result[f"{name}_relative_bias"] = np.nan
    return result


def follow_up(death_date):
    d = pd.to_datetime(death_date, errors="coerce").values.astype("datetime64[D]")
    has_date = ~pd.isna(d)
    seen = has_date & (d <= sv.T_ADMIN_END)
    end = np.full(len(d), sv.T_ADMIN_END, dtype="datetime64[D]")
    end[has_date] = np.minimum(d[has_date], sv.T_ADMIN_END)
    days = (end - sv.T0).astype("timedelta64[D]").astype(float)
    return days, seen.astype(float)


def weighted_auc(pred, event, w):
    order = np.argsort(pred, kind="mergesort")
    pred, event, w = pred[order], event[order], w[order]
    group = np.cumsum(np.r_[True, pred[1:] != pred[:-1]]) - 1
    n_groups = int(group[-1]) + 1
    w1 = np.bincount(group, weights=w * event, minlength=n_groups)
    w0 = np.bincount(group, weights=w * (1 - event), minlength=n_groups)
    below = np.r_[0.0, np.cumsum(w0)[:-1]]
    if w1.sum() <= 0 or w0.sum() <= 0:
        return np.nan
    return float((np.sum(w1 * below) + 0.5 * np.sum(w1 * w0))
                 / (w1.sum() * w0.sum()))


def brier(pred, event, w):
    total = w.sum()
    if total <= 0:
        return np.nan
    return float(np.sum(w * (pred - event) ** 2) / total)


class _Fenwick:

    def __init__(self, n):
        self.n = n
        self.t = np.zeros(n + 1, dtype=float)

    def add(self, i, v):
        i += 1
        while i <= self.n:
            self.t[i] += v
            i += i & (-i)

    def prefix(self, i):
        s = 0.0
        while i > 0:
            s += self.t[i]
            i -= i & (-i)
        return s


def _death_pair_mass(risk, time, w):
    if len(risk) == 0:
        return 0.0, 0.0, 0.0
    rank = pd.Series(risk).rank(method="dense").to_numpy().astype(int) - 1
    n_ranks = int(rank.max()) + 1
    order = np.argsort(-time, kind="mergesort")
    time, w, rank = time[order], w[order], rank[order]
    tree = _Fenwick(n_ranks)
    equal = np.zeros(n_ranks, dtype=float)
    concordant = tied = total = 0.0
    i, n = 0, len(time)
    while i < n:
        j = i
        while j < n and time[j] == time[i]:
            j += 1
        inserted = tree.prefix(n_ranks)
        for k in range(i, j):
            below = tree.prefix(rank[k])
            concordant += w[k] * below
            tied += w[k] * equal[rank[k]]
            total += w[k] * inserted
        for k in range(i, j):
            tree.add(rank[k], w[k])
            equal[rank[k]] += w[k]
        i = j
    return concordant, tied, total


def _cross_pair_mass(risk_a, w_a, risk_b, w_b):
    if len(risk_a) == 0 or len(risk_b) == 0:
        return 0.0, 0.0, 0.0
    order = np.argsort(risk_b, kind="mergesort")
    rb, wb = risk_b[order], w_b[order]
    cum = np.r_[0.0, np.cumsum(wb)]
    lo = np.searchsorted(rb, risk_a, side="left")
    hi = np.searchsorted(rb, risk_a, side="right")
    return (float(np.sum(w_a * cum[lo])),
            float(np.sum(w_a * (cum[hi] - cum[lo]))),
            float(w_a.sum() * w_b.sum()))


def harrell_c(risk, time, event, w):
    ok = (np.isfinite(risk) & np.isfinite(time) & np.isfinite(event)
          & np.isfinite(w) & (w > 0))
    risk, time, event, w = risk[ok], time[ok], event[ok], w[ok]
    if len(risk) == 0 or event.sum() == 0:
        return np.nan
    alive = event == 0
    t_end = time.max()
    if alive.any() and not np.all(time[alive] == t_end):
        raise RuntimeError(
            "The concordance index here counts every death against every "
            "survivor, which holds because censoring is administrative and "
            "everyone alive leaves on the last day of follow-up. This "
            "population has survivors leaving earlier.")
    died = ~alive
    early = died & (time < t_end)
    concordant, tied, total = _cross_pair_mass(risk[early], w[early],
                                               risk[alive], w[alive])
    c, t, n = _death_pair_mass(risk[died], time[died], w[died])
    concordant, tied, total = concordant + c, tied + t, total + n
    if total <= 0:
        return np.nan
    return float((concordant + 0.5 * tied) / total)


def performance_rows(df, cause, models):
    rows = []
    time = seen = None
    if cause not in sv.COMPETING_CAUSES and "dod_deaths" in df.columns:
        time, seen = follow_up(df["dod_deaths"])
    w_all = (df["w"].to_numpy(dtype=float) if "w" in df.columns
             else np.ones(len(df)))
    for label, source in models:
        for hz in sv.HORIZONS_YEARS:
            col = sv.risk_column(cause, source, hz)
            if col not in df.columns:
                continue
            pred = df[col].to_numpy(dtype=float)
            event = df[f"event_{hz}y"].to_numpy(dtype=float)
            ok = (np.isfinite(pred) & np.isfinite(event) & np.isfinite(w_all)
                  & (w_all > 0))
            c_index = (np.nan if time is None
                       else harrell_c(pred[ok], time[ok], seen[ok], w_all[ok]))
            p, e, w = pred[ok], event[ok], w_all[ok]
            rows.append({
                "Cause": CAUSE_LABELS[cause],
                "Model": label,
                "Horizon": f"{int(hz)}-year",
                "Records": int(len(p)),
                "Deaths": int(e.sum()),
                "Weighted AUC": round(weighted_auc(p, e, w), 4),
                "Weighted c-index": (np.nan if not np.isfinite(c_index)
                                     else round(c_index, 4)),
                "Weighted Brier score": round(brier(p, e, w), 6),
            })
    return rows
