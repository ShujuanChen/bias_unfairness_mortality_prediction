"""The threshold sweep and the two parity measures, for all-cause mortality.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reporting as rp
from survival import risk_column

REFERENCE = "pmr"

_STEP = 0.005
THRESHOLDS = np.round(np.arange(_STEP, 1.0 + _STEP / 2, _STEP), 4)
CAPACITY = 0.05


def capacity_thresholds(pred, w, capacities):
    pred = np.asarray(pred, float)
    w = np.asarray(w, float)
    capacities = np.asarray(capacities, float)
    total = w.sum()
    if total <= 0:
        return np.full(capacities.shape, np.inf)
    order = np.argsort(-pred, kind="stable")
    cum = np.cumsum(w[order]) / total
    idx = np.searchsorted(cum, capacities, side="left")
    idx = np.clip(idx, 0, len(order) - 1)
    return pred[order[idx]]


def _prepare(df, cause, years, arms):
    sources = [REFERENCE] + [s for s, _ in arms]
    pred = {}
    for s in sources:
        col = risk_column(cause, s, years)
        if col not in df.columns:
            raise SystemExit(f"{cause} {years}y: no {col} column.")
        pred[s] = df[col].to_numpy(float)

    event = df[f"event_{years}y"].to_numpy(float)
    w = df["w"].to_numpy(float)
    ok = np.isfinite(event) & np.isfinite(w) & (w > 0)
    for p in pred.values():
        ok &= np.isfinite(p)

    spec = rp.level_codes(df)
    for s in spec:
        s["codes"] = s["codes"][ok]
    return spec, event[ok], w[ok], {s: p[ok] for s, p in pred.items()}


def sweep(df, cause, years, arms):
    spec, event, w, pred = _prepare(df, cause, years, arms)
    wy = w * event
    sources = [REFERENCE] + [s for s, _ in arms]

    base = {s["var"]: (rp.accumulate(s["codes"], len(s["levels"]), w),
                       rp.accumulate(s["codes"], len(s["levels"]), wy))
            for s in spec}

    rows = []
    for tau in THRESHOLDS:
        flagged = {s: pred[s] >= float(tau) for s in sources}
        unflagged = {}
        for s in spec:
            k = len(s["levels"])
            for src in sources:
                unflagged[(src, s["var"])] = base[s["var"]][1] - rp.accumulate(
                    s["codes"], k, wy, flagged[src])
        for s in spec:
            wg, dg = base[s["var"]]
            ur = unflagged[(REFERENCE, s["var"])]
            for j, lv in enumerate(s["levels"]):
                row = {"cause": cause, "horizon_years": years,
                       "variable": s["var"], "stratum": s["title"],
                       "level": str(lv), "level_index": j,
                       "level_label": rp.level_label(s["var"], lv),
                       "threshold": float(tau),
                       "population_count": float(wg[j]),
                       "deaths_count": float(dg[j]),
                       "deaths_per_mil": float(dg[j] / wg[j] * 1e6)}
                for src, arm in arms:
                    diff = unflagged[(src, s["var"])][j] - ur[j]
                    row[f"additional_unflagged_death_{arm}"] = float(diff / wg[j] * 1e6)
                    row[f"false_negative_rate_difference_{arm}"] = float(diff / dg[j] * 100)
                rows.append(row)
    return pd.DataFrame(rows)


def parity(df, cause, years, arms):
    spec, event, w, pred = _prepare(df, cause, years, arms)
    wy = w * event
    sources = [REFERENCE] + [s for s, _ in arms]
    base = {s["var"]: (rp.accumulate(s["codes"], len(s["levels"]), w),
                       rp.accumulate(s["codes"], len(s["levels"]), wy))
            for s in spec}

    cal = []
    for src in sources:
        for s in spec:
            k = len(s["levels"])
            wg, dg = base[s["var"]]
            wp = rp.accumulate(s["codes"], k, w * pred[src])
            for j, lv in enumerate(s["levels"]):
                obs, prd = dg[j] / wg[j], wp[j] / wg[j]
                sel = s["codes"] == j
                cal.append({"cause": cause, "horizon_years": years,
                            "variable": s["var"], "stratum": s["title"],
                            "level": str(lv), "level_index": j,
                            "level_label": rp.level_label(s["var"], lv),
                            "source": src, "n": int(sel.sum()),
                            "population_count": float(wg[j]),
                            "deaths_count": float(dg[j]),
                            "observed_risk_per_mil": float(obs * 1e6),
                            "predicted_risk_per_mil": float(prd * 1e6),
                            "relative_bias_pc": float((prd - obs) / obs * 100)})

    par = []
    caps = {s: float(capacity_thresholds(pred[s], w, [CAPACITY])[0])
            for s in sources}
    for src in sources:
        flag = pred[src] >= caps[src]
        for s in spec:
            k = len(s["levels"])
            wg, dg = base[s["var"]]
            nf = rp.accumulate(s["codes"], k, w, flag)
            nfd = rp.accumulate(s["codes"], k, wy, flag)
            for j, lv in enumerate(s["levels"]):
                par.append({"cause": cause, "horizon_years": years,
                            "variable": s["var"], "stratum": s["title"],
                            "level": str(lv), "level_index": j,
                            "level_label": rp.level_label(s["var"], lv),
                            "source": src, "capacity": CAPACITY,
                            "threshold": caps[src],
                            "flagged_per_mil": float(nf[j] / wg[j] * 1e6),
                            "ppv_pc": float(nfd[j] / nf[j] * 100) if nf[j] > 0 else np.nan,
                            "sensitivity_pc": float(nfd[j] / dg[j] * 100)})

    return pd.DataFrame(cal), pd.DataFrame(par)


def level_range(frame, value, by):
    g = frame.groupby(by + ["source"])[value].agg(["min", "max"]).reset_index()
    g["range"] = g["max"] - g["min"]
    return g
