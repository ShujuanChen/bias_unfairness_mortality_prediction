#!/usr/bin/env python3
"""Convert cause-specific hazards into cumulative incidence.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "common"))
from tabular import read_table, write_table

import survival as sv
from paths import path_temp

COMPETING = sv.COMPETING_CAUSES
HORIZONS = sv.HORIZONS_DAYS

CHUNK = 200_000


def _model_dirs(cause: str, source: str):
    if source == "pmr":
        return [sv.get_model_dir(cause, source, f) for f in sv.FOLDS]
    return [sv.get_model_dir(cause, source, None)]


def step_interp(times: np.ndarray, values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    if len(times) == 0:
        return np.zeros_like(grid, dtype=float)
    idx = np.searchsorted(times, grid, side="right") - 1
    return np.where(idx >= 0, values[np.clip(idx, 0, len(values) - 1)], 0.0).astype(float)


def load_cause(cause: str, source: str, grid: np.ndarray):
    ids, lp, model_of_row, baselines = [], [], [], []
    for k, d in enumerate(_model_dirs(cause, source)):
        preds = np.load(d / "pmr_predictions.npz")
        if "lp" not in preds:
            raise RuntimeError(
                f"{d / 'pmr_predictions.npz'} has no linear predictor. Retrain "
                "this (cause, source)."
            )
        base = np.load(d / "baseline_hazard.npz")
        ids.append(preds["pmr_id"].astype(int))
        lp.append(preds["lp"].astype(float))
        model_of_row.append(np.full(len(preds["pmr_id"]), k, dtype=int))
        baselines.append(step_interp(base["event_times"].astype(float),
                                     base["cum_baseline_hazard"].astype(float), grid))
    out = pd.DataFrame({"pmr_id": np.concatenate(ids),
                        "lp": np.concatenate(lp),
                        "model": np.concatenate(model_of_row)})
    if out["pmr_id"].duplicated().any():
        raise RuntimeError(f"cause={cause} source={source}: an individual is scored twice.")
    return out.sort_values("pmr_id").reset_index(drop=True), np.vstack(baselines)


def cumulative_incidence(source: str) -> pd.DataFrame:
    """Cumulative incidence at each horizon, for every cause, for one source."""
    grid = np.arange(0.0, float(max(HORIZONS)) + 1.0)
    horizon_at = {hz: int(np.searchsorted(grid, float(hz))) for hz in HORIZONS}

    per_cause = {}
    for c in COMPETING:
        per_cause[c] = load_cause(c, source, grid)

    ids = per_cause[COMPETING[0]][0]["pmr_id"].to_numpy()
    for c in COMPETING[1:]:
        if not np.array_equal(ids, per_cause[c][0]["pmr_id"].to_numpy()):
            raise RuntimeError(
                f"source={source}: cause {c} covers a different set of individuals "
                "than the others. Cumulative incidence needs all causes on the "
                "same people."
            )

    n = len(ids)
    out = {"pmr_id": ids}
    for hz in HORIZONS:
        yrs = int(round(hz / 365.25))
        for c in COMPETING:
            out[f"{c}_cif_{yrs}y_{source}"] = np.empty(n)
        out[f"survival_{yrs}y_{source}"] = np.empty(n)

    for start in range(0, n, CHUNK):
        stop = min(start + CHUNK, n)
        m = stop - start
        risk = {c: np.exp(per_cause[c][0]["lp"].to_numpy()[start:stop]) for c in COMPETING}
        model_of = {c: per_cause[c][0]["model"].to_numpy()[start:stop] for c in COMPETING}

        cif = {c: np.zeros(m) for c in COMPETING}
        H_prev = {c: np.zeros(m) for c in COMPETING}
        surv_prev = np.ones(m)
        for j in range(len(grid)):
            dH = {}
            for c in COMPETING:
                H = risk[c] * per_cause[c][1][model_of[c], j]
                dH[c] = H - H_prev[c]
                H_prev[c] = H
            surv = np.exp(-sum(H_prev.values()))

            died = surv_prev - surv
            dH_total = sum(dH.values())
            share = np.divide(died, dH_total, out=np.zeros(m), where=dH_total > 0)
            for c in COMPETING:
                cif[c] += dH[c] * share
            for hz, jj in horizon_at.items():
                if j == jj:
                    yrs = int(round(hz / 365.25))
                    for c in COMPETING:
                        out[f"{c}_cif_{yrs}y_{source}"][start:stop] = cif[c]
                    out[f"survival_{yrs}y_{source}"][start:stop] = surv
            surv_prev = surv
        sv.log_step(f"  {source}: {stop:,} of {n:,}")

    for hz in HORIZONS:
        yrs = int(round(hz / 365.25))
        closure = out[f"survival_{yrs}y_{source}"].copy()
        for c in COMPETING:
            closure = closure + out[f"{c}_cif_{yrs}y_{source}"]

        if not np.allclose(closure, 1.0, atol=1e-8):
            raise RuntimeError(
                f"source={source} horizon={yrs}y: cumulative incidences plus "
                f"survival range from {closure.min():.9f} to {closure.max():.9f} "
                "rather than summing to one. The cause set is not exhaustive, or "
                "a cause is missing individuals."
            )
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sources", default="ukb,pmr",
                    help="comma-separated training sources to convert")
    args = ap.parse_args()
    sources = [s.strip() for s in args.sources.split(",") if s.strip()]

    cif = None
    for source in sources:
        sv.log_step(f"cumulative incidence: {source}")
        df = cumulative_incidence(source)
        cif = df if cif is None else cif.merge(df, on="pmr_id", how="inner")

    path = path_temp("predictions", "cumulative_incidence.csv")
    write_table(cif, path)
    sv.log_step(f"wrote {path} (rows={len(cif)})")

    for cause in COMPETING:
        target = path_temp("predictions", f"{cause}_individual_predictions.csv", create=False)
        if not target.exists():
            raise RuntimeError(f"Missing {target}. Run stage 3 for every cause first.")
        preds = read_table(target)
        cols = {f"{cause}_cif_{int(round(hz / 365.25))}y_{s}":
                f"{s}_cif_{int(round(hz / 365.25))}y"
                for hz in HORIZONS for s in sources}

        preds = preds.drop(columns=[c for c in cols.values() if c in preds.columns])
        merged = preds.merge(cif[["pmr_id"] + list(cols)].rename(columns=cols),
                             on="pmr_id", how="left")
        if merged[list(cols.values())].isna().any().any():
            raise RuntimeError(
                f"{cause}: some individuals have no cumulative incidence. The "
                "prediction table and the fitted models cover different people."
            )
        write_table(merged, target)
        sv.log_step(f"added cumulative incidence to {target}")


if __name__ == "__main__":
    main()
