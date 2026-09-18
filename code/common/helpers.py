"""Constants and shared readers for the reporting stages."""

from pathlib import Path
import numpy as np

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
