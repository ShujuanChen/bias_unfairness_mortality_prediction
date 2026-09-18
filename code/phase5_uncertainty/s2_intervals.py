#!/usr/bin/env python3
"""Intervals on the bias and on the correction, from the finished replicates.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from tabular import read_table, write_table
from paths import path_results, path_temp, setting              # noqa: E402
from reporting import as_level                                  # noqa: E402

PHASE = "phase5_uncertainty"

POINT_COLUMNS = ("observed_per_mil", "ukb_per_mil", "pmr_per_mil")

POINT_ESTIMATE_SOURCE = {
    "cause": "bias_correction_by_cause.csv",
    "subgroup": "bias_correction_by_subgroup.csv",
}
POINT_ESTIMATE_COLUMNS = {
    "bias": "ukb_bias_per_mil",
    "observed_per_mil": "observed_per_mil",
    "ukb_per_mil": "ukb_predicted_per_mil",
    "pmr_per_mil": "pmr_predicted_per_mil",
}
WEIGHTED_POINT_COLUMNS = {
    "rate": "wukb_predicted_per_mil",
    "correction": "bias_correction_pct",
}


def attach_point_estimates(frame, shape, keys):
    path = Path(path_temp(POINT_ESTIMATE_SOURCE[shape], create=False))
    if not path.exists():
        raise RuntimeError(
            f"No {path.name} in this run's results. Phase 4 writes the point "
            "estimates that phase 5 reports its intervals against, so run "
            "phase 4 before this stage.")
    point = read_table(path)

    wanted = dict(POINT_ESTIMATE_COLUMNS)
    if "weight_source" in point.columns and point["weight_source"].nunique() == 1:
        source = str(point["weight_source"].iloc[0])
        for stem, col in WEIGHTED_POINT_COLUMNS.items():
            if f"{source}_{stem}_mean" in frame.columns and col in point.columns:
                wanted[f"{source}_{stem}"] = col

    missing = [c for c in wanted.values() if c not in point.columns]
    if missing:
        raise RuntimeError(
            f"{path.name} is missing the columns phase 5 reports against: "
            + ", ".join(missing))

    point = point[keys + list(wanted.values())].copy()
    point.columns = keys + [f"{stem}_point" for stem in wanted]

    frame = frame.copy()
    matched = [f"__key{i}" for i in range(len(keys))]
    for key, col in zip(keys, matched):
        frame[col] = as_level(frame[key])
        point[col] = as_level(point[key])

    left = set(map(tuple, frame[matched].to_numpy()))
    right = set(map(tuple, point[matched].to_numpy()))
    only_interval = sorted(left - right)
    if only_interval:
        raise RuntimeError(
            f"{len(only_interval)} row(s) have an interval and no point estimate "
            f"in {path.name}, the first being {only_interval[0]}. An interval "
            "belongs to a reported quantity, so this means the two did not come "
            "from the same analysis.")

    only_point = sorted(right - left)
    if only_point:
        raise RuntimeError(
            f"{len(only_point)} stratum row(s) are in {path.name} and have no "
            f"interval, the first being {only_point[0]}. The two tables are "
            "built from the same strata, so a mismatch means they did not come "
            "from the same analysis.")

    merged = frame.merge(point.drop(columns=keys), on=matched, how="left",
                         validate="one_to_one").drop(columns=matched)

    ordered = list(keys) + [c for c in ("cause_label", "n_replicates")
                            if c in merged.columns]
    for stem in wanted:
        for suffix in ("point", "mean", "lo", "hi"):
            col = f"{stem}_{suffix}"
            if col in merged.columns and col not in ordered:
                ordered.append(col)
    ordered += [c for c in merged.columns if c not in ordered]
    return merged[ordered]


def _weighted_columns(df: pd.DataFrame) -> list:
    return [c for c in df.columns
            if c.startswith("ukbw_") and c.endswith("_per_mil")]


def load_replicates(root: Path) -> pd.DataFrame:
    name = "rates.csv"
    frames = []
    for d in sorted(root.glob("replicate_*")):
        path = d / name
        if not path.exists():
            continue
        f = read_table(path)
        f["replicate"] = int(d.name.split("_")[1])
        frames.append(f)
    if not frames:
        raise RuntimeError(
            f"No replicate {name} tables under {root}. Run the replicates first.")
    return pd.concat(frames, ignore_index=True)


def summarise(reps: pd.DataFrame, level: float) -> pd.DataFrame:
    lo_q, hi_q = (1 - level) / 2, 1 - (1 - level) / 2
    weighted = _weighted_columns(reps)
    keys = ["cause", "horizon_years", "variable", "level"]

    records = []
    for key, g in reps.groupby(keys, sort=False):
        row = dict(zip(keys, key if isinstance(key, tuple) else (key,)))
        row["n_replicates"] = len(g)

        observed = g["observed_per_mil"].to_numpy(dtype=float)
        unweighted = g["ukb_per_mil"].to_numpy(dtype=float)

        # The uncorrected bias, as a distribution over replicates.
        bias = unweighted - observed
        row["bias_mean"] = float(bias.mean())
        row["bias_lo"] = float(np.quantile(bias, lo_q))
        row["bias_hi"] = float(np.quantile(bias, hi_q))

        for col in POINT_COLUMNS:
            v = g[col].to_numpy(dtype=float)
            row[f"{col}_mean"] = float(v.mean())
            row[f"{col}_lo"] = float(np.quantile(v, lo_q))
            row[f"{col}_hi"] = float(np.quantile(v, hi_q))

        for col in weighted:
            v = g[col].to_numpy(dtype=float)
            key_name = col[:-len("_per_mil")]
            row[f"{key_name}_rate_mean"] = float(v.mean())
            row[f"{key_name}_rate_lo"] = float(np.quantile(v, lo_q))
            row[f"{key_name}_rate_hi"] = float(np.quantile(v, hi_q))

            denominator = np.abs(unweighted - observed)
            correction = np.where(denominator != 0,
                                  100.0 * (denominator - np.abs(v - observed))
                                  / denominator,
                                  np.nan)
            good = np.isfinite(correction)
            if good.sum() < len(correction):
                row[f"{key_name}_correction_dropped"] = int((~good).sum())
            correction = correction[good]
            if len(correction):
                row[f"{key_name}_correction_mean"] = float(correction.mean())
                row[f"{key_name}_correction_lo"] = float(np.quantile(correction, lo_q))
                row[f"{key_name}_correction_hi"] = float(np.quantile(correction, hi_q))
                row[f"{key_name}_correction_median"] = float(np.median(correction))
        records.append(row)

    return pd.DataFrame(records)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--level", type=float,
                    default=float(setting("phase5", "interval_level")))
    ap.add_argument("--bootstrap-dir", default="bootstrap_multiplier",
                    help="replicate directory under temp/")
    args = ap.parse_args()

    root = Path(path_temp(args.bootstrap_dir, create=False))
    reps = load_replicates(root)
    n = reps["replicate"].nunique()

    needed = int(np.ceil(2 / ((1 - args.level) / 2)))
    if n < needed:
        print(f"WARNING: {n} replicate(s) support a {args.level:.0%} percentile "
              f"interval poorly. Each tail holds {n * (1 - args.level) / 2:.1f} "
              f"replicates against the {needed} this level wants, so the bounds "
              f"are the extreme replicates rather than estimated quantiles.",
              file=sys.stderr)

    if "variable" not in reps.columns:
        reps["variable"], reps["level"] = "overall", "all"
    summary = summarise(reps, args.level)

    overall = summary[summary["variable"] == "overall"].drop(
        columns=["variable", "level"])
    subgroup = summary[summary["variable"] != "overall"]

    for frame, shape, keys, name in (
            (overall, "cause", ["cause", "horizon_years"],
             "bootstrap_intervals_by_cause"),
            (subgroup, "subgroup", ["cause", "horizon_years", "variable", "level"],
             "bootstrap_intervals_by_subgroup")):
        if len(frame):
            frame = attach_point_estimates(frame, shape, keys)
        frame.insert(len(keys), "interval_level", args.level)
        out = Path(path_results(PHASE)) / f"{name}.xlsx"
        out.parent.mkdir(parents=True, exist_ok=True)
        write_table(frame, out)
        print(f"{n} replicates, {len(frame)} rows, wrote {out}")

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "phase4_correction"))
    import s1_results as correction

    df = correction.load_cause(correction.ALL_CAUSE)
    correction.correction_by_cause()
    if df is not None:
        correction.correction_by_stratum(df)
    else:
        print("no assembled all-cause predictions, stratum figure not redrawn")


if __name__ == "__main__":
    main()
