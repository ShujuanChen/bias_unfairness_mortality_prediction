#!/usr/bin/env python3
"""Collect one prediction per individual per training source.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "common"))
from tabular import write_table

import survival as dl
from cohort import prepare_data
from paths import path_temp


def _load_block(d: Path) -> pd.DataFrame:
    path = d / "pmr_predictions.npz"
    if not path.exists():
        raise RuntimeError(f"Missing model outputs at {d}. Train this (cause, source) first.")
    a = np.load(path)
    out = {"pmr_id": a["pmr_id"].astype(int)}
    out.update({f"risk_{y}y": a[f"risk_{y}y"].astype(float)
                for y in dl.HORIZONS_YEARS})
    return pd.DataFrame(out)


def collect(cause: str, source: str, expected: set) -> pd.DataFrame:
    """One prediction per individual for a single source."""
    if source == "pmr":
        preds = pd.concat([_load_block(dl.get_model_dir(cause, source, f)) for f in dl.FOLDS],
                          ignore_index=True)
        if preds["pmr_id"].duplicated().any():
            raise RuntimeError(
                f"cause={cause} source={source}: fold blocks overlap, so an "
                "individual was scored by more than one model. The evaluation "
                "would be partly in sample."
            )
    else:
        preds = _load_block(dl.get_model_dir(cause, source, None))
        if preds["pmr_id"].duplicated().any():
            raise RuntimeError(f"cause={cause} source={source}: duplicate identifiers.")

    got = set(preds["pmr_id"].astype(int))
    if got != expected:
        raise RuntimeError(
            f"cause={cause} source={source}: predictions cover {len(got)} "
            f"individuals against an evaluation cohort of {len(expected)}; "
            f"missing {len(expected - got)}, unexpected {len(got - expected)}."
        )
    return preds.set_index("pmr_id")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cause", required=True, choices=dl.OUTCOMES)
    ap.add_argument("--sources", default="ukb,pmr",
                    help="comma-separated training sources to assemble")
    args = ap.parse_args()

    sources = [s.strip() for s in args.sources.split(",") if s.strip()]
    unknown = [s for s in sources if s not in dl.TRAINING_SOURCES]
    if unknown:
        raise RuntimeError(f"Unknown training source(s): {unknown}")

    pmr = prepare_data(args.cause)["pmr"]
    keep = (["pmr_id", "w", "dod_deaths", "RGN11CD"]
            + [f"event_{y}y" for y in dl.HORIZONS_YEARS])
    out = pmr[[c for c in keep + dl.RHS_COVARIATES + dl.TARGET_CARRIED_COLUMNS
               if c in pmr.columns]].copy()
    out["outcome"] = args.cause
    expected = set(out["pmr_id"].astype(int))

    for source in sources:
        preds = collect(args.cause, source, expected)
        indexed = out.set_index("pmr_id")
        for y in dl.HORIZONS_YEARS:
            indexed[f"{source}_pred_{y}y"] = preds[f"risk_{y}y"]
        out = indexed.reset_index()
        dl.log_step(f"{source}: {len(preds)} predictions assembled")

    path = path_temp("predictions", f"{args.cause}_individual_predictions.csv")
    write_table(out, path)
    dl.log_step(f"wrote {path} (rows={len(out)}, sources={','.join(sources)})")


if __name__ == "__main__":
    main()
