#!/usr/bin/env python3
"""Observed and predicted rates for one replicate and a chosen set of sources.
"""

from __future__ import annotations

import argparse
import gc
import sys
from pathlib import Path


import pandas as pd

COMMON = Path(__file__).resolve().parents[1] / "common"
sys.path.insert(0, str(COMMON))
from tabular import write_table
from helpers import ALL_CAUSE, CAUSES, compute_rates, load_cause
from paths import setting
from reporting import STRATA, as_level
from survival import risk_column

HORIZONS = list(setting("phase1", "horizons_years"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sources", required=True,
                    help="comma-separated, for example ukb,ukbw_hse_superlearner")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    sources = [s.strip() for s in a.sources.split(",") if s.strip()]

    rows = []
    for cause in CAUSES:
        df = load_cause(cause)
        if df is None:
            print(f"no assembled predictions for {cause}, skipped")
            continue
        for hz in HORIZONS:
            cols = {}
            for s in sources:
                c = risk_column(cause, s, hz)
                if c in df.columns:
                    cols[s] = c
                else:
                    print(f"{cause} {hz}y: {c} not assembled, skipped")
            def emit(frame, variable, level):
                r = compute_rates(frame, f"event_{hz}y", cols)
                row = {"cause": cause, "horizon_years": hz,
                       "variable": variable, "level": level,
                       "n_evaluated": int(r["n"]), "n_events": int(r["n_events"]),
                       "observed_per_mil": r["observed_rate"] * 1e6}
                for s in cols:
                    row[f"{s}_per_mil"] = r[f"{s}_pred_rate"] * 1e6
                rows.append(row)

            emit(df, "overall", "all")
            if cause != ALL_CAUSE:
                continue
            for var, levels in STRATA:
                if var not in df.columns:
                    continue
                vals = as_level(df[var])
                for level in levels:
                    mask = vals == level
                    if not mask.any():
                        continue
                    emit(df[mask], var, level)
        del df
        gc.collect()

    out = pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    write_table(out, a.out)
    pd.set_option("display.width", 200)
    print(out.to_string(index=False))
    print(f"\nWrote {a.out}")


if __name__ == "__main__":
    main()
