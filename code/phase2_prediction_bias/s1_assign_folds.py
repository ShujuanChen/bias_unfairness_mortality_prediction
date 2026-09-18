#!/usr/bin/env python3
"""Build the PMR fold assignment for one cause.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "common"))

import folds
import survival as sv
from cohort import prepare_data
from paths import require_pipeline


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cause", required=True, choices=sv.OUTCOMES)
    args = ap.parse_args()

    require_pipeline("outcome", "The outcome model")

    pmr = prepare_data(args.cause)["pmr"]
    if "pmr_id" not in pmr.columns:
        raise RuntimeError(
            "PMR is missing the stable 'pmr_id' column. It is assigned in phase 1 "
            "stage 2, before any complete-case filter, so that the split does not "
            "depend on which variables a given pipeline needs."
        )
    folds.build(args.cause,
                ids=pmr["pmr_id"].astype(int),
                status=pmr["status"].to_numpy(dtype=int),
                frame=pmr)


if __name__ == "__main__":
    main()
