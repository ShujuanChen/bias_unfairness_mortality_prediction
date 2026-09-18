#!/usr/bin/env python3
"""Draw one replicate of the multiplier bootstrap.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from tabular import read_table, write_table
from paths import analytic_ids, path_temp, pmr_analytic_ids, setting

REFERENCE_SAMPLE = ("census" if str(setting("phase3", "weight_method"))
                    .startswith("census_") else "hse")

SEED_BASE = int(setting("phase5", "seed"))

DRAW_SUBDIR = "bootstrap_multiplier"


def _replicate_seed(replicate: int) -> int:
    return SEED_BASE + int(replicate)


def _normalised_exponential(n: int, rng) -> np.ndarray:
    m = rng.exponential(scale=1.0, size=int(n))
    total = m.sum()
    if not np.isfinite(total) or total <= 0:
        raise SystemExit("The multiplier draw produced a non-positive total.")
    return m * (len(m) / total)


def _cluster_multiplier(strata: np.ndarray, clusters: np.ndarray, rng) -> np.ndarray:
    m = np.zeros(len(strata), dtype=float)
    for value in pd.unique(strata):
        rows = np.flatnonzero(strata == value)
        units = pd.unique(clusters[rows])
        per_unit = dict(zip(units, _normalised_exponential(len(units), rng)))
        m[rows] = [per_unit[c] for c in clusters[rows]]
    return m


def _stratified_multiplier(strata: np.ndarray, rng, take_in_full=()) -> np.ndarray:
    m = np.zeros(len(strata), dtype=float)
    for value in pd.unique(strata):
        idx = np.flatnonzero(strata == value)
        if value in take_in_full:
            m[idx] = 1.0
            continue
        m[idx] = _normalised_exponential(len(idx), rng)
    return m


def _summary(m: np.ndarray) -> dict:
    m = np.asarray(m, dtype=float)
    ss = float((m ** 2).sum())
    return {"mean": float(m.mean()), "variance": float(m.var(ddof=0)),
            "min": float(m.min()), "max": float(m.max()),
            "effective_sample_size": float(m.sum() ** 2 / ss) if ss > 0 else 0.0,
            "ess_over_n": float((m.sum() ** 2 / ss) / len(m)) if ss > 0 else 0.0}


def _harmonised(name: str) -> Path:
    return path_temp("harmonised", name, create=False)


def draw(replicate: int, out_root: Path) -> dict:
    rng = np.random.default_rng(_replicate_seed(replicate))

    # ---- cohort ---------------------------------------------------------------
    ukb = read_table(_harmonised("ukb_with_pmr.csv"), usecols=["participant_id"],
                     low_memory=False)
    ukb = ukb.loc[ukb["participant_id"].isin(set(analytic_ids()))].reset_index(drop=True)
    n_ukb = len(ukb)
    ukb_m = _normalised_exponential(n_ukb, rng)

    # ---- target ---------------------------------------------------------------
    pmr = read_table(_harmonised("pmr_with_ukb.csv"), usecols=["pmr_id", "w"],
                     low_memory=False)
    pmr = pmr.loc[pmr["pmr_id"].isin(set(pmr_analytic_ids()))].reset_index(drop=True)
    pmr_w = pmr["w"].to_numpy()
    pmr_m = _stratified_multiplier(pmr_w, rng, take_in_full=(1,))

    # ---- reference ------------------------------------------------------------
    reference = "census" if REFERENCE_SAMPLE == "census" else "hse"
    if reference == "hse":
        ref = read_table(_harmonised("hse.csv"), low_memory=False)
        if "psu" not in ref.columns:
            raise RuntimeError(
                "The harmonised reference sample has no `psu` column, so the "
                "perturbation cannot respect the survey design. Rebuild phase 1 "
                "stage 3 with the sampling unit carried through.")
        ref_m = _cluster_multiplier(ref["year"].to_numpy(), ref["psu"].to_numpy(), rng)
    else:
        ref = read_table(_harmonised("census.csv"), low_memory=False)
        ref_m = _normalised_exponential(len(ref), rng)

    for label, m in (("cohort", ukb_m), ("target", pmr_m), ("reference", ref_m)):
        if not np.isfinite(m).all():
            raise SystemExit(f"The {label} multipliers are not all finite.")
        if (m <= 0).any():
            raise SystemExit(
                f"{int((m <= 0).sum())} {label} multipliers are not strictly "
                "positive. Exp(1) cannot produce these, so the draw is wrong.")

    d = out_root / f"replicate_{replicate:04d}"
    d.mkdir(parents=True, exist_ok=True)
    np.save(d / "ukb_multiplicity.npy", ukb_m.astype(np.float64))
    np.save(d / "pmr_multiplicity.npy", pmr_m.astype(np.float64))
    write_table(pd.DataFrame({"participant_id": ukb["participant_id"].to_numpy(),
                              "multiplicity": ukb_m.astype(np.float64)}),
                d / "ukb_multiplicity.csv")

    ref_out = ref.copy()
    if reference == "hse":
        ref_out["weight_individual"] = (
            ref_out["weight_individual"].to_numpy(dtype=float) * ref_m)
    else:
        ref_out["bootstrap_multiplier"] = ref_m
    write_table(ref_out, d / f"{reference}.csv")

    manifest = {
        "replicate": int(replicate),
        "seed": _replicate_seed(replicate),
        "multiplier_distribution": "Exp(1), normalised to mean one within stratum",
        "cohort": {"n": int(n_ukb), "multipliers": _summary(ukb_m)},
        "target": {"n": int(len(pmr)),
                   "perturbed": int((pmr_w != 1).sum()),
                   "multipliers_sampled_strata": _summary(pmr_m[pmr_w != 1])},
        "reference": {"sample": reference, "n": int(len(ref)),
                      "sampling_units": (int(pd.unique(ref["psu"]).size)
                                         if reference == "hse" else int(len(ref))),
                      "multipliers": _summary(ref_m)},
    }
    with open(d / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    return manifest


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--replicate", type=int)
    ap.add_argument("--first", type=int)
    ap.add_argument("--last", type=int)
    args = ap.parse_args()

    out_root = path_temp(DRAW_SUBDIR, create=False)

    if args.replicate is not None:
        todo = [args.replicate]
    elif args.first is not None and args.last is not None:
        todo = list(range(args.first, args.last + 1))
    else:
        ap.error("give --replicate, or --first and --last")

    for b in todo:
        m = draw(b, out_root)
        u = m["cohort"]["multipliers"]
        p = m["target"]["multipliers_sampled_strata"]
        h = m["reference"]["multipliers"]
        print(f"replicate {b:4d}  "
              f"cohort mean {u['mean']:.4f} var {u['variance']:.4f} ess/n {u['ess_over_n']:.4f}  "
              f"target var {p['variance']:.4f}  reference var {h['variance']:.4f}")


if __name__ == "__main__":
    main()
