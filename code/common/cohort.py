"""The rows every model is fitted to, and the outcome built on them.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tabular import read_table

import survival as sv
from paths import analytic_ids, path_temp, pmr_analytic_ids, repo_root, setting

_ICD_COL = "underlying_cause_of_death_icd"
_UKB_COLUMNS = (["participant_id", "date_of_death", _ICD_COL, "RGN11CD"]
                + sv.RHS_COVARIATES)
_PMR_COLUMNS = (["pmr_id", "w", "dod_deaths", _ICD_COL, "RGN11CD"]
                + sv.RHS_COVARIATES + sv.TARGET_CARRIED_COLUMNS)


_HARMONISED_CACHE: dict = {}


def _read_harmonised(name: str, columns: list[str]) -> pd.DataFrame:
    key = (str(path_temp("harmonised", name, create=False)), tuple(columns))
    hit = _HARMONISED_CACHE.get(key)
    if hit is not None:
        return hit.copy(deep=True)
    frame = _read_harmonised_uncached(name, columns)
    _HARMONISED_CACHE[key] = frame
    return frame.copy(deep=True)


def _read_harmonised_uncached(name: str, columns: list[str]) -> pd.DataFrame:
    path = path_temp("harmonised", name, create=False)
    if not path.exists():
        raise RuntimeError(f"Missing {path}. Run phase 1 first.")
    have = read_table(path, nrows=0).columns
    missing = [c for c in columns if c not in have]
    if missing:
        raise RuntimeError(f"{path} is missing {missing}. The harmonisation has changed.")
    return read_table(path, usecols=columns, low_memory=False)


def _restrict(frame: pd.DataFrame, id_col: str, keep: set[int], label: str) -> pd.DataFrame:
    before = len(frame)
    out = frame.loc[frame[id_col].astype(int).isin(keep)].reset_index(drop=True)
    sv.log_step(f"analytic population: {len(out)} of {before} {label} rows retained")
    if len(out) != len(keep):
        raise RuntimeError(
            f"The {label} analytic population holds {len(keep)} records but only "
            f"{len(out)} were found in the harmonised file. The population and the "
            "harmonised data are out of step. Rerun phase 1."
        )
    return out


def _apply_bootstrap_multipliers(ukb, pmr):
    directory = os.environ.get("REWEIGHTING_BOOTSTRAP_DIR", "").strip()
    if not directory:
        return ukb, pmr

    root = Path(directory)
    out = []
    for frame, name in ((ukb, "ukb"), (pmr, "pmr")):
        path = root / f"{name}_multiplicity.npy"
        if not path.exists():
            raise RuntimeError(
                f"REWEIGHTING_BOOTSTRAP_DIR is set to {root} but {path.name} is "
                "not there. Draw the replicate before running any stage of it."
            )
        m = np.load(path)
        if len(m) != len(frame):
            raise RuntimeError(
                f"{path.name} has {len(m)} multipliers against a "
                f"{name.upper()} frame of {len(frame)}. The replicate was drawn "
                "against a different population."
            )
        if not (m > 0).all():
            raise RuntimeError(
                f"{path.name} holds a multiplier that is not strictly positive. "
                "A row weighing nothing is not inert: a large share of them "
                "makes one learner of the participation ensemble return nothing."
            )
        drawn = frame.copy()
        drawn["bootstrap_multiplicity"] = m.astype(float)
        if "w" in drawn.columns:
            drawn["w"] = drawn["w"].to_numpy(dtype=float) * m.astype(float)
        out.append(drawn.reset_index(drop=True))
        sv.log_step(f"bootstrap replicate {root.name}: {name.upper()} "
                    f"{len(m)} multipliers applied, largest "
                    f"{float(m.max(initial=0)):.4g}")

    return out[0], out[1]


_PREPARED_CACHE: dict = {}


def prepare_data(cause: str) -> dict:
    key = (cause, os.environ.get("REWEIGHTING_BOOTSTRAP_DIR", "").strip())
    hit = _PREPARED_CACHE.get(key)
    if hit is not None:
        return {"ukb": hit["ukb"].copy(deep=True),
                "pmr": hit["pmr"].copy(deep=True),
                "levels": hit["levels"],
                "outcome_cfg": hit["outcome_cfg"]}
    out = _prepare_data_uncached(cause)
    _PREPARED_CACHE[key] = out
    return {"ukb": out["ukb"].copy(deep=True),
            "pmr": out["pmr"].copy(deep=True),
            "levels": out["levels"],
            "outcome_cfg": out["outcome_cfg"]}


def _prepare_data_uncached(cause: str) -> dict:
    cfg = sv.read_framework_config(repo_root())
    outcome_cfg = sv.get_outcome_cfg_from_framework(cfg, cause)

    ukb_raw = sv.coerce_covariates(_read_harmonised("ukb_with_pmr.csv", _UKB_COLUMNS))
    pmr_raw = sv.coerce_covariates(_read_harmonised("pmr_with_ukb.csv", _PMR_COLUMNS))

    ukb_outcome = sv.build_survival_outcome(ukb_raw, "date_of_death", outcome_cfg)
    pmr_outcome = sv.build_survival_outcome(pmr_raw, "dod_deaths", outcome_cfg)
    ukb = pd.concat([ukb_raw.reset_index(drop=True), ukb_outcome.reset_index(drop=True)], axis=1)
    pmr = pd.concat([pmr_raw.reset_index(drop=True), pmr_outcome.reset_index(drop=True)], axis=1)

    ukb = _restrict(ukb, "participant_id", analytic_ids(), "UKB")
    pmr = _restrict(pmr, "pmr_id", pmr_analytic_ids(), "PMR")

    ukb, pmr = _apply_bootstrap_multipliers(ukb, pmr)

    sv.assert_common_support(ukb, pmr, sv.RHS_CATEGORICAL)

    for name, frame in (("UKB", ukb), ("PMR", pmr)):
        n_bad = int((frame["time_days"] < 0).sum())
        if n_bad:
            raise RuntimeError(
                f"{n_bad} {name} records have a death date before the baseline. "
                "The harmonisation stage is expected to exclude these."
            )

    for hz in sv.HORIZONS_DAYS:
        suffix = f"{int(round(hz / 365.25))}y"
        for d in (ukb, pmr):
            d[f"event_{suffix}"] = ((d["status"].to_numpy() == 1) &
                                    (d["time_days"].to_numpy(dtype=float) <= float(hz))).astype(int)

    if "w" not in pmr.columns:
        raise RuntimeError("PMR data is missing the 'w' (design weight) column.")

    levels = sv.build_combined_levels(ukb, pmr)
    return {"ukb": ukb, "pmr": pmr, "levels": levels, "outcome_cfg": outcome_cfg}


WINSORISE_PROBS = tuple(setting("phase3", "winsorise_percentiles"))


def _winsorise(w: np.ndarray) -> np.ndarray:
    lo, hi = np.quantile(w, WINSORISE_PROBS)
    out = np.clip(w, lo, hi)
    return out / out.mean()


def _apply_pending_winsorisation(weights_dir: Path, w: np.ndarray) -> np.ndarray:
    marker = weights_dir / "winsorise.txt"
    if not marker.exists() or marker.read_text().strip().splitlines()[:1] != ["1"]:
        return w
    return _winsorise(w)


def _apply_bootstrap_multiplier(out: pd.DataFrame) -> pd.DataFrame:
    if "bootstrap_multiplicity" in out.columns:
        out["w"] = out["w"].to_numpy(dtype=float) * \
            out["bootstrap_multiplicity"].to_numpy(dtype=float)
    return out


def attach_ukb_weight(ukb: pd.DataFrame, source: str) -> pd.DataFrame:
    out = ukb.copy()
    if source == "ukb":
        out["w"] = 1.0
        return _apply_bootstrap_multiplier(out)
    if not source.startswith("ukbw_"):
        raise RuntimeError(f"Unexpected UKB-derived source: {source}")

    key = source[len("ukbw_"):]
    keep = analytic_ids()
    cfg = sv.read_framework_config(repo_root())
    src = sv.get_weight_source_from_framework(cfg, key)
    path = path_temp("weights", key, "ukb_weights.csv", create=False)
    if not path.exists():
        raise RuntimeError(f"Missing weight file for {source}: {path}. Run phase 3 first.")
    wts = read_table(path)
    if src["column"] not in wts.columns:
        raise RuntimeError(f"Weight column {src['column']!r} not in {path}")
    wts = wts[["eid", src["column"]]].rename(columns={src["column"]: "w"})
    wts["eid"] = pd.to_numeric(wts["eid"], errors="raise").astype(int)

    out = out.merge(wts, left_on="participant_id", right_on="eid", how="left")
    out["w"] = pd.to_numeric(out["w"], errors="raise")

    missing = out.loc[out["w"].isna(), "participant_id"]
    if len(missing):
        raise RuntimeError(
            f"{len(missing)} of {len(out)} participants have no {source} weight, "
            f"for example {list(missing.head(5))}. Weights must cover the "
            "analytic population exactly. Re-estimate them in phase 3."
        )
    extra = set(wts["eid"]) - set(keep)
    if extra:
        raise RuntimeError(
            f"The {source} weight file covers {len(extra)} participants outside the "
            f"analytic population, for example {sorted(extra)[:5]}. The "
            "participation model was estimated on a different population from the "
            "one it is applied to. Re-estimate it in phase 3."
        )

    w = out["w"].to_numpy(dtype=float)
    winsorised = _apply_pending_winsorisation(path.parent, w)
    if winsorised is not w:
        lo, hi = (100 * p for p in WINSORISE_PROBS)
        sv.log_step(f"{source}: weights winsorised at the {lo:g} and {hi:g} "
                    f"percentiles on use, largest {w.max():.4g} to "
                    f"{winsorised.max():.4g}")
        out["w"] = winsorised

    return _apply_bootstrap_multiplier(out)
