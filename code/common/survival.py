#!/usr/bin/env python3
"""Building blocks for the deep survival mortality-prediction pipeline."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn

# ────────────────────────────────────────────────────────────────────────────
# Settings
# ────────────────────────────────────────────────────────────────────────────

from paths import setting

T0 = np.datetime64(setting("phase1", "baseline_date"))
T_ADMIN_END = np.datetime64(setting("phase1", "follow_up_end"))

HORIZONS_YEARS = list(setting("phase1", "horizons_years"))
HORIZONS_DAYS = [int(round(365.25 * x)) for x in HORIZONS_YEARS]

RHS_NUMERIC = list(setting("phase1", "predictors", "outcome", "numeric"))
RHS_CATEGORICAL = list(setting("phase1", "predictors", "outcome", "categorical"))

RHS_COVARIATES = RHS_NUMERIC + RHS_CATEGORICAL

TARGET_CARRIED_COLUMNS = ["nssec_class8"]

CATEGORY_LEVEL_ORDERS = {
    "sex": ["1", "2"],
    "ethnicity": [str(x) for x in range(1, 7)],
    "tenure": [str(x) for x in range(1, 8)],
    "household_size": ["1", "2"],
    "econstatus": ["1", "2", "3", "4"],
    "education": ["1", "2", "3", "4"],
    "ruralurban": ["1", "2", "3", "4"],
    "health": ["1", "2", "3", "4"],
    "imd_decile": [str(x) for x in range(1, 11)],
}

OUTCOMES = [
    "all_cause_mortality",
    "cancer_mortality",
    "cardiovascular_mortality",
    "respiratory_mortality",
    "digestive_mortality",
    "other_mortality",
]

_declared = list(setting("phase2", "outcomes"))
if set(_declared) != set(OUTCOMES):
    raise RuntimeError(
        "framework_config.json declares the outcomes "
        + ", ".join(sorted(_declared)) + " and this module holds "
        + ", ".join(sorted(OUTCOMES)) + ". They must be the same six.")

COMPETING_CAUSES = [
    "cancer_mortality",
    "cardiovascular_mortality",
    "respiratory_mortality",
    "digestive_mortality",
    "other_mortality",
]

TRAINING_SOURCES = [
    "ukb",
    "ukbw_hse_superlearner",
    "ukbw_hse_lassologit",
    "ukbw_census_superlearner",
    "pmr",
    "ukbw_hse_raking",
    "ukbw_hse_entropy_balancing",
]

_declared_weights = ["ukbw_" + k for k in setting("phase2", "weight_sources")]
if sorted(_declared_weights) != sorted(s for s in TRAINING_SOURCES
                                       if s.startswith("ukbw_")):
    raise RuntimeError(
        "framework_config.json declares the weight sources "
        + ", ".join(sorted(_declared_weights)) + " and this module trains "
        + ", ".join(sorted(s for s in TRAINING_SOURCES if s.startswith("ukbw_")))
        + ". They must be the same set.")

def source_seed_term(source: str) -> int:
    return 100 * (TRAINING_SOURCES.index(source) + 1)

FOLDS = list(range(int(setting("phase2", "splits", "folds"))))
SEED = int(setting("phase2", "splits", "seed"))

TRAINING_SCHEDULE = dict(setting("phase2", "deep_surv"))

for _field, _only in (("baseline_estimator", "weighted_breslow"),
                      ("activation", "relu"),
                      ("optimizer", "adamw"),
                      ("ties", "efron")):
    if TRAINING_SCHEDULE[_field] != _only:
        raise RuntimeError(
            f"framework_config.json sets phase2.deep_surv.{_field} to "
            f"{TRAINING_SCHEDULE[_field]!r}. This pipeline implements "
            f"{_only!r} and nothing else.")

# ────────────────────────────────────────────────────────────────────────────
# Logging
# ────────────────────────────────────────────────────────────────────────────

def log_step(*parts):
    from datetime import datetime
    print(f"[{datetime.now().isoformat(timespec='seconds')}] " + "".join(str(p) for p in parts), flush=True)

# ────────────────────────────────────────────────────────────────────────────
# Coercion
# ────────────────────────────────────────────────────────────────────────────

def coerce_text_series(x):
    if pd.api.types.is_numeric_dtype(x):
        x = pd.to_numeric(x, errors="coerce").round().astype("Int64")
    s = x.astype("string").str.strip()
    s = s.replace({"": pd.NA, "NA": pd.NA, "nan": pd.NA, "<NA>": pd.NA})
    return s

def coerce_covariates(df):
    d = df.copy()
    if "age" in d.columns:
        d["age"] = pd.to_numeric(d["age"], errors="raise")
    for col in RHS_CATEGORICAL:
        if col == "imd_decile":
            v = pd.to_numeric(d[col], errors="raise")
            d[col] = v.round().astype("Int64").astype("string")
            d[col] = d[col].replace({"<NA>": pd.NA})
        else:
            d[col] = coerce_text_series(d[col])
    if "RGN11CD" in d.columns:
        d["RGN11CD"] = coerce_text_series(d["RGN11CD"])
    return d

# ────────────────────────────────────────────────────────────────────────────
# ICD matching for outcome construction
# ────────────────────────────────────────────────────────────────────────────

def normalise_icd_code(x):
    s = pd.Series(x, copy=False).astype("string").str.upper().str.strip()
    s = s.replace({"": pd.NA, "NA": pd.NA})
    s = s.str.replace(r"[^A-Z0-9]", "", regex=True)
    s = s.replace({"": pd.NA})
    return s

def normalise_icd_spec(x):
    s = str(x).strip().upper()
    s = "".join(ch for ch in s if ch.isalnum() or ch == "-")
    if not s:
        return None
    return s

def expand_icd_specs_for_matching(specs):
    specs = [normalise_icd_spec(x) for x in specs]
    specs = [x for x in specs if x]
    out = {"icd10_root3": set(), "icd10_exact4": set()}

    for spec in specs:
        if len(spec) == 3 and spec[0].isalpha() and spec[1:].isdigit():
            out["icd10_root3"].add(spec)
            continue
        if len(spec) == 7 and spec[0].isalpha() and spec[3] == "-" and spec[4].isalpha():
            left, right = spec.split("-")
            if left[0] != right[0]:
                raise RuntimeError(f"ICD-10 root ranges must share the same letter prefix: {spec}")
            lo = int(left[1:]); hi = int(right[1:])
            if lo > hi:
                raise RuntimeError(f"Invalid ICD-10 range: {spec}")
            for value in range(lo, hi + 1):
                out["icd10_root3"].add(f"{left[0]}{value:02d}")
            continue
        if len(spec) == 4 and spec[0].isalpha() and spec[1:].isdigit():
            out["icd10_exact4"].add(spec)
            continue
        if len(spec) == 9 and spec[0].isalpha() and spec[4] == "-" and spec[5].isalpha():
            left, right = spec.split("-")
            if left[0] != right[0]:
                raise RuntimeError(f"ICD-10 exact ranges must share the same letter prefix: {spec}")
            lo = int(left[1:]); hi = int(right[1:])
            if lo > hi:
                raise RuntimeError(f"Invalid ICD-10 exact range: {spec}")
            for value in range(lo, hi + 1):
                out["icd10_exact4"].add(f"{left[0]}{value:03d}")
            continue
        raise RuntimeError(f"Unsupported ICD spec: {spec}")

    return out

def icd_vector_matches_plan(codes, plan):
    norm = normalise_icd_code(codes)
    out = pd.Series(False, index=norm.index)

    if plan["icd10_root3"]:
        idx = norm.str.match(r"^[A-Z][0-9]{2}", na=False)
        out.loc[idx] = out.loc[idx] | norm.loc[idx].str.slice(0, 3).isin(plan["icd10_root3"])
    if plan["icd10_exact4"]:
        idx = norm.str.match(r"^[A-Z][0-9]{3}", na=False)
        out.loc[idx] = out.loc[idx] | norm.loc[idx].str.slice(0, 4).isin(plan["icd10_exact4"])
    return out.fillna(False)

def get_harmonised_death_icd_columns(df, match_scope):
    underlying = "underlying_cause_of_death_icd"
    if underlying not in df.columns:
        raise RuntimeError(f"Missing underlying cause-of-death column: {underlying}")
    if match_scope != "underlying":
        raise RuntimeError(f"Only 'underlying' match scope is supported in this pipeline; got {match_scope!r}")
    return [underlying]

def build_outcome_death_match(df, outcome_cfg):
    if outcome_cfg["type"] == "all_cause":
        return pd.Series(True, index=df.index)

    def match(specs):
        plan = expand_icd_specs_for_matching(list(specs))
        hit = pd.Series(False, index=df.index)
        for col in get_harmonised_death_icd_columns(df, outcome_cfg["match_scope"]):
            hit = hit | icd_vector_matches_plan(df[col], plan)
        return hit

    if outcome_cfg["type"] == "residual":
        named = pd.Series(False, index=df.index)
        for specs in outcome_cfg["excludes_icd10"]:
            named = named | match(specs)
        return ~named

    return match(outcome_cfg.get("icd10", []))

def build_survival_outcome(df, death_date_col, outcome_cfg):
    death_date = pd.to_datetime(df[death_date_col], errors="coerce").values.astype("datetime64[D]")
    matched_target = build_outcome_death_match(df, outcome_cfg).to_numpy(dtype=bool)
    death_any_in_followup = (~pd.isna(death_date)) & (death_date <= T_ADMIN_END)
    target_death_in_followup = death_any_in_followup & matched_target

    censor_date = np.full(len(df), T_ADMIN_END, dtype="datetime64[D]")
    has_date = ~pd.isna(death_date)
    censor_date[has_date] = np.minimum(death_date[has_date], T_ADMIN_END)

    time_days = (censor_date - T0).astype("timedelta64[D]").astype(int)

    return pd.DataFrame(
        {
            "censor_date": censor_date.astype("datetime64[D]"),
            "time_days": time_days.astype(int),
            "status": target_death_in_followup.astype(int),
            "death_any_in_followup": death_any_in_followup.astype(int),
            "target_death_in_followup": target_death_in_followup.astype(int),
        },
        index=df.index,
    )

# ────────────────────────────────────────────────────────────────────────────
# Common-support harmonisation
# ────────────────────────────────────────────────────────────────────────────

def assert_common_support(ukb_df, pmr_df, vars_):
    support_levels = {}
    problems = []

    for var in vars_:
        ukb_vals = set(pd.Series(ukb_df[var].dropna().astype(str)).unique().tolist())
        pmr_vals = set(pd.Series(pmr_df[var].dropna().astype(str)).unique().tolist())
        only_ukb = sorted(ukb_vals - pmr_vals)
        only_pmr = sorted(pmr_vals - ukb_vals)

        if only_ukb or only_pmr:
            n_ukb = int(ukb_df[var].astype(str).isin(only_ukb).sum()) if only_ukb else 0
            n_pmr = int(pmr_df[var].astype(str).isin(only_pmr).sum()) if only_pmr else 0
            problems.append(
                f"  {var}: training-only {only_ukb} ({n_ukb} rows), "
                f"evaluation-only {only_pmr} ({n_pmr} rows)"
            )
        support_levels[var] = sorted(ukb_vals & pmr_vals)

    if problems:
        raise RuntimeError(
            "Categorical support differs between the training and evaluation "
            "sources:\n" + "\n".join(problems) +
            "\nHarmonisation is expected to produce identical support. Correct "
            "the variable mapping rather than excluding these individuals: "
            "dropping them here would change the evaluation population after it "
            "has been fixed."
        )

    log_step(f"common support verified across {len(vars_)} categorical "
             "predictors; no individual excluded")
    return support_levels

def build_combined_levels(ukb_df, pmr_df):
    levels = {}
    for col in RHS_CATEGORICAL:
        present = set(pd.concat([ukb_df[col], pmr_df[col]], axis=0).dropna().astype(str).unique().tolist())
        ordered = [x for x in CATEGORY_LEVEL_ORDERS[col] if x in present]
        leftovers = sorted(present.difference(ordered))
        levels[col] = ordered + leftovers
    return levels

# ────────────────────────────────────────────────────────────────────────────
# Preprocessing (one-hot + standardisation)
# ────────────────────────────────────────────────────────────────────────────

@dataclass
class Preprocessor:
    numeric_features: list
    categorical_features: list
    categorical_levels: dict
    numeric_means: dict
    numeric_sds: dict
    feature_names: list

    def transform(self, df):
        parts = []
        feature_names = []
        for col in self.numeric_features:
            x = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
            mean = self.numeric_means[col]
            sd = self.numeric_sds[col]
            if not np.isfinite(sd) or sd <= 0:
                sd = 1.0
            z = (x - mean) / sd
            parts.append(z.reshape(-1, 1))
            feature_names.append(col)
        for col in self.categorical_features:
            vals = df[col].astype("string").to_numpy()
            for lev in self.categorical_levels[col]:
                parts.append((vals == lev).astype(float).reshape(-1, 1))
                feature_names.append(f"{col}__{lev}")
        X = np.hstack(parts).astype(np.float32) if parts else np.empty((len(df), 0), dtype=np.float32)
        return X, feature_names

def serialise_preprocessor(p: Preprocessor) -> dict:
    return {
        "numeric_features": list(p.numeric_features),
        "categorical_features": list(p.categorical_features),
        "categorical_levels": {k: list(v) for k, v in p.categorical_levels.items()},
        "numeric_means": {k: float(v) for k, v in p.numeric_means.items()},
        "numeric_sds": {k: float(v) for k, v in p.numeric_sds.items()},
        "feature_names": list(p.feature_names),
    }

def fit_preprocessor(train_df, combined_levels) -> Preprocessor:
    numeric_means, numeric_sds = {}, {}
    for col in RHS_NUMERIC:
        x = pd.to_numeric(train_df[col], errors="coerce").to_numpy(dtype=float)
        mean = float(np.nanmean(x))
        sd = float(np.nanstd(x))
        if not np.isfinite(sd) or sd <= 0:
            sd = 1.0
        numeric_means[col] = mean
        numeric_sds[col] = sd

    feature_names = list(RHS_NUMERIC)
    for col in RHS_CATEGORICAL:
        feature_names.extend([f"{col}__{lev}" for lev in combined_levels[col]])

    return Preprocessor(
        numeric_features=list(RHS_NUMERIC),
        categorical_features=list(RHS_CATEGORICAL),
        categorical_levels=combined_levels,
        numeric_means=numeric_means,
        numeric_sds=numeric_sds,
        feature_names=feature_names,
    )

# ────────────────────────────────────────────────────────────────────────────
# Network architecture, feed-forward, with hidden_dims read from the spec
# ────────────────────────────────────────────────────────────────────────────

class DeepSurvNet(nn.Module):
    """Feed-forward network with per-layer dropout and a linear-predictor output."""

    def __init__(self, input_dim, hidden_dims, dropout):
        super().__init__()
        layers = []
        prev_dim = int(input_dim)
        for hidden_dim in [int(x) for x in hidden_dims if int(x) > 0]:
            layers.append(nn.Linear(prev_dim, int(hidden_dim)))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(float(dropout)))
            prev_dim = int(hidden_dim)
        self.backbone = nn.Sequential(*layers)
        self.output = nn.Linear(prev_dim, 1)

    def forward(self, X):
        X = self.backbone(X)
        return self.output(X).squeeze(1)

def make_torch_model(input_dim, seed):
    """The network phase2.deep_surv declares, seeded by its caller."""
    torch.manual_seed(int(seed))
    np.random.seed(int(seed))
    return DeepSurvNet(
        input_dim=int(input_dim),
        hidden_dims=TRAINING_SCHEDULE["hidden_dims"],
        dropout=float(TRAINING_SCHEDULE["dropout"]),
    )

def tensor_from_numpy(x, dtype=torch.float32):
    return torch.from_numpy(np.asarray(x)).to(dtype=dtype)

def make_optimizer(model, lr, weight_decay):
    """The optimiser phase2.deep_surv declares."""
    return torch.optim.AdamW(model.parameters(), lr=float(lr),
                             weight_decay=float(weight_decay))

# ────────────────────────────────────────────────────────────────────────────
# Cox loss + design tensors + Breslow baseline
# ────────────────────────────────────────────────────────────────────────────

def _select_device():
    named = str(setting("phase2", "device")).strip().lower()
    if named == "cpu":
        return torch.device("cpu")
    if named != "gpu":
        raise SystemExit(f'phase2.device is "{named}", and it takes "gpu" or "cpu".')
    if not torch.cuda.is_available():
        raise SystemExit('phase2.device is "gpu" and this machine has no CUDA '
                         'device. Set it to "cpu" in framework_config.json.')
    return torch.device("cuda")

DEVICE = _select_device()

def _segment_sum(values: torch.Tensor, group_end_idx: torch.Tensor) -> torch.Tensor:
    cum = torch.cumsum(values, dim=0)[group_end_idx]
    prev = torch.cat([torch.zeros(1, dtype=cum.dtype, device=cum.device), cum[:-1]])
    return cum - prev

@dataclass
class CoxDesign:
    X_ord: torch.Tensor
    weight_ord: torch.Tensor
    event_weight_ord: torch.Tensor
    group_index: torch.Tensor
    group_end_idx: torch.Tensor
    group_event_weight: torch.Tensor
    event_group_mask: torch.Tensor
    total_event_weight: float
    event_times: np.ndarray

    def to(self, device) -> "CoxDesign":
        if device is None:
            return self
        moved = {f: getattr(self, f) for f in
                 ("total_event_weight", "event_times")}
        for f in ("X_ord", "weight_ord", "event_weight_ord", "group_index",
                  "group_end_idx", "group_event_weight", "event_group_mask"):
            moved[f] = getattr(self, f).to(device)
        return CoxDesign(**moved)

def build_cox_design(X, time_days, status, sample_weight) -> CoxDesign:
    time_days = np.asarray(time_days, dtype=float)
    status = np.asarray(status, dtype=int)
    sample_weight = np.asarray(sample_weight, dtype=np.float32)

    order = np.argsort(-time_days, kind="mergesort")
    time_ord = time_days[order]
    status_ord = status[order]
    weight_ord = sample_weight[order]
    X_ord = np.asarray(X, dtype=np.float32)[order]

    group_start = np.r_[True, time_ord[1:] != time_ord[:-1]]
    group_index = np.cumsum(group_start).astype(np.int64) - 1
    group_end_idx = np.where(np.r_[time_ord[1:] != time_ord[:-1], True])[0].astype(np.int64)
    event_weight_ord = (weight_ord * (status_ord == 1)).astype(np.float32)
    group_event_weight = np.bincount(group_index, weights=event_weight_ord, minlength=len(group_end_idx)).astype(np.float32)
    event_group_mask = group_event_weight > 0

    return CoxDesign(
        X_ord=tensor_from_numpy(X_ord, dtype=torch.float32),
        weight_ord=tensor_from_numpy(weight_ord, dtype=torch.float32),
        event_weight_ord=tensor_from_numpy(event_weight_ord, dtype=torch.float32),
        group_index=torch.from_numpy(group_index.astype(np.int64)),
        group_end_idx=torch.from_numpy(group_end_idx.astype(np.int64)),
        group_event_weight=tensor_from_numpy(group_event_weight, dtype=torch.float32),
        event_group_mask=torch.from_numpy(event_group_mask.astype(bool)),
        total_event_weight=float(np.sum(event_weight_ord)),
        event_times=time_ord[group_end_idx].astype(float),
    )

def weighted_cox_loss(model: nn.Module, design: CoxDesign) -> torch.Tensor:
    lp = model(design.X_ord)
    weights = design.weight_ord

    weighted_exp = weights * torch.exp(lp)

    cum_weighted_exp = torch.cumsum(weighted_exp, dim=0)
    S = cum_weighted_exp[design.group_end_idx]

    is_event = design.event_weight_ord > 0
    event_weighted_exp = torch.where(is_event, weighted_exp, torch.zeros_like(weighted_exp))
    T = _segment_sum(event_weighted_exp, design.group_end_idx)

    weighted_event_lp = design.event_weight_ord * lp
    group_event_eta = _segment_sum(weighted_event_lp, design.group_end_idx)
    d_g_f = _segment_sum(is_event.to(lp.dtype), design.group_end_idx)
    sum_w_event = design.group_event_weight

    mask = design.event_group_mask
    S_ev = S[mask]
    T_ev = T[mask]
    d_ev = d_g_f[mask]
    sum_w_ev = sum_w_event[mask]
    eta_w_ev = group_event_eta[mask]

    if d_ev.numel() == 0:
        return torch.zeros((), dtype=lp.dtype, device=lp.device)

    d_max = int(d_g_f.max().item())
    if d_max <= 0:
        return torch.zeros((), dtype=lp.dtype, device=lp.device)

    l_idx = torch.arange(d_max, dtype=lp.dtype, device=lp.device)
    l_over_d = l_idx.unsqueeze(0) / d_ev.unsqueeze(1)
    denom = S_ev.unsqueeze(1) - l_over_d * T_ev.unsqueeze(1)
    valid = l_idx.unsqueeze(0) < d_ev.unsqueeze(1)
    log_denom = torch.where(
        valid,
        torch.log(torch.clamp(denom, min=1e-12)),
        torch.zeros_like(denom),
    )
    sum_log_denom = log_denom.sum(dim=1)

    contrib = eta_w_ev - (sum_w_ev / d_ev) * sum_log_denom
    total = contrib.sum()
    denom_normaliser = max(design.total_event_weight, 1e-8)
    return -total / denom_normaliser

def build_weighted_breslow_baseline(model, design: CoxDesign):
    model.eval()
    design = design.to(next(model.parameters()).device)
    with torch.no_grad():
        lp = model(design.X_ord)
        log_weight = torch.log(torch.clamp(design.weight_ord, min=1e-8))
        log_risk_term = lp + log_weight
        log_cum_risk = torch.logcumsumexp(log_risk_term, dim=0)
        group_log_risk = log_cum_risk[design.group_end_idx]
        mask = design.event_group_mask
        event_times_desc = design.event_times[design.event_group_mask.cpu().numpy()].astype(float)
        delta_hazard_desc = design.group_event_weight[mask] / torch.exp(group_log_risk[mask])
        event_times = np.asarray(event_times_desc[::-1], dtype=float)
        delta_hazard = torch.flip(delta_hazard_desc, dims=[0])
        cum_hazard = torch.cumsum(delta_hazard, dim=0)
    return {"event_times": event_times,
            "cum_baseline_hazard": cum_hazard.cpu().numpy().astype(float)}

def predict_lp_torch(model, X, batch_size=32768):
    model.eval()
    device = next(model.parameters()).device
    X_t = tensor_from_numpy(X, dtype=torch.float32)
    out = []
    with torch.no_grad():
        for s in range(0, X_t.shape[0], batch_size):
            e = min(s + batch_size, X_t.shape[0])
            out.append(model(X_t[s:e].to(device)).cpu().numpy())
    return np.concatenate(out, axis=0) if out else np.empty((0,), dtype=np.float32)

def predict_risk_from_baseline(lp, baseline, horizons_days):
    lp = np.asarray(lp, dtype=float)
    event_times = np.asarray(baseline["event_times"], dtype=float)
    cum_hazard = np.asarray(baseline["cum_baseline_hazard"], dtype=float)
    out = {}
    for horizon in horizons_days:
        horizon = float(horizon)
        if len(event_times) == 0 or horizon <= 0:
            base_h = 0.0
        else:
            ix = np.searchsorted(event_times, horizon, side="right") - 1
            base_h = float(cum_hazard[ix]) if ix >= 0 else 0.0
        out[int(horizon)] = 1.0 - np.exp(-base_h * np.exp(lp))
    return out

# ────────────────────────────────────────────────────────────────────────────
# The permanent fold split, and the partitions taken within it
# ────────────────────────────────────────────────────────────────────────────

def stratified_kfold_indices(status, k, seed):
    rng = np.random.default_rng(int(seed))
    status = np.asarray(status, dtype=int)
    event_idx = np.where(status == 1)[0]
    nonevent_idx = np.where(status == 0)[0]
    rng.shuffle(event_idx)
    rng.shuffle(nonevent_idx)
    folds_event = np.array_split(event_idx, k)
    folds_nonevent = np.array_split(nonevent_idx, k)
    folds = [np.sort(np.concatenate([fe, fn])) for fe, fn in zip(folds_event, folds_nonevent)]
    return folds

def _stratified_carve(idx, status, fraction, rng):
    sub_status = np.asarray(status)[idx]
    event_pos = np.where(sub_status == 1)[0]
    nonevent_pos = np.where(sub_status == 0)[0]
    rng.shuffle(event_pos)
    rng.shuffle(nonevent_pos)
    n_event = int(round(len(event_pos) * float(fraction)))
    n_nonevent = int(round(len(nonevent_pos) * float(fraction)))
    carved_pos = np.concatenate([event_pos[:n_event], nonevent_pos[:n_nonevent]])
    rest_pos = np.concatenate([event_pos[n_event:], nonevent_pos[n_nonevent:]])
    return np.sort(idx[carved_pos]), np.sort(idx[rest_pos])

VALIDATION_FRACTION = float(setting("phase2", "splits", "validation_fraction"))
TEST_FRACTION = float(setting("phase2", "splits", "test_fraction_single_split"))
VALIDATION_FRACTION_OF_POOL = VALIDATION_FRACTION / (1.0 - 1.0 / len(FOLDS))

def stratified_train_val_test_split_within_fold(status, fold_test_idx, seed):
    val_fraction_of_pool = VALIDATION_FRACTION_OF_POOL
    rng = np.random.default_rng(int(seed))
    status = np.asarray(status, dtype=int)
    n = len(status)
    test_idx = np.sort(np.asarray(fold_test_idx, dtype=int))

    pool_mask = np.ones(n, dtype=bool)
    pool_mask[test_idx] = False
    pool_idx = np.where(pool_mask)[0]

    val_idx, train_idx = _stratified_carve(pool_idx, status, val_fraction_of_pool, rng)
    return train_idx, val_idx, test_idx

def stratified_train_val_test_split(status, seed):
    val_fraction = VALIDATION_FRACTION
    test_fraction = TEST_FRACTION
    rng = np.random.default_rng(int(seed))
    status = np.asarray(status, dtype=int)
    all_idx = np.arange(len(status))

    test_idx, pool_idx = _stratified_carve(all_idx, status, test_fraction, rng)
    val_share = float(val_fraction) / max(1e-12, 1.0 - float(test_fraction))
    val_idx, train_idx = _stratified_carve(pool_idx, status, val_share, rng)
    return train_idx, val_idx, test_idx

# ────────────────────────────────────────────────────────────────────────────
# Training with a staged learning-rate decay
# ────────────────────────────────────────────────────────────────────────────

def train_with_lr_schedule(model, train_design, val_design, seed):
    sched = TRAINING_SCHEDULE
    torch.manual_seed(int(seed))
    np.random.seed(int(seed))

    device = DEVICE
    if device.type == "cpu":
        torch.set_num_threads(1)
    model = model.to(device)
    train_design = train_design.to(device)
    if val_design is not None:
        val_design = val_design.to(device)

    base_lr = float(sched["base_learning_rate"])
    lr_floor = float(sched["lr_floor"])
    lr_decay = float(sched["lr_decay_factor"])
    epochs_per_stage = int(sched["epochs_per_stage"])
    patience = int(sched["early_stop_patience"])
    grad_clip = float(sched["grad_clip_max_norm"])
    weight_decay = float(sched["weight_decay"])

    history = {"stage": [], "lr": [], "epoch": [], "train_loss": [], "val_loss": []}
    best_val = float("inf")
    best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}

    lr = base_lr
    stage = 0
    while True:
        optimizer = make_optimizer(model, lr=lr, weight_decay=weight_decay)
        no_improve = 0
        for epoch in range(1, epochs_per_stage + 1):
            model.train()
            optimizer.zero_grad(set_to_none=True)
            loss = weighted_cox_loss(model, train_design)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=grad_clip)
            optimizer.step()

            train_loss = float(loss.detach().cpu().item())
            if val_design is not None:
                model.eval()
                with torch.no_grad():
                    vloss = float(weighted_cox_loss(model, val_design).detach().cpu().item())
            else:
                vloss = train_loss

            history["stage"].append(int(stage))
            history["lr"].append(float(lr))
            history["epoch"].append(int(epoch))
            history["train_loss"].append(train_loss)
            history["val_loss"].append(vloss)

            if vloss + 1e-6 < best_val:
                best_val = vloss
                best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
                no_improve = 0
            else:
                no_improve += 1
                if no_improve >= patience:
                    break

        stage += 1
        if lr <= lr_floor + 1e-12:
            break
        lr = max(lr / lr_decay, lr_floor)

    model.load_state_dict(best_state)
    return model, history, best_val

# ────────────────────────────────────────────────────────────────────────────
# Path / config helpers
# ────────────────────────────────────────────────────────────────────────────

def find_framework_root(start_path: Path) -> Path:
    from paths import repo_root
    return repo_root(start_path)

def read_framework_config(start_path: Path) -> dict:
    root = find_framework_root(start_path)
    with open(root / "framework_config.json") as f:
        cfg = json.load(f)
    cfg["framework_root"] = str(root)
    return cfg

def get_outcome_cfg_from_framework(cfg: dict, key: str) -> dict:
    outcomes = cfg["phase2"]["outcomes"]
    if key not in outcomes:
        raise RuntimeError(f"Outcome {key!r} not in framework_config.json")
    raw = outcomes[key]
    out = {
        "key": key,
        "label": raw["label"],
        "type": raw["type"],
        "icd10": list(raw.get("icd10", []) or []),
    }
    if out["type"] != "all_cause":
        out["match_scope"] = raw["match_scope"]
    if out["type"] == "residual":
        excludes = raw.get("excludes_icd10")
        if not excludes:
            raise RuntimeError(
                f"Outcome {key!r} is residual but names no causes to exclude. "
                "Without them it would match every death."
            )
        out["excludes_icd10"] = [list(x) for x in excludes]
    return out

def get_weight_source_from_framework(cfg: dict, key: str) -> dict:
    sources = cfg["phase2"]["weight_sources"]
    if key not in sources:
        raise RuntimeError(f"Weight source {key!r} not in framework_config.json")
    return {"key": key, **sources[key]}

def get_model_dir(cause: str, source: str, fold=None) -> Path:
    from paths import path_temp
    leaf = "single" if fold is None else f"fold_{int(fold)}"
    out = path_temp("models", cause, source, leaf, create=False)
    out.mkdir(parents=True, exist_ok=True)
    return out

def risk_column(cause: str, source: str, horizon_years: int) -> str:
    if cause in COMPETING_CAUSES:
        return f"{source}_cif_{horizon_years}y"
    return f"{source}_pred_{horizon_years}y"



