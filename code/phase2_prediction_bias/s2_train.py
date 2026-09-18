#!/usr/bin/env python3
"""Fit one deep survival model and score PMR with it.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "common"))

import folds as fold_split
import survival as dl
from cohort import attach_ukb_weight, prepare_data

OUTCOMES = dl.OUTCOMES
TRAINING_SOURCES = dl.TRAINING_SOURCES
FOLDS = dl.FOLDS
HORIZONS_DAYS = dl.HORIZONS_DAYS

def _resolve_seed(cause: str, source: str, fold: int | None) -> int:
    return (dl.SEED
            + 1000 * (OUTCOMES.index(cause) + 1)
            + dl.source_seed_term(source)
            + (0 if fold is None else int(fold)))


def train_task(cause: str, source: str, fold: int | None):
    cross_fitted = (source == "pmr")
    if cross_fitted and fold is None:
        raise RuntimeError("--fold is required for source=pmr, which is cross-fitted.")
    if not cross_fitted and fold is not None:
        raise RuntimeError(
            f"--fold is not accepted for source={source}. UKB sources train "
            "a single model because they are disjoint from the evaluation cohort."
        )

    seed = _resolve_seed(cause, source, fold)
    dl.log_step(f"=== train: cause={cause} source={source} fold={fold} seed={seed} ===")

    spec = dl.TRAINING_SCHEDULE
    dl.log_step(f"specification: hidden_dims={spec['hidden_dims']} "
                f"dropout={spec['dropout']}")

    data = prepare_data(cause)
    ukb, pmr, levels = data["ukb"], data["pmr"], data["levels"]

    if cross_fitted:
        train_df = pmr
        fold_per_row = fold_split.read(cause, pmr["pmr_id"].astype(int))
        status = train_df["status"].to_numpy(dtype=int)
        train_idx, val_idx, test_idx = dl.stratified_train_val_test_split_within_fold(
            status=status, fold_test_idx=np.where(fold_per_row == int(fold))[0], seed=seed)

    else:
        train_df = attach_ukb_weight(ukb, source)
        status = train_df["status"].to_numpy(dtype=int)
        train_idx, val_idx, test_idx = dl.stratified_train_val_test_split(
            status=status, seed=seed)

    dl.log_step(f"split sizes: train={len(train_idx)} val={len(val_idx)} test={len(test_idx)}"
                f" (events in train: {int(status[train_idx].sum())})")

    preprocessor = dl.fit_preprocessor(train_df.iloc[train_idx], levels)
    X, _ = preprocessor.transform(train_df)
    time_days = train_df["time_days"].to_numpy(dtype=float)
    weight = train_df["w"].to_numpy(dtype=np.float32)

    _ev = status == 1
    _wev = weight[train_idx][_ev[train_idx]].astype(float)
    _mass, _sq = float(_wev.sum()), float((_wev ** 2).sum())
    _ess = _mass ** 2 / _sq if _sq > 0 else 0.0
    dl.log_step(f"weighted events in train: {len(_wev)} carrying an effective "
                f"{_ess:.1f}")
    if not (np.isfinite(_mass) and _mass > 0 and _sq > 0):
        raise SystemExit(
            f"{cause} {source}: the {len(_wev)} events in the training "
            f"partition carry a total weight of {_mass:.3g}. The weighted "
            "partial likelihood has nothing to fit. This is the weight set, "
            "not the model: check the effective sample size in "
            "results/phase3_weighting/weight_summary.xlsx.")

    train_design = dl.build_cox_design(X[train_idx], time_days[train_idx],
                                       status[train_idx], weight[train_idx])
    val_design = (dl.build_cox_design(X[val_idx], time_days[val_idx],
                                      status[val_idx], weight[val_idx])
                  if len(val_idx) else None)

    model = dl.make_torch_model(X.shape[1], seed=seed)
    t0 = time.time()
    model, history, best_val = dl.train_with_lr_schedule(
        model, train_design, val_design, seed=seed)
    elapsed = time.time() - t0
    dl.log_step(f"trained: best_val_loss={best_val:.6f} elapsed={elapsed:.1f}s")

    baseline = dl.build_weighted_breslow_baseline(model, train_design)

    scored = pmr.iloc[test_idx] if cross_fitted else pmr
    X_pmr, _ = preprocessor.transform(scored)
    lp = dl.predict_lp_torch(model, X_pmr)
    risk = dl.predict_risk_from_baseline(lp, baseline, HORIZONS_DAYS)
    dl.log_step(f"scored {len(scored)} PMR rows "
                f"({'own fold only' if cross_fitted else 'full cohort'})")

    out_dir = dl.get_model_dir(cause, source, fold)
    torch.save({k: v.detach().cpu() for k, v in model.state_dict().items()},
               out_dir / "model_state.pt")
    with open(out_dir / "preprocessor.json", "w") as f:
        json.dump(dl.serialise_preprocessor(preprocessor), f, indent=2)
    np.savez_compressed(out_dir / "baseline_hazard.npz",
                        event_times=baseline["event_times"],
                        cum_baseline_hazard=baseline["cum_baseline_hazard"])
    np.savez_compressed(out_dir / "pmr_predictions.npz",
                        pmr_id=scored["pmr_id"].astype(int).to_numpy(),
                        lp=np.asarray(lp, dtype=float),
                        **{f"risk_{y}y": risk[d].astype(float)
                           for y, d in zip(dl.HORIZONS_YEARS, HORIZONS_DAYS)})
    with open(out_dir / "test_metrics.json", "w") as f:
        json.dump({
            "cause": cause, "source": source,
            "fold": None if fold is None else int(fold),
            "cross_fitted": bool(cross_fitted), "seed": int(seed),
            "n_train": int(len(train_idx)), "n_val": int(len(val_idx)),
            "n_test": int(len(test_idx)),
            "n_test_events": int(status[test_idx].sum()),
            "best_val_loss": float(best_val),
            "training_seconds": float(elapsed),
            "input_dim": int(X.shape[1]),
            "hidden_dims": list(spec["hidden_dims"]),
            "dropout": float(spec["dropout"]),
            "schedule": spec,
        }, f, indent=2)
    with open(out_dir / "history.json", "w") as f:
        json.dump(history, f, indent=2)
    dl.log_step(f"wrote {out_dir}")


def _parse_tasks(spec: str) -> list:
    """Read "cause:source[:fold]" items, comma separated, into task tuples.
    """
    tasks = []
    for item in spec.split(","):
        item = item.strip()
        if not item:
            continue
        parts = item.split(":")
        if len(parts) not in (2, 3):
            raise SystemExit(f"Task {item!r} is not cause:source or cause:source:fold")
        cause, source = parts[0], parts[1]
        fold = int(parts[2]) if len(parts) == 3 else None
        if cause not in OUTCOMES:
            raise SystemExit(f"Task {item!r} names an unknown cause {cause!r}")
        if source not in TRAINING_SOURCES:
            raise SystemExit(f"Task {item!r} names an unknown source {source!r}")
        if fold is not None and fold not in FOLDS:
            raise SystemExit(f"Task {item!r} names an unknown fold {fold}")
        tasks.append((cause, source, fold))
    if not tasks:
        raise SystemExit("--tasks is empty")
    return tasks


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cause", choices=OUTCOMES)
    ap.add_argument("--source", choices=TRAINING_SOURCES)
    ap.add_argument("--fold", type=int, choices=FOLDS,
                    help="required for --source pmr, rejected otherwise")
    ap.add_argument("--tasks", default=None,
                    help='comma-separated "cause:source" or "cause:source:fold"')
    args = ap.parse_args()

    if args.tasks:
        if args.cause or args.source or args.fold is not None:
            ap.error("--tasks replaces --cause, --source and --fold")
        tasks = _parse_tasks(args.tasks)
        dl.log_step(f"=== {len(tasks)} tasks in one process ===")
        for i, (cause, source, fold) in enumerate(tasks, start=1):
            dl.log_step(f"=== task {i} of {len(tasks)} ===")
            train_task(cause, source, fold)
        return

    if not args.cause or not args.source:
        ap.error("give --cause and --source, or --tasks")
    train_task(args.cause, args.source, args.fold)


if __name__ == "__main__":
    main()
