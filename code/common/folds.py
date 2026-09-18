"""The PMR fold assignment, shared by every model family.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tabular import read_table, write_table

import survival as sv
from paths import path_temp

GROUP = "pmr"


def assignment_path(cause: str, create: bool = True) -> Path:
    return path_temp("folds", cause, f"{GROUP}.csv", create=create)


MAX_SPLIT_ATTEMPTS = 20


def _levels_confined_to_one_fold(frame, fold_id):
    confined = []
    for col in sv.RHS_CATEGORICAL:
        if col not in frame.columns:
            continue
        values = frame[col].astype(str)
        for level in values.dropna().unique():
            rows = values == level
            n_folds = len(np.unique(fold_id[rows.to_numpy()]))
            if n_folds < 2:
                confined.append((col, level, int(rows.sum()), n_folds))
    return confined


def build(cause: str, ids, status, frame=None) -> pd.DataFrame:
    base = sv.SEED + 1000 * (sv.OUTCOMES.index(cause) + 1) + 1
    fold_id = None
    for attempt in range(MAX_SPLIT_ATTEMPTS):
        seed = base + attempt
        blocks = sv.stratified_kfold_indices(status=status,
                                             k=len(sv.FOLDS), seed=seed)
        fold_id = np.empty(len(ids), dtype=int)
        for k, idx in enumerate(blocks):
            fold_id[idx] = k
        if frame is None:
            break
        confined = _levels_confined_to_one_fold(frame, fold_id)
        if not confined:
            if attempt:
                sv.log_step(f"fold assignment for {cause} redrawn {attempt} time(s): "
                            "the earlier draws confined a level to one fold")
            break
        impossible = [c for c in confined if c[2] < 2]
        if impossible:
            raise RuntimeError(
                f"No five-fold split of {cause} can give every level a training "
                "example, because these are carried by fewer than two "
                "individuals:\n"
                + "\n".join(f"  {v} = {lev}: {n} row(s)" for v, lev, n, _ in impossible)
                + "\nThis is a property of the data rather than of the split. "
                "Phase 1 fixes common support, so a level this rare should have "
                "been excluded there."
            )
    else:
        raise RuntimeError(
            f"{MAX_SPLIT_ATTEMPTS} stratified splits of {cause} all confined a "
            "categorical level to one fold. The last attempt left:\n"
            + "\n".join(f"  {v} = {lev}: {n} rows in {f} fold"
                         for v, lev, n, f in confined)
        )

    out = pd.DataFrame({"id": np.asarray(ids).astype(int), "fold": fold_id})
    path = assignment_path(cause)
    write_table(out, path)
    sv.log_step(f"wrote {path} (seed={seed}, "
                + ", ".join(f"fold {k}: {int((fold_id == k).sum())}" for k in sv.FOLDS) + ")")
    return out


def read(cause: str, ids) -> np.ndarray:
    path = assignment_path(cause, create=False)
    if not path.exists():
        raise RuntimeError(
            f"No fold assignment at {path}. Build it with phase 2 stage 1:\n"
            f"    python code/phase2_prediction_bias/s1_assign_folds.py --cause {cause}"
        )
    table = read_table(path)
    ids = pd.Series(np.asarray(ids)).astype(int)
    lookup = dict(zip(table["id"].astype(int), table["fold"].astype(int)))
    unknown = ids[~ids.isin(lookup)]
    if len(unknown):
        raise RuntimeError(
            f"{len(unknown)} PMR identifiers are absent from {path}; it was built "
            "on a different cohort. Rerun phase 2 stage 1."
        )
    if len(lookup) != len(ids):
        raise RuntimeError(
            f"{path} assigns {len(lookup)} individuals against a cohort of "
            f"{len(ids)}. Rerun phase 2 stage 1."
        )
    return ids.map(lookup).to_numpy(dtype=int)
