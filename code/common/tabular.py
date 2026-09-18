"""Reading and writing the pipeline's tables. Mirrors code/common/tabular.R.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

_ID_COLUMNS = ("eid", "participant_id", "pmr_id", "ref_id")
_NA_VALUES = ["NA"]
_ISO_DATE = re.compile(r"^\d{4}-\d{2}-\d{2}$")


def _restore_dates(frame: pd.DataFrame) -> pd.DataFrame:
    for col in frame.columns:
        if frame[col].dtype != object:
            continue
        seen = frame[col].dropna()
        if seen.empty or not all(isinstance(v, str) for v in seen.head(1000)):
            continue
        if seen.map(lambda v: bool(_ISO_DATE.match(str(v)))).all():
            frame[col] = pd.to_datetime(frame[col], format="%Y-%m-%d")
    return frame


def read_table(path, usecols=None, nrows=None, dtype=None, **kw):
    kw.pop("na", None)
    kw.pop("low_memory", None)
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"No table at {p}")
    if p.suffix == ".xlsx":
        frame = pd.read_excel(p, usecols=list(usecols) if usecols else None,
                              nrows=int(nrows) if nrows is not None else None,
                              dtype=dtype, na_values=_NA_VALUES,
                              keep_default_na=False)
        return frame
    frame = pd.read_csv(
        p,
        usecols=list(usecols) if usecols else None,
        nrows=int(nrows) if nrows is not None else None,
        dtype=dtype,
        na_values=_NA_VALUES,
        keep_default_na=False,
        low_memory=False,
        **kw,
    )
    if dtype is None:
        frame = _restore_dates(frame)
    widen = {c: "int64" for c in _ID_COLUMNS
             if c in frame.columns and str(frame[c].dtype).startswith("int")}
    if widen:
        frame = frame.astype(widen)
    return frame


def write_table(frame: pd.DataFrame, path, index: bool = False) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    out = frame if index else frame.reset_index(drop=True)
    if p.suffix == ".xlsx":
        out.to_excel(p, index=index, na_rep="NA")
    else:
        out.to_csv(p, index=index, na_rep="NA", date_format="%Y-%m-%d")
    return p
