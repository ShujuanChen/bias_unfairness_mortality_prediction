"""Path resolution for data/, temp/ and results/. Mirrors code/common/paths.R.
"""

from __future__ import annotations

import os
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tabular import read_table


def repo_root(start: Path | None = None) -> Path:
    p = Path(start or __file__).resolve()
    if p.is_file():
        p = p.parent
    for candidate in [p, *p.parents]:
        if (candidate / "framework_config.json").is_file():
            return candidate
    raise RuntimeError(f"framework_config.json not found above {p}")


_SETTINGS: dict | None = None


def setting(*keys: str):
    """One declared setting, by its path through framework_config.json."""
    global _SETTINGS
    if _SETTINGS is None:
        with open(repo_root() / "framework_config.json") as fh:
            _SETTINGS = json.load(fh)
    node = _SETTINGS
    for key in keys:
        if not isinstance(node, dict) or key not in node:
            raise RuntimeError("framework_config.json is missing " + ".".join(keys))
        node = node[key]
    return node


def run_id() -> str:
    rid = os.environ.get("REWEIGHTING_RUN_ID", "").strip()
    if not rid:
        raise RuntimeError(
            "No run identifier. Stage scripts are launched by run_pipeline.R, "
            "which stamps each invocation with a timestamp. To run a stage "
            "alone, set REWEIGHTING_RUN_ID to the run you want to write into.")
    return rid


def temp_root(root: Path | None = None) -> Path:
    override = os.environ.get("REWEIGHTING_TEMP_ROOT", "").strip()
    if override:
        return Path(override)
    return (root or repo_root()) / "temp" / run_id()


def path_temp(*parts: str, root: Path | None = None, create: bool = True) -> Path:
    """Intermediates. Parent directories created on demand."""
    p = temp_root(root).joinpath(*parts)
    if create:
        p.parent.mkdir(parents=True, exist_ok=True)
    return p


def path_results(*parts: str, root: Path | None = None, create: bool = True) -> Path:
    p = (root or repo_root()).joinpath("results", run_id(), *parts)
    if create:
        p.parent.mkdir(parents=True, exist_ok=True)
    return p


def analytic_ids() -> set[int]:
    path = path_temp("population", "ukb_population_membership.csv", create=False)
    if not path.exists():
        raise RuntimeError(
            f"Missing the analytic population at {path}. Build it with phase 1."
        )
    frame = read_table(path, usecols=["participant_id", "in_unified"])
    flag = frame["in_unified"]
    if flag.dtype == object or str(flag.dtype) == "string":
        flag = flag.astype(str).str.strip().str.upper() == "TRUE"
    return {int(v) for v in frame.loc[flag.astype(bool), "participant_id"]}


def require_pipeline(name: str, what: str) -> None:
    path = path_temp("population", "pipelines.txt", create=False)
    if not path.exists():
        raise RuntimeError(
            "No population/pipelines.txt in this run. Run phase 1 stage 5 first.")
    built = [line.strip() for line in path.read_text().split("\n") if line.strip()]
    if name not in built:
        raise RuntimeError(
            f"{what} needs the {name} pipeline, but this run's population was "
            f"built for {', '.join(built)}. Rebuild phase 1 with --pipelines "
            f"including {name}.")


def pmr_analytic_ids() -> set[int]:
    path = path_temp("population", "pmr_analytic_ids.csv", create=False)
    if not path.exists():
        raise RuntimeError(
            f"Missing the PMR analytic population at {path}. Build it with phase 1."
        )
    return {int(v) for v in read_table(path, usecols=["pmr_id"])["pmr_id"]}
