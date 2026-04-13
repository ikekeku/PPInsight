"""
provenance – run-level metadata for docking score traceability.

Every time ``collect_scores`` processes a docking output directory, a
**run_id** is generated and attached to every score row.  A companion
JSON sidecar (``<scores_file>.provenance.json``) maps each ``run_id``
to rich metadata about the run:

* engine name and version (when available)
* source directory
* collection timestamp (ISO-8601)
* hostname
* engine-specific parameters extracted from config files

This lets users (and future-you) answer: *"Where did this number come
from?"* — especially when the same pair is docked multiple times with
different settings.

The ``run_id`` column is added to the unified TSV but is **not** a
metric.  Downstream tools (``load_scores``, ``available_metrics``)
already skip non-metric columns.
"""

from __future__ import annotations

import datetime
import json
import os
import platform
import re
from typing import Any

# ---------------------------------------------------------------------------
# Run-ID generation
# ---------------------------------------------------------------------------

def make_run_id(engine: str,
                pair: tuple[str, str] | None = None,
                timestamp: datetime.datetime | None = None) -> str:
    """Build a short, human-readable run identifier.

    Format::

        <engine>_<pA>_<pB>_<YYYYMMDDTHHMMSS>

    or ``<engine>_<YYYYMMDDTHHMMSS>`` when *pair* is not supplied.

    Examples
    --------
    >>> make_run_id("lightdock", ("2UUY_rec", "2UUY_lig"))
    'lightdock_2UUY_rec_2UUY_lig_20260402T152300'
    """
    ts = timestamp or datetime.datetime.now()
    ts_str = ts.strftime("%Y%m%dT%H%M%S")
    parts = [engine]
    if pair:
        parts.extend(pair)
    parts.append(ts_str)
    return "_".join(parts)


# ---------------------------------------------------------------------------
# Engine-specific metadata extractors
# ---------------------------------------------------------------------------

def _extract_lightdock_meta(run_dir: str) -> dict[str, Any]:
    """Read ``setup.json`` from a LightDock simulation directory."""
    setup_path = os.path.join(run_dir, "setup.json")
    if not os.path.isfile(setup_path):
        return {"note": "setup.json not found"}
    with open(setup_path) as fh:
        return json.load(fh)


def _extract_haddock_meta(run_dir: str) -> dict[str, Any]:
    """Extract version, timestamp, and step list from a HADDOCK ``log``."""
    log_path = os.path.join(run_dir, "log")
    meta: dict[str, Any] = {}
    if not os.path.isfile(log_path):
        meta["note"] = "log file not found"
        return meta

    with open(log_path) as fh:
        for line in fh:
            # "Starting HADDOCK3 v2025.9.1 on 2025-11-17 12:54:00"
            _HADDOCK_START = (
                r"Starting HADDOCK3\s+(v[\d.]+)"
                r"\s+on\s+(\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2}:\d{2})"
            )
            m = re.search(_HADDOCK_START, line)
            if m:
                meta["haddock_version"] = m.group(1)
                meta["run_started"] = m.group(2)
            # "Python 3.9.25 (main, ...)"
            if line.strip().startswith("Python "):
                meta["python_version"] = line.strip()
            # Step list: "[...] Reading instructions step 0_topoaa"
            m2 = re.search(r"Reading instructions step (\S+)", line)
            if m2:
                meta.setdefault("steps", []).append(m2.group(1))
    return meta


def _extract_rosetta_meta(run_dir: str) -> dict[str, Any]:
    """Read Rosetta flag files from *run_dir* (or its parent)."""
    import glob as _glob

    meta: dict[str, Any] = {}
    # Rosetta flag files are typically named flag_*
    flag_files = sorted(
        _glob.glob(os.path.join(run_dir, "flag_*"))
        + _glob.glob(os.path.join(run_dir, "..", "flag_*"))
    )
    for fp in flag_files:
        name = os.path.basename(fp)
        with open(fp) as fh:
            meta[name] = fh.read().strip()
    if not flag_files:
        meta["note"] = "no flag_* files found"
    return meta


_ENGINE_EXTRACTORS: dict[str, Any] = {
    "lightdock": _extract_lightdock_meta,
    "haddock": _extract_haddock_meta,
    "rosetta": _extract_rosetta_meta,
}


def extract_run_metadata(engine: str, run_dir: str) -> dict[str, Any]:
    """Build a provenance record for one docking run.

    Returns a dict suitable for JSON serialisation::

        {
            "engine": "lightdock",
            "source_dir": "/abs/path/to/simulation",
            "collected_at": "2026-04-02T15:23:00",
            "hostname": "okik-mbp",
            "engine_meta": { ... }   # engine-specific config/params
        }
    """
    extractor = _ENGINE_EXTRACTORS.get(engine)
    engine_meta = extractor(run_dir) if extractor else {}

    return {
        "engine": engine,
        "source_dir": os.path.abspath(run_dir),
        "collected_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "hostname": platform.node(),
        "engine_meta": engine_meta,
    }


# ---------------------------------------------------------------------------
# Sidecar I/O
# ---------------------------------------------------------------------------

def sidecar_path(scores_path: str) -> str:
    """Return the provenance sidecar path for a given scores file.

    ``scores.tsv`` → ``scores.tsv.provenance.json``
    """
    return scores_path + ".provenance.json"


def write_sidecar(scores_path: str,
                  provenance: dict[str, dict]) -> str:
    """Write (or merge into) the provenance sidecar.

    *provenance* maps ``run_id → metadata_dict``.  If the sidecar
    already exists, new entries are merged (existing entries are
    preserved — they are never overwritten).

    Returns the sidecar file path.
    """
    out = sidecar_path(scores_path)
    existing: dict = {}
    if os.path.isfile(out):
        with open(out) as fh:
            try:
                existing = json.load(fh)
            except json.JSONDecodeError:
                existing = {}

    # Merge — new entries added, old entries untouched
    existing.update(provenance)

    with open(out, "w") as fh:
        json.dump(existing, fh, indent=2, default=str)
    return out


def read_sidecar(scores_path: str) -> dict[str, dict]:
    """Read the provenance sidecar, returning an empty dict if absent."""
    out = sidecar_path(scores_path)
    if not os.path.isfile(out):
        return {}
    with open(out) as fh:
        return json.load(fh)
