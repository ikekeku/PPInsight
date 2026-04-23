"""
collect_scores – aggregate docking results into a unified scores file.

Usage::

    collect_scores examples/haddock3/run1-test \
        examples/lightdock/simulation -o scores.tsv
    collect_scores path/to/rosetta_run \
        --model-label rosetta -o scores.tsv

Each positional argument is a directory produced by a supported
docking pipeline.  The tool auto-detects the docking engine and
extracts relevant scores, then writes a unified TSV/CSV that
:func:`ppinsight.visualizer.load_scores` can consume.

Unified output columns
----------------------
Required: ``model``, ``score_type``, ``score_value``, ``proteinA``,
``proteinB``

Optional traceability columns added when available:
``run_id``, ``pose_id``, ``output_path``, ``source_file``,
``pose_rank``
"""

import argparse
import datetime
import glob
import os
import re
import sys

import pandas as pd

from ppinsight.provenance import (
    extract_run_metadata,
    make_run_id,
    write_sidecar,
)


def _resolve_optional_path(base_dir: str, raw_path) -> str:
    """Resolve *raw_path* relative to *base_dir* when possible.

    Returns an empty string when the path is missing or placeholder-like.
    """
    if raw_path is None:
        return ""
    text = str(raw_path).strip()
    if not text or text in {"-", "nan", "None"}:
        return ""
    if os.path.isabs(text):
        return os.path.abspath(os.path.normpath(text))
    return os.path.abspath(os.path.normpath(os.path.join(base_dir, text)))


def _pose_id_from_ref(ref, fallback: str) -> str:
    """Build a stable pose identifier from a path-like reference."""
    if ref is None:
        return fallback
    text = str(ref).strip()
    if not text or text in {"-", "nan", "None"}:
        return fallback
    stem = os.path.splitext(os.path.basename(text))[0]
    return stem or fallback


def _rosetta_output_path(base_dir: str, description: str) -> str:
    """Find the most likely Rosetta decoy PDB for *description*."""
    candidates = [
        os.path.join(base_dir, f"docked_{description}.pdb"),
        os.path.join(base_dir, f"{description}.pdb"),
        os.path.join(base_dir, f"decoy_{description}.pdb"),
    ]
    for candidate in candidates:
        if os.path.isfile(candidate):
            return os.path.abspath(candidate)
    return ""

# ---------------------------------------------------------------------------
# HADDOCK parser
# ---------------------------------------------------------------------------

def _parse_haddock(run_dir: str,
                   pair: tuple[str, str] | None = None,
                   label: str = "haddock",
                   *,
                   no_clusters: bool = False) -> pd.DataFrame:
    """Parse HADDOCK results from a run directory.

    By default (``no_clusters=False``), this parser follows the HADDOCK best
    practice and looks for **cluster-level** results first:

    1. ``clustfcc`` or ``clustrmsd`` step outputs (``cluster.tsv``)
    2. Falls back to the per-model ``capri_ss.tsv`` from the last
       ``caprieval`` step.

    Pass ``no_clusters=True`` to skip cluster files and always use
    per-model ``capri_ss.tsv`` (useful for debugging or when clustering
    was intentionally disabled in the HADDOCK config).

    Returns a DataFrame with unified columns.
    """
    pA, pB = pair or ("", "")

    # ------------------------------------------------------------------
    # Try cluster-level results first (standard HADDOCK best practice)
    # ------------------------------------------------------------------
    if not no_clusters:
        cluster_rows = _parse_haddock_clusters(run_dir, pair=pair, label=label)
        if cluster_rows is not None:
            return cluster_rows

    # ------------------------------------------------------------------
    # Fall back to per-model capri_ss.tsv
    # ------------------------------------------------------------------
    pattern = os.path.join(run_dir, "**", "capri_ss.tsv")
    hits = sorted(glob.glob(pattern, recursive=True))
    if not hits:
        raise FileNotFoundError(
            f"No capri_ss.tsv found in {run_dir}. Is this a HADDOCK run directory?"
        )
    # Use the last caprieval stage (highest numbered)
    capri_path = os.path.abspath(hits[-1])
    df = pd.read_csv(capri_path, sep="\t")
    df.columns = [c.strip().lower() for c in df.columns]

    # Metrics we want to lift into the unified file
    metric_cols = [
        c for c in ("score", "dockq", "irmsd", "fnat", "lrmsd")
        if c in df.columns
    ]
    if not metric_cols:
        raise ValueError(
            f"capri_ss.tsv at {capri_path} has no recognized "
            "score columns"
        )

    rows: list[dict] = []
    capri_dir = os.path.dirname(capri_path)
    for row_idx, (_, row) in enumerate(df.iterrows(), start=1):
        pose_path = _resolve_optional_path(capri_dir, row.get("model"))
        pose_id = _pose_id_from_ref(row.get("model"), f"{label}_pose_{row_idx:04d}")
        pose_rank = int(row.get("caprieval_rank", row_idx))
        for metric in metric_cols:
            rows.append({
                "model": label,
                "score_type": metric,
                "score_value": float(row[metric]),
                "proteinA": pA,
                "proteinB": pB,
                "pose_id": pose_id,
                "output_path": pose_path,
                "source_file": capri_path,
                "pose_rank": pose_rank,
            })
    return pd.DataFrame(rows)


def _parse_haddock_clusters(run_dir: str,
                            pair: tuple[str, str] | None = None,
                            label: str = "haddock") -> pd.DataFrame | None:
    """Parse HADDOCK cluster-level results (``clustfcc`` / ``clustrmsd``).

    HADDOCK's ``clustfcc`` or ``clustrmsd`` modules write cluster summary
    files.  After clustering, the ``caprieval`` step that follows reports
    per-cluster statistics in its ``capri_ss.tsv`` — the ``model`` column
    will reference cluster centres instead of individual poses.

    If no cluster-aware caprieval output is found, returns *None* so the
    caller can fall back.
    """
    pA, pB = pair or ("", "")

    # Heuristic: if there is a clustfcc or clustrmsd directory, the
    # caprieval that *follows* it will contain cluster-level scores.
    cluster_step_dirs = sorted(
        glob.glob(os.path.join(run_dir, "*_clustfcc"))
        + glob.glob(os.path.join(run_dir, "*_clustrmsd"))
    )
    if not cluster_step_dirs:
        return None

    # Find the caprieval that comes *after* the clustering step.
    # HADDOCK steps are numbered: 7_clustfcc → 8_seletopclusts → 9_caprieval
    # We want the highest-numbered caprieval that sits after a clust step.
    all_capri = sorted(glob.glob(os.path.join(run_dir, "*_caprieval", "capri_ss.tsv")))
    if not all_capri:
        return None

    # Parse the last caprieval (post-clustering)
    capri_path = os.path.abspath(all_capri[-1])
    df = pd.read_csv(capri_path, sep="\t")
    df.columns = [c.strip().lower() for c in df.columns]

    metric_cols = [
        c for c in ("score", "dockq", "irmsd", "fnat", "lrmsd")
        if c in df.columns
    ]
    if not metric_cols:
        return None

    rows: list[dict] = []
    capri_dir = os.path.dirname(capri_path)
    for row_idx, (_, row) in enumerate(df.iterrows(), start=1):
        pose_path = _resolve_optional_path(capri_dir, row.get("model"))
        pose_id = _pose_id_from_ref(
            row.get("model"),
            f"{label}_cluster_pose_{row_idx:04d}",
        )
        pose_rank = int(row.get("caprieval_rank", row_idx))
        for metric in metric_cols:
            rows.append({
                "model": label,
                "score_type": metric,
                "score_value": float(row[metric]),
                "proteinA": pA,
                "proteinB": pB,
                "pose_id": pose_id,
                "output_path": pose_path,
                "source_file": capri_path,
                "pose_rank": pose_rank,
            })
    result = pd.DataFrame(rows)
    # Tag rows so downstream code knows these are cluster-level results
    result["source"] = "cluster"
    return result


# ---------------------------------------------------------------------------
# LightDock parser
# ---------------------------------------------------------------------------

# Regex used by the original parser to extract the last token from each
# gso_*.out line — currently unused but kept for reference.
_LIGHTDOCK_SCORE_RE = re.compile(r"\S+\s*$")  # last whitespace-separated token


def _parse_lightdock(sim_dir: str,
                     pair: tuple[str, str] | None = None,
                     label: str = "lightdock") -> pd.DataFrame:
    """Parse LightDock ``gso_<step>.out`` files from *sim_dir*.

    Uses the highest-step ``gso_*.out`` file in each ``swarm_*`` folder and
    takes the *Scoring* column (last column).

    Returns a DataFrame with unified columns.
    """
    swarm_dirs = sorted(glob.glob(os.path.join(sim_dir, "swarm_*")))
    if not swarm_dirs:
        raise FileNotFoundError(
            f"No swarm_* directories in {sim_dir}. Is this a LightDock simulation?"
        )

    rows: list[dict] = []
    pA, pB = pair or ("", "")

    for swarm in swarm_dirs:
        gso_files = sorted(glob.glob(os.path.join(swarm, "gso_*.out")))
        if not gso_files:
            continue
        # highest step number — gso_100.out beats gso_10.out
        def _step_num(p):
            m = re.search(r"gso_(\d+)\.out$", p)
            return int(m.group(1)) if m else -1

        last_gso = os.path.abspath(max(gso_files, key=_step_num))
        swarm_name = os.path.basename(swarm)
        pose_idx = 0
        with open(last_gso) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                # last whitespace-separated token is the scoring value
                parts = line.split()
                try:
                    score = float(parts[-1])
                except (ValueError, IndexError):
                    continue
                pose_path = os.path.join(swarm, f"lightdock_{pose_idx}.pdb")
                output_path = (
                    os.path.abspath(pose_path)
                    if os.path.isfile(pose_path)
                    else ""
                )
                rows.append({
                    "model": label,
                    "score_type": "luciferin_score",
                    "score_value": score,
                    "proteinA": pA,
                    "proteinB": pB,
                    "pose_id": f"{swarm_name}:lightdock_{pose_idx}",
                    "output_path": output_path,
                    "source_file": last_gso,
                    "pose_rank": pose_idx,
                    "swarm": swarm_name,
                })
                pose_idx += 1

    if not rows:
        raise ValueError(f"Could not extract any scores from {sim_dir}")
    return pd.DataFrame(rows)


def _parse_lightdock_clusters(sim_dir: str,
                              pair: tuple[str, str] | None = None,
                              label: str = "lightdock") -> pd.DataFrame:
    """Parse LightDock **cluster representatives** from *sim_dir*.

    After running ``lgd_cluster_bsas.py`` in each swarm directory, a file
    called ``cluster.repr`` is created.  Each line has the format::

        cluster_id : population : best_scoring : glowworm_id : pdb_file

    This parser reads those files and returns **one row per cluster
    representative** — the best-scoring pose from each structural cluster.

    Why clustering matters
    ~~~~~~~~~~~~~~~~~~~~~~
    LightDock generates hundreds of poses per swarm.  Many of these are
    near-duplicates (structurally very similar).  The BSAS clustering
    groups them by RMSD (default cutoff 4 Å) so that each cluster is a
    distinct binding mode.  The **representative** of each cluster is the
    pose with the highest score.  Using cluster representatives instead of
    all raw poses avoids counting the same binding mode many times and
    gives a much more meaningful "best score per distinct mode".

    If no ``cluster.repr`` files are found, falls back to the standard
    per-glowworm parser (:func:`_parse_lightdock`).

    Returns a DataFrame with unified columns, with an extra
    ``cluster_id`` and ``cluster_pop`` (population) column.
    """
    swarm_dirs = sorted(glob.glob(os.path.join(sim_dir, "swarm_*")))
    if not swarm_dirs:
        raise FileNotFoundError(
            f"No swarm_* directories in {sim_dir}. Is this a LightDock simulation?"
        )

    rows: list[dict] = []
    pA, pB = pair or ("", "")
    found_any = False

    for swarm in swarm_dirs:
        repr_file = os.path.join(swarm, "cluster.repr")
        if not os.path.isfile(repr_file):
            continue
        found_any = True
        swarm_name = os.path.basename(swarm)
        with open(repr_file) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split(":")
                if len(parts) < 5:
                    continue
                try:
                    cluster_id = int(parts[0])
                    population = int(parts[1])
                    score = float(parts[2])
                except (ValueError, IndexError):
                    continue
                pose_ref = parts[4].strip() if len(parts) >= 5 else ""
                output_path = _resolve_optional_path(swarm, pose_ref)
                pose_id = (
                    f"{swarm_name}:"
                    f"{_pose_id_from_ref(pose_ref, f'cluster_{cluster_id}')}"
                )
                rows.append({
                    "model": label,
                    "score_type": "luciferin_score",
                    "score_value": score,
                    "proteinA": pA,
                    "proteinB": pB,
                    "pose_id": pose_id,
                    "output_path": output_path,
                    "source_file": os.path.abspath(repr_file),
                    "pose_rank": cluster_id,
                    "swarm": swarm_name,
                    "cluster_id": cluster_id,
                    "cluster_pop": population,
                })

    if not found_any:
        # No cluster files → fall back to the raw glowworm parser
        return _parse_lightdock(sim_dir, pair=pair, label=label)

    if not rows:
        raise ValueError(f"All cluster.repr files in {sim_dir} were empty")
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Score aggregation
# ---------------------------------------------------------------------------

def aggregate_scores(
    scores_df: pd.DataFrame,
    strategy: str = "best",
    n: int = 5,
) -> pd.DataFrame:
    """Reduce many poses per (model, pair) down to a single representative score.

    After docking, each engine produces many individual pose scores for one
    protein pair — sometimes hundreds.  For pair-level comparison (e.g. "does
    pair A-B score better than pair C-D?") you need a single number.

    **Direction-awareness**: ``"best"`` and ``"topN_mean"`` respect score
    direction.  For higher-is-better metrics (e.g. LightDock luciferin)
    the highest scores are selected; for lower-is-better metrics (e.g.
    HADDOCK score, Rosetta I_sc) the lowest scores are selected.  The
    direction is looked up from :data:`ppinsight.visualizer.METRIC_METADATA`.

    Strategies
    ----------
    * **best** — keep only the single best score per group.
    * **topN_mean** — average the top *n* scores.
    * **median** — take the median of all scores.
    * **mean** — take the mean of all scores.

    Parameters
    ----------
    scores_df : DataFrame
        Unified scores.
    strategy : str
        One of ``"best"``, ``"topN_mean"``, ``"median"``, ``"mean"``.
    n : int
        How many top scores to average when strategy is ``"topN_mean"``.

    Returns
    -------
    DataFrame
        One row per (model, score_type, proteinA, proteinB) combination.
    """
    # Lazy import to avoid circular dependency (visualizer imports us too).
    from ppinsight.visualizer import get_metric_direction

    STRATEGIES = ("best", "topN_mean", "median", "mean")
    if strategy not in STRATEGIES:
        raise ValueError(f"strategy must be one of {STRATEGIES}, got '{strategy}'")

    df = scores_df.copy()
    df["score_value"] = pd.to_numeric(df["score_value"], errors="coerce")

    group_cols = ["model", "score_type"]
    if "run_id" in df.columns:
        group_cols.append("run_id")
    if "proteinA" in df.columns:
        group_cols.append("proteinA")
    if "proteinB" in df.columns:
        group_cols.append("proteinB")
    # Carry labels through if present
    extra_cols = [c for c in ("label", "family") if c in df.columns]
    group_cols_full = group_cols + extra_cols

    def _agg(g):
        # Determine sort direction from the score_type of this group.
        # When include_groups=False, the group key columns are stripped
        # from `g`, but they're available in `g.name` (a tuple matching
        # group_cols_full).  score_type is always at index 1.
        group_key = g.name if isinstance(g.name, tuple) else (g.name,)
        st_idx = (
            group_cols_full.index("score_type")
            if "score_type" in group_cols_full
            else -1
        )
        st = str(group_key[st_idx]) if st_idx >= 0 and st_idx < len(group_key) else ""
        # For lower-is-better metrics we sort ascending so .iloc[0] / .head()
        # give the best (lowest) values; for higher-is-better we sort descending.
        higher_better = get_metric_direction(st)
        vals = g["score_value"].dropna().sort_values(ascending=not higher_better)
        if strategy == "best":
            return vals.iloc[0] if len(vals) else float("nan")
        elif strategy == "topN_mean":
            return vals.head(n).mean()
        elif strategy == "median":
            return vals.median()
        else:  # mean
            return vals.mean()

    grouped = df.groupby(group_cols_full, as_index=False)
    try:
        result = grouped.apply(
            lambda g: pd.Series({"score_value": _agg(g)}),
            include_groups=False,
        )
    except TypeError:
        # pandas < 2.2 does not support include_groups
        result = grouped.apply(lambda g: pd.Series({"score_value": _agg(g)}))
    # Flatten if needed
    if isinstance(result.columns, pd.MultiIndex):
        result.columns = ["_".join(c).strip("_") for c in result.columns]
    return result


# ---------------------------------------------------------------------------
# Rosetta parser
# ---------------------------------------------------------------------------

def _parse_rosetta(out_dir: str,
                   pair: tuple[str, str] | None = None,
                   label: str = "rosetta") -> pd.DataFrame:
    """Parse Rosetta docking results from *out_dir*.

    Supports **three** output formats (checked in priority order):

    1. **Clustered scores CSV** (``clustered_scores.csv``) — written by
       :func:`ppinsight.rosetta.analyze.cluster_and_rank` after hierarchical
       Cα-RMSD clustering.  Only **cluster representatives** (the best-
       scoring member of each cluster) are returned.

    2. **Native Rosetta ``.sc`` score files** — space-separated tables
       produced by ``-out:file:scorefile``.  Key metrics extracted:

       - ``I_sc`` (interface score — primary RosettaDock quality metric)
       - ``total_score``
       - ``Irms`` (interface RMSD, if present)
       - ``rms`` (ligand RMSD, if present)
       - ``Fnat`` (fraction of native contacts, if present)

    3. **PPInsight ``docking_scores.csv``** — the legacy CSV written by
       the PyRosetta wrapper (``ppinsight.rosetta.dock``).  Falls back to
       this when no ``.sc`` file is found.

    Returns a DataFrame with unified columns.
    """
    # ----- try clustered_scores.csv first (best practice) -----
    clustered_path = os.path.join(out_dir, "clustered_scores.csv")
    if os.path.isfile(clustered_path):
        return _parse_rosetta_clustered(clustered_path, pair=pair, label=label)

    # ----- try .sc files (native Rosetta) -----
    sc_files = sorted(glob.glob(os.path.join(out_dir, "*.sc")))
    if sc_files:
        return _parse_rosetta_scorefile(sc_files, pair=pair, label=label)

    # ----- fall back to legacy PPInsight CSV -----
    csv_path = os.path.abspath(os.path.join(out_dir, "docking_scores.csv"))
    if not os.path.isfile(csv_path):
        raise FileNotFoundError(
            f"No clustered_scores.csv, .sc score file, or docking_scores.csv "
            f"in {out_dir}. Is this a Rosetta output directory?"
        )

    df = pd.read_csv(csv_path)
    df.columns = [c.strip().lower() for c in df.columns]

    rows: list[dict] = []
    pA, pB = pair or ("", "")
    score_col = "score" if "score" in df.columns else df.columns[-1]
    for row_idx, (_, row) in enumerate(df.iterrows(), start=1):
        run_ref = row.get("run", row_idx)
        pose_id = f"run_{int(run_ref):04d}" if str(run_ref).isdigit() else str(run_ref)
        rows.append({
            "model": label,
            "score_type": "interface_score",
            "score_value": float(row[score_col]),
            "proteinA": pA,
            "proteinB": pB,
            "pose_id": pose_id,
            "output_path": "",
            "source_file": csv_path,
            "pose_rank": row_idx,
        })
    return pd.DataFrame(rows)


def _parse_rosetta_clustered(clustered_csv: str,
                             pair: tuple[str, str] | None = None,
                             label: str = "rosetta") -> pd.DataFrame:
    """Parse a ``clustered_scores.csv`` produced by decoy clustering.

    Only **cluster representatives** (``is_top_of_cluster == True``) are
    returned — one row per distinct structural cluster.  This follows
    the RosettaDock best practice of reporting results per binding mode
    rather than per individual decoy.

    Recognised metric columns (case-insensitive):
    ``i_sc``, ``total_score``, ``irms``, ``rms``, ``fnat``,
    ``dg_separated``, ``cluster_size``.
    """
    _COL_MAP = {
        "i_sc": "i_sc",
        "total_score": "total_score",
        "irms": "irms",
        "rms": "rms",
        "fnat": "fnat",
        "dg_separated": "dg_separated",
        "cluster_size": "cluster_size",
    }

    clustered_csv = os.path.abspath(clustered_csv)
    df = pd.read_csv(clustered_csv)
    df.columns = [c.strip().lower() for c in df.columns]

    # Only keep cluster representatives
    if "is_top_of_cluster" in df.columns:
        df = df[df["is_top_of_cluster"].astype(bool)]

    pA, pB = pair or ("", "")
    rows: list[dict] = []
    rosetta_dir = os.path.dirname(clustered_csv)

    for row_idx, (_, row) in enumerate(df.iterrows(), start=1):
        description = str(row.get("description", "")).strip()
        pose_id = description or f"rosetta_pose_{row_idx:04d}"
        output_path = (
            _rosetta_output_path(rosetta_dir, description)
            if description else ""
        )
        for csv_col, our_name in _COL_MAP.items():
            if csv_col in df.columns:
                try:
                    v = float(row[csv_col])
                except (ValueError, TypeError):
                    continue
                rows.append({
                    "model": label,
                    "score_type": our_name,
                    "score_value": v,
                    "proteinA": pA,
                    "proteinB": pB,
                    "pose_id": pose_id,
                    "output_path": output_path,
                    "source_file": clustered_csv,
                    "pose_rank": row_idx,
                    "cluster_id": (
                        int(row["cluster"])
                        if "cluster" in df.columns else None
                    ),
                    "cluster_rank": (
                        int(row["cluster_rank"])
                        if "cluster_rank" in df.columns else None
                    ),
                    "cluster_size": (
                        int(row["cluster_size"])
                        if "cluster_size" in df.columns else None
                    ),
                })

    if not rows:
        raise ValueError(
            f"No recognized score columns in {clustered_csv}. "
            f"Columns found: {list(df.columns)}"
        )
    result = pd.DataFrame(rows)
    result["source"] = "cluster"
    return result


def _parse_rosetta_scorefile(sc_files: list[str],
                             pair: tuple[str, str] | None = None,
                             label: str = "rosetta") -> pd.DataFrame:
    """Parse one or more Rosetta ``.sc`` score files.

    Rosetta ``.sc`` files are space-separated, header on the second line
    (first line is ``SEQUENCE:``).  We extract every recognized metric.

    Recognized columns (case-insensitive):
    ``I_sc``, ``total_score``, ``Irms``, ``rms``, ``Fnat``
    """
    # Mapping from .sc column name (lowered) to our unified score_type
    _COL_MAP = {
        "i_sc": "i_sc",
        "total_score": "total_score",
        "irms": "irms",
        "rms": "rms",
        "fnat": "fnat",
    }

    all_rows: list[dict] = []
    pA, pB = pair or ("", "")

    for sc_path in sc_files:
        sc_path = os.path.abspath(sc_path)
        header = None
        row_idx = 0
        with open(sc_path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                # First non-blank line starting with SEQUENCE: is a preamble
                if line.upper().startswith("SEQUENCE:"):
                    continue
                parts = line.split()
                # Header line starts with SCORE:
                if parts[0].upper() == "SCORE:":
                    if header is None:
                        header = [c.lower() for c in parts[1:]]
                        continue
                    # Data row
                    vals = parts[1:]
                    if header and len(vals) == len(header):
                        row_idx += 1
                        row_dict = dict(zip(header, vals, strict=False))
                        description = str(row_dict.get("description", "")).strip()
                        pose_id = description or f"rosetta_pose_{row_idx:04d}"
                        output_path = _rosetta_output_path(
                            os.path.dirname(sc_path), description,
                        )
                        for sc_col, our_name in _COL_MAP.items():
                            if sc_col in row_dict:
                                try:
                                    v = float(row_dict[sc_col])
                                except ValueError:
                                    continue
                                all_rows.append({
                                    "model": label,
                                    "score_type": our_name,
                                    "score_value": v,
                                    "proteinA": pA,
                                    "proteinB": pB,
                                    "pose_id": pose_id,
                                    "output_path": output_path,
                                    "source_file": sc_path,
                                    "pose_rank": row_idx,
                                })

    if not all_rows:
        raise ValueError(
            f"No recognized scores in .sc files: {sc_files}"
        )
    return pd.DataFrame(all_rows)


# ---------------------------------------------------------------------------
# Auto-detection
# ---------------------------------------------------------------------------

def detect_engine(directory: str) -> str:
    """Guess which docking engine produced *directory*.

    Delegates to :func:`ppinsight.registry.detect_engine` so that any
    registered plugin's detector is checked automatically.

    Returns one of the registered engine names (e.g. ``'haddock'``,
    ``'lightdock'``, ``'rosetta'``).

    Raises ``ValueError`` if the directory doesn't match any known layout.
    """
    from ppinsight.registry import detect_engine as _detect
    return _detect(directory)


def _get_parser(engine: str):
    """Return the parser callable for *engine* from the registry."""
    from ppinsight.registry import get as _get_plugin
    plugin = _get_plugin(engine)
    if plugin.parser is None:
        raise ValueError(f"Engine '{engine}' has no parser registered.")
    return plugin.parser



def collect(directories: list[str],
            labels: list[str] | None = None,
            pair: tuple[str, str] | None = None,
            use_clusters: bool = False,
            no_haddock_clusters: bool = False) -> tuple[pd.DataFrame, dict]:
    """Parse multiple docking output directories into a unified DataFrame.

    Parameters
    ----------
    directories : list[str]
        Paths to docking output directories.
    labels : list[str] | None
        Optional per-directory model label overrides.  If *None*, the
        detected engine name is used (``'haddock'``, ``'lightdock'``,
        ``'rosetta'``).
    pair : tuple[str, str] | None
        ``(proteinA, proteinB)`` to fill in every row.
    use_clusters : bool
        For LightDock directories, use ``cluster.repr`` files instead of
        raw ``gso_*.out``.  Falls back to raw if no cluster files exist.
    no_haddock_clusters : bool
        If *True*, skip automatic HADDOCK cluster parsing and always use
        per-model ``capri_ss.tsv``.  By default (False), HADDOCK cluster
        results are parsed automatically (HADDOCK best practice).

    Returns
    -------
    (pd.DataFrame, dict)
        the unified scores DataFrame and a provenance dict mapping each
        ``run_id`` to its metadata.  The DataFrame retains lightweight
        row-level traceability columns such as ``run_id``, ``pose_id``,
        ``output_path``, and ``source_file`` when the engine output
        exposes them.
    """
    if labels and len(labels) != len(directories):
        raise ValueError(
            f"Got {len(labels)} labels but {len(directories)} directories"
        )

    now = datetime.datetime.now()
    frames: list[pd.DataFrame] = []
    provenance: dict[str, dict] = {}

    for i, d in enumerate(directories):
        engine = detect_engine(d)
        label = labels[i] if labels else engine

        # --- parse scores as before ---
        if use_clusters and engine == "lightdock":
            frame = _parse_lightdock_clusters(d, pair=pair, label=label)
        elif engine == "haddock":
            frame = _parse_haddock(d, pair=pair, label=label,
                                   no_clusters=no_haddock_clusters)
        else:
            parser = _get_parser(engine)
            frame = parser(d, pair=pair, label=label)

        # --- provenance: one run_id per (directory, engine) ---
        rid = make_run_id(
            engine,
            pair=pair,
            timestamp=now + datetime.timedelta(seconds=i),
        )
        provenance[rid] = extract_run_metadata(engine, d)

        frame = frame.copy()
        frame["run_id"] = rid

        frames.append(frame)

    if not frames:
        raise ValueError("No directories processed")
    return pd.concat(frames, ignore_index=True), provenance


def annotate_with_labels(scores_df: pd.DataFrame,
                         pairs_df: pd.DataFrame) -> pd.DataFrame:
    """Add ``label`` and ``family`` columns to *scores_df* by matching pairs.

    The merge is tried in both orderings of (proteinA, proteinB) because
    different engines may list the proteins in either order.  Case-
    insensitive matching is used because protein names often differ in
    capitalisation across data sources.

    *pairs_df* must have columns ``proteinA``, ``proteinB``, ``label`` and
    optionally ``family``.  Matching is case-insensitive on protein names.

    Returns a copy of *scores_df* with the extra columns (unmatched rows
    get ``label='unknown'``).
    """
    df = scores_df.copy()
    pf = pairs_df[["proteinA", "proteinB", "label"]].copy()
    if "family" in pairs_df.columns:
        pf["family"] = pairs_df["family"]

    # Normalise for matching (case-insensitive, whitespace-stripped)
    df["_pA"] = df["proteinA"].str.strip().str.upper()
    df["_pB"] = df["proteinB"].str.strip().str.upper()
    pf["_pA"] = pf["proteinA"].str.strip().str.upper()
    pf["_pB"] = pf["proteinB"].str.strip().str.upper()

    # Merge on (proteinA, proteinB) — first try direct ordering
    merged = df.merge(
        pf[["_pA", "_pB", "label"] + (["family"] if "family" in pf.columns else [])],
        on=["_pA", "_pB"], how="left", suffixes=("", "_anno"),
    )
    # Try reversed pair for rows that didn't match the first time.
    # This handles the case where pairs_df has (A, B) but scores_df has (B, A).
    unmatched = (
        merged["label_anno"].isna()
        if "label_anno" in merged.columns
        else merged["label"].isna()
    )
    if unmatched.any():
        rev = pf.rename(columns={"_pA": "_pB", "_pB": "_pA"})
        merged2 = df[unmatched].merge(
            rev[
                ["_pA", "_pB", "label"]
                + (["family"] if "family" in rev.columns else [])
            ],
            on=["_pA", "_pB"], how="left", suffixes=("", "_rev"),
        )
        label_col = "label_anno" if "label_anno" in merged.columns else "label"
        if "label_rev" in merged2.columns:
            merged.loc[unmatched, label_col] = merged2["label_rev"].values
        if "family" in pf.columns:
            fam_col = "family_anno" if "family_anno" in merged.columns else "family"
            if "family_rev" in merged2.columns:
                merged.loc[unmatched, fam_col] = merged2["family_rev"].values

    # Clean up
    label_src = "label_anno" if "label_anno" in merged.columns else "label"
    df["label"] = merged[label_src].fillna("unknown")
    if "family" in pf.columns:
        fam_src = "family_anno" if "family_anno" in merged.columns else "family"
        if fam_src in merged.columns:
            df["family"] = merged[fam_src].fillna("")
    df = df.drop(columns=["_pA", "_pB"], errors="ignore")
    return df


def _has_nonempty_pair_context(scores_df: pd.DataFrame) -> bool:
    """Return True when collected rows include usable protein pair names."""
    required = {"proteinA", "proteinB"}
    if not required.issubset(scores_df.columns):
        return False

    pa = scores_df["proteinA"].fillna("").astype(str).str.strip()
    pb = scores_df["proteinB"].fillna("").astype(str).str.strip()
    return bool((pa.ne("") & pb.ne("")).any())


def _slugify_filename_part(text: str) -> str:
    """Return a filesystem-safe token for default output filenames."""
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", text.strip())
    cleaned = cleaned.strip("._")
    return cleaned or "run"


def _next_available_file_path(path: str) -> str:
    """Return *path* if free, else append ``_N`` before the extension."""
    if not os.path.exists(path):
        return path

    base, ext = os.path.splitext(path)
    idx = 1
    while True:
        candidate = f"{base}_{idx}{ext}"
        if not os.path.exists(candidate):
            return candidate
        idx += 1


def _default_output_path(directories: list[str]) -> str:
    """Build a descriptive default output path from input directory names."""
    out_dir = os.path.join("data", "output", "scores")
    stems = [
        _slugify_filename_part(os.path.basename(os.path.normpath(d)))
        for d in directories
    ]
    if len(stems) == 1:
        stem = f"scores_{stems[0]}"
    else:
        stem = f"scores_{stems[0]}_plus_{len(stems) - 1}_runs"
    return _next_available_file_path(os.path.join(out_dir, f"{stem}.tsv"))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    """CLI entrypoint for collecting docking scores into a unified TSV/CSV."""
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate docking scores from one or more engine output "
            "directories into a single unified scores file."
        ),
    )
    parser.add_argument(
        "directories",
        nargs="+",
        help=(
            "Paths to docking output directories.  The producing engine is "
            "auto-detected from directory contents."
        ),
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help=(
            "Output file path.  If omitted, collect writes an auto-named "
            "TSV in data/output/scores/ based on input run directory names "
            "and appends a numeric suffix when needed to avoid overwriting.  "
            "Extension determines the format: .tsv → tab-separated, "
            ".csv → comma-separated.  Parent directories are created "
            "automatically.  This unified file is the input for "
            "'ppinsight compare'."
        ),
    )
    parser.add_argument(
        "--labels",
        nargs="*",
        default=None,
        help=(
            "Model display labels for each directory (default: auto-detected "
            "engine name, e.g. 'lightdock', 'haddock').  Must match the "
            "number of directories if given.  Use this when you have "
            "multiple runs of the same engine and need distinct names "
            "(e.g. 'haddock_ambig' vs 'haddock_noambig')."
        ),
    )
    parser.add_argument(
        "--pair",
        default=None,
        help=(
            "Protein pair as 'proteinA:proteinB'.  Fills the proteinA and "
            "proteinB columns in the output.  Required whenever collected "
            "rows do not already include protein names (common for single "
            "run directories).  Keep this set even when using --pairs so "
            "label annotation can match rows correctly."
        ),
    )
    parser.add_argument(
        "--pairs",
        default=None,
        help=(
            "Path to a pairs CSV/TSV (from 'ppinsight parse') to annotate "
            "each score row with an 'interaction' or 'non-interaction' label.  "
            "Enables --classify and --plot-type roc in 'ppinsight compare'.  "
            "This flag does not infer proteinA/proteinB values from "
            "directories; it only labels rows that already have pair names.  "
            "Omit when ground-truth labels are unavailable."
        ),
    )
    parser.add_argument(
        "--use-clusters",
        action="store_true",
        help=(
            "For LightDock directories, read cluster representative files "
            "(cluster.repr) instead of raw gso_*.out swarm outputs.  "
            "Cluster representatives give one score per distinct binding "
            "mode per swarm, which is more meaningful and less noisy.  "
            "Enable when LightDock post-processing (clustering) has been "
            "run; omit for raw-pose analysis or when clustering was skipped."
        ),
    )
    parser.add_argument(
        "--no-haddock-clusters",
        action="store_true",
        help=(
            "Disable automatic HADDOCK cluster-level parsing.  By default, "
            "cluster-level results from clustfcc/clustrmsd caprieval steps "
            "are preferred over per-model capri_ss.tsv — this is HADDOCK "
            "best practice because cluster-averaged scores are more robust.  "
            "Use this flag only when you want to analyse raw per-model "
            "scores (e.g. for debugging or per-structure DockQ evaluation)."
        ),
    )
    parser.add_argument(
        "--agg",
        choices=["best", "topN_mean", "median", "mean"],
        default=None,
        help=(
            "Reduce many poses to one score per (model, pair, score_type).  "
            "'best' keeps only the top-ranked score (direction-aware — the "
            "pipeline knows whether lower or higher is better for each "
            "metric).  'topN_mean' averages the top N (see --agg-n).  "
            "'median' and 'mean' use all poses.  Use 'best' or 'topN_mean' "
            "for production benchmarks; omit for full-distribution analysis."
        ),
    )
    parser.add_argument(
        "--agg-n",
        type=int,
        default=5,
        help=(
            "How many top scores to average for --agg topN_mean (default: 5).  "
            "Higher values smooth out stochastic noise but may dilute the "
            "signal if only one or two poses are near-native."
        ),
    )

    args = parser.parse_args(argv)

    pair: tuple[str, str] | None = None
    if args.pair:
        parts = args.pair.split(":")
        if len(parts) != 2:
            print("ERROR: --pair must be 'proteinA:proteinB'", file=sys.stderr)
            sys.exit(2)
        pair = (parts[0], parts[1])

    try:
        df, provenance = collect(
            args.directories, labels=args.labels, pair=pair,
            use_clusters=args.use_clusters,
            no_haddock_clusters=args.no_haddock_clusters,
        )
    except FileNotFoundError as exc:
        msg = str(exc)
        print(f"ERROR: {msg}", file=sys.stderr)
        # Suggest helpful flags based on the error message
        if "capri_ss.tsv" in msg:
            print("Hint: if this is a HADDOCK run without caprieval, try "
                  "--no-haddock-clusters to parse per-model scores instead.",
                  file=sys.stderr)
        elif "swarm_" in msg and args.use_clusters:
            print("Hint: --use-clusters was set but no swarm directories "
                  "were found.  Verify the directory is a LightDock "
                  "simulation output.", file=sys.stderr)
        sys.exit(1)
    except ValueError as exc:
        msg = str(exc)
        print(f"ERROR: {msg}", file=sys.stderr)
        if "labels" in msg.lower():
            print("Hint: --labels must have exactly one entry per directory, "
                  "or omit it to auto-detect engine names.", file=sys.stderr)
        sys.exit(1)

    # Annotate with interaction labels if a pairs file is given
    if args.pairs:
        if not _has_nonempty_pair_context(df):
            print(
                "ERROR: --pairs requires populated proteinA/proteinB columns "
                "in collected scores.",
                file=sys.stderr,
            )
            print(
                "Hint: for single-run collection, pass --pair "
                "proteinA:proteinB.  --pairs adds labels to existing pair "
                "names; it does not infer pair names from directories.",
                file=sys.stderr,
            )
            sys.exit(2)
        pairs_sep = "\t" if args.pairs.endswith(".tsv") else ","
        pairs_df = pd.read_csv(args.pairs, sep=pairs_sep)
        df = annotate_with_labels(df, pairs_df)

    # Aggregate if requested
    if args.agg:
        df = aggregate_scores(df, strategy=args.agg, n=args.agg_n)

    # Write output
    out = args.output or _default_output_path(args.directories)
    if args.output is None:
        print(f"No --output provided; auto-selected: {out}")
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    sep = "\t" if out.endswith(".tsv") else ","
    df.to_csv(out, sep=sep, index=False)
    print(f"Wrote {len(df)} score rows to {out}")

    # Write provenance sidecar
    if provenance:
        sc = write_sidecar(out, provenance)
        print(f"Wrote provenance to {sc}")

    # Always print a summary so the user can sanity-check score ranges.
    print("\n── Summary ──")
    summary = (
        df.groupby(["model", "score_type"])["score_value"]
        .agg(["count", "mean", "min", "max"])
        .reset_index()
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
