"""
quality – DockQ-based quality assessment for docking predictions.

Computes CAPRI-standard quality metrics (DockQ, Fnat, I-RMSD, L-RMSD)
by comparing predicted docked complexes against a known native structure.

This module wraps the `DockQ <https://github.com/wallnerlab/DockQ>`_
library and adds:

* CAPRI quality classification (incorrect / acceptable / medium / high)
* Batch evaluation of entire docking run directories
* Integration with the PPInsight unified scores format

Usage::

    # Single model evaluation
    from ppinsight.quality import evaluate_complex
    result = evaluate_complex("model.pdb", "native.pdb")
    print(result)
    # {'DockQ': 0.85, 'fnat': 0.9, 'iRMSD': 1.2, 'LRMSD': 2.1,
    #  'capri_class': 'high', ...}

    # Batch evaluation of a directory
    from ppinsight.quality import evaluate_directory
    df = evaluate_directory("docking_output/", "native.pdb")

CLI::

    ppinsight quality model.pdb native.pdb
    ppinsight quality docking_output/ native.pdb --engine lightdock
    ppinsight quality scores.tsv native.pdb -o scores_quality.tsv
"""

import argparse
import glob
import os
import sys
from typing import Any

import pandas as pd

from ppinsight.utils import find_column, resolve_traceability_path

# ---------------------------------------------------------------------------
# DockQ import (optional – installed via ``pip install ppinsight[quality]``)
#
# DockQ is a third-party tool for computing CAPRI-standard quality metrics.
# It is distributed separately and currently requires a fork for numpy ≥ 2
# compatibility (upstream PR #61).  We make the import optional so that:
#   1. ``import ppinsight`` never fails due to a missing DockQ install.
#   2. Pure-Python helpers (classify_capri, capri_summary, etc.) remain
#      usable without DockQ.
#   3. Functions that truly need DockQ (evaluate_complex, evaluate_directory)
#      accept injectable _load_fn/_run_fn overrides for testing, and only
#      check for DockQ when the real (non-injected) functions are used.
# ---------------------------------------------------------------------------

_DOCKQ_AVAILABLE = False
_DOCKQ_IMPORT_ERROR: str | None = None

try:
    from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
    _DOCKQ_AVAILABLE = True
except ImportError as _exc:
    _DOCKQ_IMPORT_ERROR = (
        f"DockQ is not installed or failed to import ({_exc}).  "
        "Install it with:  pip install 'ppinsight[quality]'  "
        "or:  pip install --no-binary=DockQ --no-cache-dir "
        "'dockq @ git+https://github.com/nrontsis/DockQ.git@update-to-numpy>2'"
    )

    # Stub functions so the module can still be imported.
    # Calling them raises ImportError with install instructions.
    def load_PDB(*a, **kw):  # noqa: N802
        raise ImportError(_DOCKQ_IMPORT_ERROR)

    def run_on_all_native_interfaces(*a, **kw):
        raise ImportError(_DOCKQ_IMPORT_ERROR)


def _require_dockq() -> None:
    """Raise a helpful ImportError if DockQ is not available.

    Called at the top of functions that need the real DockQ backend
    (i.e. when no injectable test doubles are provided).
    """
    if not _DOCKQ_AVAILABLE:
        raise ImportError(_DOCKQ_IMPORT_ERROR)


# ---------------------------------------------------------------------------
# CAPRI classification thresholds
# ---------------------------------------------------------------------------

#: CAPRI quality thresholds on the DockQ score (Basu & Wallner 2016).
#: A model is classified by the *highest* category whose threshold it meets.
CAPRI_THRESHOLDS: dict[str, float] = {
    "high": 0.80,
    "medium": 0.49,
    "acceptable": 0.23,
}

#: Ordered from best to worst for display / sorting.
CAPRI_ORDER: list[str] = ["high", "medium", "acceptable", "incorrect"]


def classify_capri(dockq: float) -> str:
    """Classify a DockQ score into a CAPRI quality category.

    Parameters
    ----------
    dockq : float
        DockQ score in [0, 1].

    Returns
    -------
    str
        One of ``'high'``, ``'medium'``, ``'acceptable'``, ``'incorrect'``.
    """
    for label, threshold in CAPRI_THRESHOLDS.items():
        if dockq >= threshold:
            return label
    return "incorrect"


# ---------------------------------------------------------------------------
# Core evaluation
# ---------------------------------------------------------------------------

def evaluate_complex(
    model_path: str,
    native_path: str,
    *,
    chain_map: dict[str, str] | None = None,
    capri_peptide: bool = False,
    _load_fn=None,
    _run_fn=None,
) -> dict[str, Any]:
    """Evaluate a docked model against a native structure using DockQ.

    Parameters
    ----------
    model_path : str
        Path to the predicted docked complex PDB/mmCIF.
    native_path : str
        Path to the native (reference) complex PDB/mmCIF.
    chain_map : dict[str, str] | None
        Mapping of model chain IDs to native chain IDs.  If *None*,
        DockQ will attempt automatic chain mapping.
    capri_peptide : bool
        Use CAPRI peptide thresholds (shorter interface distances).
    _load_fn, _run_fn
        Injectable overrides for testing (default: DockQ functions).

    Returns
    -------
    dict
        Keys include ``DockQ``, ``fnat``, ``iRMSD``, ``LRMSD``,
        ``capri_class``, ``F1``, ``clashes``, ``model_path``,
        ``native_path``, and per-interface details.
    """
    loader = _load_fn or load_PDB
    runner = _run_fn or run_on_all_native_interfaces

    # Only require the real DockQ backend when no test doubles are injected.
    # This lets tests pass _load_fn/_run_fn fakes without needing DockQ,
    # while real callers get a clear ImportError with install instructions.
    if _load_fn is None or _run_fn is None:
        _require_dockq()

    model = loader(model_path)
    native = loader(native_path)

    kwargs: dict[str, Any] = {}
    if chain_map is not None:
        kwargs["chain_map"] = chain_map
    kwargs["capri_peptide"] = capri_peptide

    # DockQ returns (result_dict, total_dockq_score).
    # result_dict is keyed by interface label (e.g. "AB") and each value
    # is a dict with DockQ, fnat, iRMSD, LRMSD, F1, clashes, etc.
    result_dict, total_dockq = runner(model, native, **kwargs)

    # Build a flat summary — most docking benchmarks have a single
    # interface (receptor–ligand), so we surface its metrics at the
    # top level for convenience.
    summary: dict[str, Any] = {
        "model_path": model_path,
        "native_path": native_path,
        "DockQ": float(total_dockq),
        "capri_class": classify_capri(float(total_dockq)),
    }

    # Extract per-interface metrics
    interfaces = {}
    for iface_key, metrics in result_dict.items():
        interfaces[iface_key] = dict(metrics)
        # For the overall summary, use the first interface's detailed metrics
        if "fnat" not in summary:
            for field in ("fnat", "iRMSD", "LRMSD", "F1", "clashes",
                          "nat_correct", "nat_total", "fnonnat",
                          "nonnat_count", "model_total"):
                if field in metrics:
                    summary[field] = float(metrics[field])

    summary["interfaces"] = interfaces
    summary["n_interfaces"] = len(interfaces)

    return summary


def evaluate_directory(
    run_dir: str,
    native_path: str,
    *,
    engine: str | None = None,
    chain_map: dict[str, str] | None = None,
    glob_pattern: str | None = None,
    capri_peptide: bool = False,
    _load_fn=None,
    _run_fn=None,
) -> pd.DataFrame:
    """Evaluate all docked models in a directory against a native structure.

    Parameters
    ----------
    run_dir : str
        Directory containing docked model PDB files.
    native_path : str
        Path to the native (reference) complex PDB/mmCIF.
    engine : str | None
        Docking engine name — used to determine file layout.
        If *None*, tries auto-detection.
    chain_map : dict[str, str] | None
        Chain mapping (model → native).  Passed to :func:`evaluate_complex`.
    glob_pattern : str | None
        Custom glob pattern for finding model PDB files.  If *None*, an
        engine-appropriate pattern is used.
    capri_peptide : bool
        Use CAPRI peptide thresholds.
    _load_fn, _run_fn
        Injectable overrides for testing.

    Returns
    -------
    pd.DataFrame
        One row per model with columns ``model_path``, ``DockQ``,
        ``fnat``, ``iRMSD``, ``LRMSD``, ``capri_class``, …
    """
    model_files = _find_model_files(run_dir, engine=engine,
                                     glob_pattern=glob_pattern)
    if not model_files:
        raise FileNotFoundError(
            f"No model PDB files found in {run_dir}"
            + (f" for engine={engine}" if engine else "")
        )

    rows: list[dict[str, Any]] = []
    for mpath in sorted(model_files):
        try:
            result = evaluate_complex(
                mpath, native_path,
                chain_map=chain_map,
                capri_peptide=capri_peptide,
                _load_fn=_load_fn,
                _run_fn=_run_fn,
            )
            # Flatten: drop the nested interfaces dict for the DataFrame
            # (each interface's details are still in the full result dict
            # returned by evaluate_complex, but tabular output is simpler).
            row = {k: v for k, v in result.items() if k != "interfaces"}
            rows.append(row)
        except Exception as exc:
            # Record the failure instead of aborting the whole batch.
            # DockQ can fail on malformed PDBs, chain mismatches, etc.
            # Downstream code can filter on capri_class != 'error'.
            rows.append({
                "model_path": mpath,
                "native_path": native_path,
                "DockQ": float("nan"),
                "capri_class": "error",
                "error": str(exc),
            })

    return pd.DataFrame(rows)


def _find_model_files(
    run_dir: str,
    engine: str | None = None,
    glob_pattern: str | None = None,
) -> list[str]:
    """Find docked-model PDB files in *run_dir*.

    Parameters
    ----------
    run_dir : str
        Root directory of docking output.
    engine : str | None
        If given, use engine-specific search patterns.
    glob_pattern : str | None
        Override pattern.

    Returns
    -------
    list[str]
        Sorted list of absolute paths to PDB files.
    """
    if glob_pattern:
        hits = glob.glob(os.path.join(run_dir, glob_pattern), recursive=True)
        return sorted(hits)

    # Auto-detect engine if not given, using the same registry-based
    # detector that collect_scores uses.  Silently fall back to a generic
    # recursive search if detection fails (e.g. unfamiliar directory layout).
    if engine is None:
        try:
            from ppinsight.collect_scores import detect_engine
            engine = detect_engine(run_dir)
        except (ValueError, ImportError):
            engine = None

    # Engine-specific search patterns mirror each engine's output layout.
    if engine == "lightdock":
        # LightDock writes one PDB per pose in each swarm directory.
        return sorted(
            glob.glob(os.path.join(run_dir, "swarm_*", "lightdock_*.pdb"))
        )
    elif engine == "haddock":
        # HADDOCK writes refined models in numbered step directories.
        # Prefer the most-refined stage: emref > flexref > rigidbody.
        for step in ("*_emref", "*_flexref", "*_rigidbody"):
            hits = glob.glob(os.path.join(run_dir, step, "*.pdb"))
            if hits:
                return sorted(hits)
        # Fallback: any PDB anywhere in the run tree
        return sorted(
            glob.glob(os.path.join(run_dir, "**", "*.pdb"), recursive=True)
        )
    elif engine == "rosetta":
        # RosettaDock outputs decoy PDBs into output_files/ by default.
        hits = glob.glob(os.path.join(run_dir, "output_files", "*.pdb"))
        if not hits:
            hits = glob.glob(os.path.join(run_dir, "*.pdb"))
        return sorted(hits)
    else:
        # Unknown engine: recursively find all PDB files as a best guess.
        return sorted(
            glob.glob(os.path.join(run_dir, "**", "*.pdb"), recursive=True)
        )


def _guess_sep(path: str) -> str:
    """Infer CSV vs TSV separator from a file extension."""
    return "\t" if str(path).lower().endswith(".tsv") else ","


def _looks_like_scores_file(path: str) -> bool:
    """Return True when *path* looks like a unified scores table."""
    return os.path.splitext(path)[1].lower() in {".csv", ".tsv"}


def _print_quality_summary(df: pd.DataFrame, *, stream=sys.stdout) -> None:
    """Print the standard CAPRI summary for *df*."""
    counts = capri_summary(df)
    rate = capri_success_rate(df)
    print(f"\n── CAPRI Quality Summary ({len(df)} models) ──", file=stream)
    for cat in CAPRI_ORDER + (["error"] if "error" in counts else []):
        n = counts.get(cat, 0)
        pct = 100 * n / len(df) if len(df) > 0 else 0
        print(f"  {cat:12s}: {n:4d}  ({pct:5.1f}%)", file=stream)
    print(f"  Success rate: {rate:.1%}", file=stream)

    valid = df.dropna(subset=["DockQ"]).sort_values("DockQ", ascending=False)
    if len(valid) > 0:
        print("\n── Top-5 models by DockQ ──", file=stream)
        for _, row in valid.head(5).iterrows():
            path = str(row.get("model_path", ""))
            print(
                f"  {os.path.basename(path):30s}  "
                f"DockQ={row['DockQ']:.4f}  "
                f"({row['capri_class']})",
                file=stream,
            )


def evaluate_scores_dataframe(
    scores_df: pd.DataFrame,
    native_path: str,
    *,
    scores_source: str | None = None,
    engine: str | None = None,
    chain_map: dict[str, str] | None = None,
    capri_peptide: bool = False,
    _load_fn=None,
    _run_fn=None,
) -> pd.DataFrame:
    """Evaluate one DockQ result per unique output_path in *scores_df*."""
    output_col = find_column(scores_df.columns, "output_path")
    if output_col is None:
        raise ValueError(
            "Unified scores file has no 'output_path' column. Re-run "
            "ppinsight collect from run directories so pose paths are "
            "recorded, or use single-model/directory mode instead."
        )

    model_col = find_column(scores_df.columns, "model")
    protein_a_col = find_column(scores_df.columns, "proteinA", "proteina")
    protein_b_col = find_column(scores_df.columns, "proteinB", "proteinb")
    run_id_col = find_column(scores_df.columns, "run_id")
    pose_id_col = find_column(scores_df.columns, "pose_id")
    base_dir = (
        os.path.dirname(os.path.abspath(scores_source))
        if scores_source else None
    )

    pose_table = scores_df.copy()
    if engine is not None and model_col is not None:
        pose_table = pose_table[
            pose_table[model_col].fillna("").astype(str).str.lower()
            == engine.lower()
        ]

    pose_table = pose_table.copy()
    pose_table["_resolved_output_path"] = pose_table[output_col].apply(
        lambda value: resolve_traceability_path(value, base_dir=base_dir)
    )
    pose_table = pose_table[
        pose_table["_resolved_output_path"].astype(str).str.strip().ne("")
    ]

    if pose_table.empty:
        raise ValueError(
            "Unified scores file has no usable output_path values. Re-run "
            "ppinsight collect from run directories so each pose keeps a "
            "resolvable model path."
        )

    pose_table = pose_table.drop_duplicates("_resolved_output_path")

    rows: list[dict[str, Any]] = []
    for _, pose_row in pose_table.iterrows():
        resolved_path = str(pose_row["_resolved_output_path"])
        metadata: dict[str, Any] = {"output_path": resolved_path}
        for logical_name, column_name in (
            ("model", model_col),
            ("proteinA", protein_a_col),
            ("proteinB", protein_b_col),
            ("run_id", run_id_col),
            ("pose_id", pose_id_col),
        ):
            if not column_name:
                continue
            value = pose_row.get(column_name)
            if pd.isna(value):
                continue
            text = str(value).strip()
            if not text or text.lower() in {"nan", "none"}:
                continue
            metadata[logical_name] = value

        if "pose_id" not in metadata:
            metadata["pose_id"] = os.path.splitext(os.path.basename(resolved_path))[0]

        try:
            result = evaluate_complex(
                resolved_path,
                native_path,
                chain_map=chain_map,
                capri_peptide=capri_peptide,
                _load_fn=_load_fn,
                _run_fn=_run_fn,
            )
            row = {k: v for k, v in result.items() if k != "interfaces"}
        except Exception as exc:
            row = {
                "model_path": resolved_path,
                "native_path": native_path,
                "DockQ": float("nan"),
                "capri_class": "error",
                "error": str(exc),
            }
        rows.append({**metadata, **row})

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Convenience: CAPRI summary statistics
# ---------------------------------------------------------------------------

def capri_summary(df: pd.DataFrame) -> dict[str, int]:
    """Count models in each CAPRI quality category.

    Parameters
    ----------
    df : pd.DataFrame
        Must have a ``capri_class`` column (from :func:`evaluate_directory`).

    Returns
    -------
    dict[str, int]
        Counts keyed by CAPRI category, including ``'error'`` if any
        models failed evaluation.
    """
    counts = df["capri_class"].value_counts().to_dict()
    # Ensure all standard categories are present
    for cat in CAPRI_ORDER:
        counts.setdefault(cat, 0)
    return counts


def capri_success_rate(df: pd.DataFrame) -> float:
    """Fraction of models that are at least *acceptable* quality.

    Parameters
    ----------
    df : pd.DataFrame
        Must have a ``capri_class`` column.

    Returns
    -------
    float
        Success rate in [0, 1].  NaN if the DataFrame is empty.
    """
    if len(df) == 0:
        return float("nan")
    acceptable_or_better = df["capri_class"].isin(
        ["acceptable", "medium", "high"]
    )
    return acceptable_or_better.sum() / len(df)


# ---------------------------------------------------------------------------
# Integration: add quality metrics to unified scores
# ---------------------------------------------------------------------------

def add_quality_to_scores(
    scores_df: pd.DataFrame,
    quality_df: pd.DataFrame,
    model_label: str = "",
) -> pd.DataFrame:
    """Append DockQ quality rows to a unified scores DataFrame.

    Converts the wide-format *quality_df* (from :func:`evaluate_directory`)
    into long-format rows compatible with the ``collect_scores`` unified
    schema (``model``, ``score_type``, ``score_value``, ``proteinA``,
    ``proteinB``).  When possible, run/pose traceability is preserved via
    ``run_id``, ``pose_id``, and ``output_path``.

    Parameters
    ----------
    scores_df : pd.DataFrame
        Existing unified scores (may be empty).
    quality_df : pd.DataFrame
        Output from :func:`evaluate_directory`.
    model_label : str
        Value for the ``model`` column (e.g. ``"lightdock"``).

    Returns
    -------
    pd.DataFrame
        Concatenated scores with added DockQ/quality rows.
    """
    quality_metrics = {
        "DockQ": "quality_dockq",
        "fnat": "quality_fnat",
        "iRMSD": "quality_irmsd",
        "LRMSD": "quality_lrmsd",
    }
    rows: list[dict] = []

    score_columns = {
        "model": find_column(scores_df.columns, "model") or "model",
        "score_type": find_column(scores_df.columns, "score_type") or "score_type",
        "score_value": find_column(scores_df.columns, "score_value") or "score_value",
        "proteinA": (
            find_column(scores_df.columns, "proteinA", "proteina")
            or "proteinA"
        ),
        "proteinB": (
            find_column(scores_df.columns, "proteinB", "proteinb")
            or "proteinB"
        ),
        "run_id": find_column(scores_df.columns, "run_id") or "run_id",
        "pose_id": find_column(scores_df.columns, "pose_id") or "pose_id",
        "output_path": find_column(scores_df.columns, "output_path") or "output_path",
    }
    quality_columns = {
        "model": find_column(quality_df.columns, "model"),
        "proteinA": find_column(quality_df.columns, "proteinA", "proteina"),
        "proteinB": find_column(quality_df.columns, "proteinB", "proteinb"),
        "run_id": find_column(quality_df.columns, "run_id"),
        "pose_id": find_column(quality_df.columns, "pose_id"),
        "output_path": find_column(quality_df.columns, "output_path"),
        "model_path": find_column(quality_df.columns, "model_path"),
    }

    def _single_nonempty(column_name: str | None) -> str:
        if not column_name or column_name not in scores_df.columns:
            return ""
        series = scores_df[column_name]
        vals = [
            str(v)
            for v in series.dropna().unique()
            if str(v).strip() and str(v).strip().lower() not in {"nan", "none"}
        ]
        return vals[0] if len(vals) == 1 else ""

    def _quality_value(qrow: pd.Series, key: str) -> str:
        column_name = quality_columns[key]
        if not column_name:
            return ""
        value = qrow.get(column_name)
        if pd.isna(value):
            return ""
        text = str(value).strip()
        if not text or text.lower() in {"nan", "none"}:
            return ""
        return text

    inferred_model = model_label or _single_nonempty(score_columns["model"])
    inferred_protein_a = _single_nonempty(score_columns["proteinA"])
    inferred_protein_b = _single_nonempty(score_columns["proteinB"])
    inferred_run_id = _single_nonempty(score_columns["run_id"])

    for _, qrow in quality_df.iterrows():
        model_path = _quality_value(qrow, "model_path")
        output_path = _quality_value(qrow, "output_path") or model_path
        pose_id = _quality_value(qrow, "pose_id")
        if not pose_id and output_path:
            pose_id = os.path.splitext(os.path.basename(output_path))[0]
        model_value = _quality_value(qrow, "model") or inferred_model or "quality"
        protein_a_value = _quality_value(qrow, "proteinA") or inferred_protein_a
        protein_b_value = _quality_value(qrow, "proteinB") or inferred_protein_b
        run_id_value = _quality_value(qrow, "run_id") or inferred_run_id
        for metric, metric_name in quality_metrics.items():
            val = qrow.get(metric)
            if pd.notna(val):
                row = {
                    score_columns["model"]: model_value,
                    score_columns["score_type"]: metric_name,
                    score_columns["score_value"]: float(val),
                    score_columns["proteinA"]: protein_a_value,
                    score_columns["proteinB"]: protein_b_value,
                }
                if run_id_value:
                    row[score_columns["run_id"]] = run_id_value
                if pose_id:
                    row[score_columns["pose_id"]] = pose_id
                if output_path:
                    row[score_columns["output_path"]] = output_path
                rows.append(row)

    if not rows:
        return scores_df

    base_scores = scores_df
    score_type_col = score_columns["score_type"]
    if score_type_col in base_scores.columns:
        existing_metric_names = {name.lower() for name in quality_metrics.values()}
        base_scores = base_scores[
            ~base_scores[score_type_col].astype(str).str.lower().isin(existing_metric_names)
        ].copy()

    quality_long = pd.DataFrame(rows)
    return pd.concat([base_scores, quality_long], ignore_index=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    """CLI entry point for DockQ quality evaluation."""
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate docking quality by comparing predicted complexes "
            "against a native structure using DockQ. The input can be a "
            "single model PDB, a docking run directory, or a unified "
            "scores TSV/CSV with an output_path column."
        ),
    )
    parser.add_argument(
        "model",
        help=(
            "Path to a docked model PDB, a docking output directory, or a "
            "unified scores TSV/CSV with a non-empty output_path column.  "
            "In directory mode, the engine determines how model PDB files "
            "are found (auto-detected, or set with --engine). In scores-file "
            "mode, DockQ/CAPRI evaluation rows are appended back into the "
            "unified scores output."
        ),
    )
    parser.add_argument(
        "native",
        help=(
            "Path to the native (reference) complex PDB.  This is the "
            "experimentally determined structure against which docked "
            "models are compared."
        ),
    )
    parser.add_argument(
        "--engine",
        choices=["lightdock", "haddock", "rosetta"],
        default=None,
        help=(
            "Docking engine (for directory mode — determines the glob "
            "pattern used to find model PDB files).  Auto-detected from "
            "directory contents if not given.  Set explicitly when "
            "auto-detection fails (e.g. non-standard directory layout).  "
            "In scores-file mode, only rows whose model matches this engine "
            "are evaluated."
        ),
    )
    parser.add_argument(
        "--chain-map",
        default=None,
        help=(
            "Chain mapping as 'modelA:nativeA,modelB:nativeB' (e.g. "
            "'A:A,B:B').  Required when model and native PDBs use "
            "different chain IDs.  If not given, DockQ attempts automatic "
            "chain mapping, which works for most standard complexes."
        ),
    )
    parser.add_argument(
        "--glob",
        default=None,
        dest="glob_pattern",
        help=(
            "Custom glob pattern for finding model PDB files in a "
            "directory (e.g. '*.pdb' or 'refined_*.pdb').  Use when the "
            "default engine-specific pattern doesn't match your file "
            "naming convention."
        ),
    )
    parser.add_argument(
        "--capri-peptide",
        action="store_true",
        help=(
            "Use CAPRI peptide thresholds instead of standard protein "
            "thresholds.  Enable when the ligand is a short peptide "
            "(< ~30 residues) — peptide docking uses looser L-RMSD "
            "cutoffs because small ligands amplify RMSD differences."
        ),
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help=(
            "Write results to a CSV/TSV file (format inferred from "
            "extension). In scores-file mode this writes an augmented "
            "unified scores file with DockQ/CAPRI evaluation rows appended. "
            "Otherwise it writes the per-model quality table."
        ),
    )
    parser.add_argument(
        "--summary",
        action="store_true",
        help=(
            "Print a CAPRI classification summary (counts and percentages "
            "of incorrect/acceptable/medium/high quality models, plus "
            "success rate and top-5 models).  Enabled by default when "
            "--output is not set."
        ),
    )

    args = parser.parse_args(argv)

    # Fail fast with a clear message if DockQ is not installed.
    # This avoids a cryptic traceback when the user runs the CLI
    # without the [quality] extra.
    if not _DOCKQ_AVAILABLE:
        print(
            f"ERROR: {_DOCKQ_IMPORT_ERROR}",
            file=sys.stderr,
        )
        print(
            "Hint: install the quality extra with one of:\n"
            "  pip install 'ppinsight[quality]'\n"
            "  pip install -e '.[quality]'  # from a local checkout",
            file=sys.stderr,
        )
        sys.exit(1)

    # Parse chain map (e.g. "A:A,B:B" → {"A": "A", "B": "B"})
    chain_map = None
    if args.chain_map:
        chain_map = {}
        for pair in args.chain_map.split(","):
            parts = pair.strip().split(":")
            if len(parts) != 2:
                print(
                    f"ERROR: Invalid chain-map entry '{pair}'. "
                    f"Expected 'modelChain:nativeChain'.",
                    file=sys.stderr,
                )
                sys.exit(2)
            chain_map[parts[0].strip()] = parts[1].strip()

    # Single scores file, single structure file, or directory?
    if os.path.isfile(args.model) and _looks_like_scores_file(args.model):
        try:
            scores_df = pd.read_csv(args.model, sep=_guess_sep(args.model))
        except FileNotFoundError:
            print(f"ERROR: scores file not found: {args.model}", file=sys.stderr)
            sys.exit(1)

        quality_df = evaluate_scores_dataframe(
            scores_df,
            args.native,
            scores_source=args.model,
            engine=args.engine,
            chain_map=chain_map,
            capri_peptide=args.capri_peptide,
        )
        merged_scores = add_quality_to_scores(scores_df, quality_df)

        if args.output:
            os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
            merged_scores.to_csv(args.output, sep=_guess_sep(args.output), index=False)
            print(
                (
                    f"Appended DockQ/CAPRI rows for {len(quality_df)} poses "
                    f"to {args.output}"
                ),
                file=sys.stderr,
            )
        else:
            merged_scores.to_csv(sys.stdout, sep="\t", index=False)

        if args.summary or not args.output:
            _print_quality_summary(quality_df, stream=sys.stderr)

    elif os.path.isfile(args.model):
        result = evaluate_complex(
            args.model, args.native,
            chain_map=chain_map,
            capri_peptide=args.capri_peptide,
        )
        # Print results
        print(f"Model:       {result['model_path']}")
        print(f"Native:      {result['native_path']}")
        print(f"DockQ:       {result['DockQ']:.4f}")
        print(f"CAPRI class: {result['capri_class']}")
        print(f"Fnat:        {result.get('fnat', 'N/A')}")
        print(f"I-RMSD:      {result.get('iRMSD', 'N/A')}")
        print(f"L-RMSD:      {result.get('LRMSD', 'N/A')}")
        print(f"F1:          {result.get('F1', 'N/A')}")
        print(f"Clashes:     {result.get('clashes', 'N/A')}")

    elif os.path.isdir(args.model):
        df = evaluate_directory(
            args.model, args.native,
            engine=args.engine,
            chain_map=chain_map,
            glob_pattern=args.glob_pattern,
            capri_peptide=args.capri_peptide,
        )

        if args.output:
            sep = "\t" if args.output.endswith(".tsv") else ","
            df.to_csv(args.output, sep=sep, index=False)
            print(f"Wrote {len(df)} quality evaluations to {args.output}")

        if args.summary or not args.output:
            _print_quality_summary(df)

    else:
        print(f"ERROR: {args.model} is not a file or directory", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
