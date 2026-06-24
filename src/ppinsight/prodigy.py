"""prodigy – PRODIGY binding-affinity wrapper for PPInsight.

Wraps the ``prodigy-prot`` package to score docked PDB structures with
PRODIGY and append the predicted ΔG (binding free energy, kcal/mol) and
Kd (dissociation constant, M) to a unified scores DataFrame.

Install the optional dependency with::

    pip install 'ppinsight[prodigy]'

CLI usage::

    ppinsight prodigy scores.tsv --output scores_prodigy.tsv
    ppinsight prodigy scores.tsv --top-n 10 \
        --metric score --engine HADDOCK
    ppinsight prodigy older_scores.tsv --pdb-dir pdbs/ --output scored.tsv
"""

from __future__ import annotations

import inspect
import os
import sys
import warnings
from copy import deepcopy
from pathlib import Path
from typing import Any

import pandas as pd

from ppinsight.utils import find_column, resolve_traceability_path
from ppinsight.visualizer import get_metric_direction


def _require_prodigy() -> Any:
    """Import and return the prodigy_prot Prodigy class.

    Raises a helpful :exc:`ImportError` if ``prodigy-prot`` is not installed.
    """
    try:
        # Older prodigy-prot releases exposed Prodigy under modules.prodigy.
        from prodigy_prot.modules.prodigy import Prodigy  # type: ignore[import]

        return Prodigy
    except ImportError:
        try:
            # Newer prodigy-prot releases expose Prodigy in predict_IC.
            from prodigy_prot.predict_IC import Prodigy  # type: ignore[import]

            return Prodigy
        except ImportError as exc:
            raise ImportError(
                "prodigy-prot is required for binding-affinity scoring but is "
                "not installed.  Install it with:\n\n"
                "    pip install 'ppinsight[prodigy]'\n"
                "    pip install -e '.[prodigy]'  # from a local checkout\n\n"
                "or directly:\n\n"
                "    pip install prodigy-prot\n"
            ) from exc


def _build_prodigy_runner(
    prodigy_class: Any,
    structure: Any,
    structure_name: str,
    selection: str,
    temperature: float,
) -> Any:
    """Instantiate ``Prodigy`` across old and new prodigy-prot APIs."""
    try:
        parameter_names = list(inspect.signature(prodigy_class).parameters)
    except (TypeError, ValueError):
        parameter_names = []

    if "struct_obj" in parameter_names:
        kwargs: dict[str, Any] = {
            "struct_obj": structure,
            "selection": selection,
            "temp": temperature,
        }
        if "name" in parameter_names:
            kwargs["name"] = structure_name
        return prodigy_class(**kwargs)

    if len(parameter_names) >= 4:
        return prodigy_class(structure, structure_name, selection, temp=temperature)

    if len(parameter_names) >= 3:
        return prodigy_class(structure, selection, temp=temperature)

    return prodigy_class(structure, selection=selection, temp=temperature)


def _parse_pdb(pdb_path: str | Path) -> Any:
    """Parse a PDB file and return a Bio.PDB structure."""
    try:
        from Bio.PDB import PDBParser  # type: ignore[import]
    except ImportError as exc:
        raise ImportError(
            "Biopython is required for PDB parsing.  Install with:\n"
            "    conda install biopython\n"
        ) from exc

    parser = PDBParser(QUIET=True)
    path = Path(pdb_path)
    return parser.get_structure(path.stem, str(path))


def _filter_structure_for_prodigy(structure: Any) -> Any:
    """Return a protein-only structure suitable for PRODIGY contact typing."""
    try:
        from Bio.PDB.Polypeptide import is_aa  # type: ignore[import]
        from Bio.PDB.Structure import Structure  # type: ignore[import]
    except ImportError:
        return structure

    if not isinstance(structure, Structure):
        return structure

    filtered = deepcopy(structure)
    for model in filtered:
        for chain in list(model):
            for residue in list(chain):
                keep_residue = (
                    residue.id[0] == " " and is_aa(residue, standard=True)
                )
                if not keep_residue:
                    chain.detach_child(residue.id)
            if not list(chain):
                model.detach_child(chain.id)
    return filtered


def score_pdb(
    pdb_path: str | Path,
    chains: list[str] | None = None,
    temperature: float = 25.0,
) -> dict[str, Any]:
    """Score a single PDB file with PRODIGY.

    Parameters
    ----------
    pdb_path : str or Path
        Path to a PDB file containing a protein-protein complex.
    chains : list of str or None
        Chain IDs to use.  If ``None``, all chains are used.
    temperature : float
        Temperature in °C for Kd prediction (default 25 °C).

    Returns
    -------
    dict
        Keys: ``prodigy_ddg`` (kcal/mol), ``prodigy_kd`` (M),
        ``nis_a``, ``nis_c``, ``n_contacts``, plus PRODIGY contact-type
        bin counts (``CC``, ``CP``, ``AC``, ``AA``, ``PP``, ``AP``).
        On failure the dict contains an ``"error"`` key and NaN values
        for ``prodigy_ddg`` and ``prodigy_kd``.
    """
    pdb_path = Path(pdb_path)

    if not pdb_path.exists():
        return {
            "prodigy_ddg": float("nan"),
            "prodigy_kd": float("nan"),
            "error": f"file_not_found: {pdb_path}",
        }

    Prodigy = _require_prodigy()

    try:
        structure = _parse_pdb(pdb_path)
        structure = _filter_structure_for_prodigy(structure)
    except Exception as exc:
        return {
            "prodigy_ddg": float("nan"),
            "prodigy_kd": float("nan"),
            "error": f"parse_error: {exc}",
        }

    if chains:
        selection = list(chains)
    else:
        seen: set[str] = set()
        unique_chains: list[str] = []
        for model in structure:
            for ch in model.get_chains():
                if ch.id not in seen:
                    seen.add(ch.id)
                    unique_chains.append(ch.id)
        selection = unique_chains

    try:
        runner = _build_prodigy_runner(
            Prodigy,
            structure,
            pdb_path.stem,
            selection,
            temperature,
        )
        runner.predict(distance_cutoff=5.5, acc_threshold=0.05)
    except ValueError as exc:
        msg = str(exc).lower()
        if "no contacts" in msg or "no inter-chain" in msg:
            return {
                "prodigy_ddg": float("nan"),
                "prodigy_kd": float("nan"),
                "error": "no_contacts",
            }
        return {
            "prodigy_ddg": float("nan"),
            "prodigy_kd": float("nan"),
            "error": f"prodigy_error: {exc}",
        }
    except Exception as exc:
        return {
            "prodigy_ddg": float("nan"),
            "prodigy_kd": float("nan"),
            "error": f"prodigy_error: {exc}",
        }

    bins = getattr(runner, "bins", {})
    return {
        "prodigy_ddg": (
            float(runner.ba_val) if runner.ba_val is not None else float("nan")
        ),
        "prodigy_kd": (
            float(runner.kd_val) if runner.kd_val is not None else float("nan")
        ),
        "nis_a": (
            float(runner.nis_a) if runner.nis_a is not None else float("nan")
        ),
        "nis_c": (
            float(runner.nis_c) if runner.nis_c is not None else float("nan")
        ),
        "n_contacts": sum(bins.values()) if bins else 0,
        "CC": float(bins.get("CC", 0)),
        "CP": float(bins.get("CP", 0)),
        "AC": float(bins.get("AC", 0)),
        "AA": float(bins.get("AA", 0)),
        "PP": float(bins.get("PP", 0)),
        "AP": float(bins.get("AP", 0)),
    }


def score_directory(
    pdb_dir: str | Path,
    pattern: str = "*.pdb",
    chains: list[str] | None = None,
    temperature: float = 25.0,
    verbose: bool = False,
) -> pd.DataFrame:
    """Score all PDB files matching *pattern* in *pdb_dir*.

    Parameters
    ----------
    pdb_dir : str or Path
        Directory containing PDB files.
    pattern : str
        Glob pattern (default ``"*.pdb"``).
    chains : list of str or None
        Chain IDs (passed to :func:`score_pdb`).
    temperature : float
        Temperature in °C for Kd prediction.
    verbose : bool
        Print progress to stderr.

    Returns
    -------
    DataFrame
        One row per PDB file.
    """
    pdb_dir = Path(pdb_dir)
    pdb_files = sorted(pdb_dir.glob(pattern))
    if not pdb_files:
        warnings.warn(
            f"No files matching '{pattern}' found in {pdb_dir}",
            stacklevel=2,
        )
        return pd.DataFrame()

    rows = []
    for pdb in pdb_files:
        if verbose:
            print(f"  Scoring {pdb.name} …", file=sys.stderr)
        result = dict(score_pdb(pdb, chains=chains, temperature=temperature))
        result["pdb"] = pdb.name
        rows.append(result)

    df = pd.DataFrame(rows)
    cols = ["pdb"] + [c for c in df.columns if c != "pdb"]
    return df[cols].reset_index(drop=True)


def _guess_sep(path: str | Path) -> str:
    """Infer CSV vs TSV separator from a file extension."""
    return "\t" if str(path).lower().endswith(".tsv") else ","


def _resolve_pose_paths(
    scores_df: pd.DataFrame,
    *,
    scores_source: str | Path | None = None,
    pdb_dir: str | Path | None = None,
) -> pd.Series:
    """Resolve one model path per row using output_path or pdb_dir+pdb."""
    output_col = find_column(scores_df.columns, "output_path")
    pdb_col = find_column(scores_df.columns, "pdb")
    base_dir = Path(scores_source).resolve().parent if scores_source else None
    pdb_root = Path(pdb_dir).resolve() if pdb_dir is not None else None

    resolved_paths: list[str] = []
    for _, row in scores_df.iterrows():
        resolved = ""
        if output_col is not None:
            resolved = resolve_traceability_path(
                row.get(output_col),
                base_dir=base_dir,
            )

        if not resolved and pdb_col is not None and pdb_root is not None:
            pdb_name = row.get(pdb_col)
            if pd.notna(pdb_name):
                text = str(pdb_name).strip()
                if text and text.lower() not in {"-", "nan", "none"}:
                    resolved = str((pdb_root / text).resolve())

        resolved_paths.append(resolved)

    return pd.Series(resolved_paths, index=scores_df.index, dtype="object")


def add_prodigy_to_scores(
    scores_df: pd.DataFrame,
    pdb_dir: str | Path | None = None,
    model_col: str = "model",
    top_n: int | None = None,
    metric: str | None = None,
    engine: str | None = None,
    chains: list[str] | None = None,
    temperature: float = 25.0,
    verbose: bool = False,
    scores_source: str | Path | None = None,
) -> pd.DataFrame:
    """Append PRODIGY ΔG and Kd columns to a unified scores DataFrame.

    Parameters
    ----------
    scores_df : DataFrame
        Unified scores (long format).  Prefer a non-empty ``output_path``
        column from ``ppinsight collect``.
    pdb_dir : str or Path or None
        Fallback directory for scores tables that only store a ``pdb``
        filename column instead of full ``output_path`` values.
    model_col : str
        Column that identifies the docking engine.
    top_n : int or None
        If set, only score the top *top_n* poses per model+pair group.
        Requires *metric* to be specified so that the selection direction
        (best = highest or lowest value) can be determined from
        :data:`~ppinsight.visualizer.METRIC_METADATA`.
    metric : str or None
        The ``score_type`` value to rank poses by when *top_n* is set.
        The sort direction is determined automatically via
        :func:`~ppinsight.visualizer.get_metric_direction`:
        ``nlargest`` for higher-is-better metrics (e.g. ``luciferin_score``),
        ``nsmallest`` for lower-is-better metrics (e.g. ``score`` for HADDOCK).
        Must be provided when *top_n* is set.
    engine : str or None
        If set, only score rows where *model_col* equals *engine*.
    chains : list of str or None
        Chain IDs (passed to :func:`score_pdb`).
    temperature : float
        Temperature in °C for Kd prediction.
    verbose : bool
        Print progress to stderr.
    scores_source : str or Path or None
        Path to the scores file. Used to resolve relative ``output_path``
        values when present.

    Returns
    -------
    DataFrame
        Copy of *scores_df* with ``prodigy_ddg`` and ``prodigy_kd`` added.
    """
    df = scores_df.copy()
    model_column = (
        find_column(df.columns, model_col)
        or find_column(df.columns, "model")
    )
    score_type_col = find_column(df.columns, "score_type")
    score_value_col = find_column(df.columns, "score_value")

    if "prodigy_ddg" not in df.columns:
        df["prodigy_ddg"] = float("nan")
    if "prodigy_kd" not in df.columns:
        df["prodigy_kd"] = float("nan")

    mask = pd.Series([True] * len(df), index=df.index)
    if engine is not None and model_column is not None:
        mask &= df[model_column].fillna("").astype(str).str.lower() == engine.lower()

    subset = df[mask].copy()
    subset["_resolved_model_path"] = _resolve_pose_paths(
        subset,
        scores_source=scores_source,
        pdb_dir=pdb_dir,
    )
    if subset["_resolved_model_path"].astype(str).str.strip().eq("").all():
        warnings.warn(
            "scores_df has no usable model paths. Prefer a non-empty "
            "'output_path' column from ppinsight collect, or provide "
            "--pdb-dir together with a 'pdb' column for external scores files.",
            stacklevel=2,
        )
        return df

    if top_n is not None:
        if metric is None:
            raise ValueError(
                "metric must be specified when top_n is set.  "
                "Pass the score_type to rank poses by (e.g. 'score' for HADDOCK, "
                "'luciferin_score' for LightDock) so that the correct sort direction "
                "(nsmallest vs nlargest) can be determined automatically."
            )
        # Filter to rows for the chosen metric only, so that mixed long-format
        # files (multiple score_type rows per pose) do not distort ranking.
        if score_type_col is not None:
            ranking_rows = subset[subset[score_type_col] == metric]
        else:
            ranking_rows = subset
        higher_is_better = get_metric_direction(metric)
        if score_value_col is None:
            raise ValueError(
                "scores_df has no score_value column, so top_n ranking cannot "
                "determine the best poses."
            )
        group_cols = [
            col
            for col in (
                model_column,
                find_column(ranking_rows.columns, "proteinA", "proteina"),
                find_column(ranking_rows.columns, "proteinB", "proteinb"),
            )
            if col is not None
        ]
        if not group_cols:
            raise ValueError(
                "scores_df must contain a model column to use top_n pose ranking."
            )
        grouped = ranking_rows.groupby(group_cols, group_keys=False)
        try:
            if higher_is_better:
                top_rows = grouped.apply(
                    lambda g: g.nlargest(top_n, score_value_col),
                    include_groups=False,
                )
            else:
                top_rows = grouped.apply(
                    lambda g: g.nsmallest(top_n, score_value_col),
                    include_groups=False,
                )
        except TypeError:
            # Pandas < 2.2 does not support include_groups.
            if higher_is_better:
                top_rows = grouped.apply(lambda g: g.nlargest(top_n, score_value_col))
            else:
                top_rows = grouped.apply(lambda g: g.nsmallest(top_n, score_value_col))
        # Restrict scoring to only the PDB files that appear in the top poses.
        top_paths = set(top_rows["_resolved_model_path"].dropna().unique())
        subset = subset[subset["_resolved_model_path"].isin(top_paths)]

    scored: dict[str, dict[str, float]] = {}
    for model_path in subset["_resolved_model_path"].dropna().unique():
        path_text = str(model_path).strip()
        if not path_text:
            continue
        if verbose:
            print(f"  Scoring {os.path.basename(path_text)} …", file=sys.stderr)
        scored[path_text] = score_pdb(
            path_text, chains=chains, temperature=temperature
        )

    for idx in subset.index:
        model_path = subset.at[idx, "_resolved_model_path"]
        if model_path in scored:
            r = scored[model_path]
            df.at[idx, "prodigy_ddg"] = r.get("prodigy_ddg", float("nan"))
            df.at[idx, "prodigy_kd"] = r.get("prodigy_kd", float("nan"))

    return df


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _build_parser():
    import argparse

    parser = argparse.ArgumentParser(
        prog="ppinsight prodigy",
        description=(
            "Score docked PDB structures with PRODIGY and append predicted\n"
            "binding free energy (\u0394G) and dissociation constant (Kd) to a\n"
            "unified scores TSV/CSV file.\n\n"
            "Requires: pip install 'ppinsight[prodigy]'"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "scores",
        help=(
            "Path to unified scores TSV/CSV. Prefer files from ppinsight "
            "collect with a non-empty output_path column. Files that only "
            "store bare pdb filenames still work with --pdb-dir."
        ),
    )
    parser.add_argument(
        "--pdb-dir",
        metavar="DIR",
        help=(
            "Fallback directory containing docked PDB files. Use this only "
            "for scores tables that have a pdb column but no output_path."
        ),
    )
    parser.add_argument(
        "--output",
        "-o",
        default=None,
        help="Output file path.  If omitted, writes to stdout.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=None,
        metavar="N",
        help=(
            "Only score the top N poses per model+pair group.  "
            "Requires --metric.  Sort direction is determined automatically: "
            "nlargest for higher-is-better metrics (e.g. luciferin_score), "
            "nsmallest for lower-is-better metrics (e.g. score for HADDOCK)."
        ),
    )
    parser.add_argument(
        "--metric",
        default=None,
        metavar="SCORE_TYPE",
        help=(
            "score_type value to rank poses by when --top-n is used "
            "(e.g. 'score' for HADDOCK, 'luciferin_score' for LightDock).  "
            "Required with --top-n."
        ),
    )
    parser.add_argument(
        "--engine",
        default=None,
        help="Only score poses from this docking engine (model name).",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=25.0,
        metavar="°C",
        help="Temperature for Kd prediction (default: 25.0 °C).",
    )
    parser.add_argument(
        "--summary",
        action="store_true",
        help="Print a per-engine summary table (mean ΔG) to stderr.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print per-file progress to stderr.",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    """CLI entry point for ``ppinsight prodigy``."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    _require_prodigy()

    if args.top_n is not None and args.metric is None:
        parser.error("--metric is required when --top-n is used")

    sep = _guess_sep(args.scores)
    try:
        scores_df = pd.read_csv(args.scores, sep=sep)
    except FileNotFoundError:
        print(
            f"ERROR: scores file not found: {args.scores}", file=sys.stderr
        )
        sys.exit(1)

    result = add_prodigy_to_scores(
        scores_df,
        pdb_dir=args.pdb_dir,
        top_n=args.top_n,
        metric=args.metric,
        engine=args.engine,
        temperature=args.temperature,
        verbose=args.verbose,
        scores_source=args.scores,
    )

    if args.summary:
        ddg_col = result["prodigy_ddg"]
        if not ddg_col.isna().all() and "model" in result.columns:
            summary = (
                result.groupby("model")["prodigy_ddg"]
                .agg(["mean", "std", "count"])
                .rename(
                    columns={
                        "mean": "mean_ddg",
                        "std": "std_ddg",
                        "count": "n_scored",
                    }
                )
                .sort_values("mean_ddg")
            )
            print("\nPRODIGY summary (ΔG kcal/mol):", file=sys.stderr)
            print(summary.to_string(), file=sys.stderr)
            print(file=sys.stderr)

    out_sep = "\t" if (args.output or "").lower().endswith(".tsv") else ","
    if args.output:
        os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
        result.to_csv(args.output, sep=out_sep, index=False)
        print(f"Saved to {args.output}", file=sys.stderr)
    else:
        result.to_csv(sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    main()
