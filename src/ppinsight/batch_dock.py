"""
batch_dock – run docking pipelines for every pair in a pairs file.

This is the main batch-mode entry point.  It reads a ``pairs.csv`` (as
produced by :mod:`ppinsight.parse_pairs`) and, for each row, runs the
requested docking engines.  It writes a **batch results table** that
records which (pair, engine) runs succeeded and where their output
directories were written.

Use :mod:`ppinsight.collect_scores` after batch docking to turn those
engine output directories into a unified ``scores.tsv`` that the
visualizer can consume.

Usage::

    # Parse the annotation table into pairs.csv first
    parse_pairs "RTK Interactome.tsv" -o pairs.csv

    # Then run batch docking (lightdock only, dry-run to verify)
    batch_dock pairs.csv --engines lightdock --dry-run

    # Full run with all engines
    batch_dock pairs.csv --engines lightdock haddock --pdb-dir pdb_files/ \
        -o batch_results.csv

    # Limit to first 5 pairs for a quick test
    batch_dock pairs.csv --engines lightdock --limit 5 -o batch_results.csv
"""

import argparse
import os
import time

import pandas as pd

from ppinsight import registry

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _resolve_pdb_for_protein(name: str, pdb_dir: str | None = None) -> str | None:
    """Try to find a PDB file for *name* in *pdb_dir* or the repo.

    Returns the path if found, else ``None``.
    """
    from ppinsight.utils import resolve_input_path
    try:
        return resolve_input_path(name, search_root=pdb_dir)
    except FileNotFoundError:
        return None


def _fetch_pdb_for_uniprot(gene_name: str, pdb_dir: str) -> str | None:
    """Use protein_fetch to download a PDB for *gene_name*.

    This is a best-effort lookup: UniProt gene names ≠ accession IDs, so
    we search UniProt for the gene name and use the first hit.

    Currently returns None (auto-fetch not yet implemented).
    Users should pre-fetch PDBs with ``ppinsight fetch`` or supply a
    ``--pdb-dir`` containing the needed files.

    Returns the PDB path if successful, else ``None``.
    """
    # TODO: implement auto-fetch via protein_fetch module
    return None



# ---------------------------------------------------------------------------
# Core batch function
# ---------------------------------------------------------------------------

def batch_dock(
    pairs_df: pd.DataFrame,
    engines: list[str],
    pdb_dir: str | None = None,
    output_root: str | None = None,
    limit: int | None = None,
    dry_run: bool = False,
) -> pd.DataFrame:
    """Run docking for every pair in *pairs_df*.

    Parameters
    ----------
    pairs_df : DataFrame
        Must have columns ``proteinA``, ``proteinB``.  Optional columns:
        ``label`` (interaction/non-interaction), ``family``, ``references``.
    engines : list[str]
        Docking engines to run (``"lightdock"``, ``"haddock"``, ``"rosetta"``).
    pdb_dir : str | None
        Directory containing pre-fetched PDB files.
    output_root : str | None
        Root directory for docking outputs.
    limit : int | None
        Only process the first *limit* pairs (useful for testing).
    dry_run : bool
        If True, don't actually run docking — just report what would happen.

    Returns
    -------
    pd.DataFrame
        A results table with columns: ``proteinA``, ``proteinB``, ``label``,
        ``engine``, ``output_dir``, ``status``.
    """
    if output_root is None:
        from ppinsight.utils import _project_root
        output_root = os.path.join(_project_root(), "data", "output")

    # Normalize column names for matching — handles whitespace and casing
    # differences from different spreadsheet exports.
    norm_cols = {c.lower().replace(" ", "") for c in pairs_df.columns}
    if not {"proteina", "proteinb"}.issubset(norm_cols):
        raise ValueError(
            f"pairs_df must have 'proteinA' and 'proteinB' columns. "
            f"Got: {list(pairs_df.columns)}"
        )

    df = pairs_df.copy()
    # Normalize column names
    df.columns = [c.strip() for c in df.columns]

    if limit:
        df = df.head(limit)

    results: list[dict] = []
    total = len(df)

    for i, row in df.iterrows():
        pA = str(row["proteinA"]).strip()
        pB = str(row["proteinB"]).strip()
        label = str(row.get("label", "")).strip()
        family = str(row.get("family", "")).strip()

        print(f"\n[{i+1}/{total}] {pA} vs {pB} (label={label})")

        if dry_run:
            for eng in engines:
                print(f"  [dry-run] Would run {eng}")
                results.append({
                    "proteinA": pA, "proteinB": pB, "label": label,
                    "family": family, "engine": eng,
                    "output_dir": "(dry-run)", "status": "dry_run",
                })
            continue

        # Resolve PDB files
        rec_pdb = _resolve_pdb_for_protein(pA, pdb_dir)
        lig_pdb = _resolve_pdb_for_protein(pB, pdb_dir)

        if not rec_pdb:
            print(f"  WARNING: No PDB found for {pA} — skipping")
            for eng in engines:
                results.append({
                    "proteinA": pA, "proteinB": pB, "label": label,
                    "family": family, "engine": eng,
                    "output_dir": "", "status": "pdb_missing_A",
                })
            continue
        if not lig_pdb:
            print(f"  WARNING: No PDB found for {pB} — skipping")
            for eng in engines:
                results.append({
                    "proteinA": pA, "proteinB": pB, "label": label,
                    "family": family, "engine": eng,
                    "output_dir": "", "status": "pdb_missing_B",
                })
            continue

        print(f"  Receptor: {rec_pdb}")
        print(f"  Ligand:   {lig_pdb}")

        for eng in engines:
            try:
                plugin = registry.get(eng)
            except KeyError:
                print(f"  WARNING: Unknown engine '{eng}' — skipping")
                results.append({
                    "proteinA": pA, "proteinB": pB, "label": label,
                    "family": family, "engine": eng,
                    "output_dir": "", "status": "unknown_engine",
                })
                continue

            runner = plugin.runner
            if runner is None:
                print(f"  WARNING: Engine '{eng}' has no runner — skipping")
                results.append({
                    "proteinA": pA, "proteinB": pB, "label": label,
                    "family": family, "engine": eng,
                    "output_dir": "", "status": "no_runner",
                })
                continue

            print(f"  Running {eng}...")
            t0 = time.time()
            out_dir = runner(rec_pdb, lig_pdb, output_root, label)
            elapsed = time.time() - t0

            status = "ok" if out_dir else "failed"
            print(f"  [{eng}] {status} ({elapsed:.1f}s)")
            results.append({
                "proteinA": pA, "proteinB": pB, "label": label,
                "family": family, "engine": eng,
                "output_dir": out_dir or "", "status": status,
            })

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    """CLI: batch-run docking for all pairs in a pairs file."""
    parser = argparse.ArgumentParser(
        description=(
            "Run docking pipelines for every pair in a pairs file.  "
            "Produces a batch results table listing each run directory."
        ),
    )
    parser.add_argument(
        "pairs",
        nargs="?",
        default=None,
        help=(
            "Path to a pairs CSV/TSV (from 'ppinsight parse' or manual "
            "creation).  Each row must have proteinA and proteinB columns.  "
            "An optional 'label' column is preserved in the output."
        ),
    )
    parser.add_argument(
        "--engines",
        nargs="+",
        default=["lightdock"],
        help=(
            "Docking engine(s) to run (default: lightdock).  Specify "
            "multiple engines to benchmark them head-to-head on the same "
            "pairs (e.g. --engines lightdock haddock rosetta).  Each engine "
            "must be registered in the plugin registry."
        ),
    )
    parser.add_argument(
        "--list-engines",
        action="store_true",
        help=(
            "List all registered docking engines (name, runner/parser "
            "availability, description) and exit.  Use this to discover "
            "which engines are installed and available."
        ),
    )
    parser.add_argument(
        "--pdb-dir",
        default=None,
        help=(
            "Directory containing pre-fetched PDB files.  Batch mode uses "
            "this directory (and local path resolution) to find receptor/"
            "ligand structures.  Automatic UniProt/PDB fetch fallback is "
            "not implemented here yet, so pre-fetch with 'ppinsight fetch' "
            "or provide --pdb-dir for reliable runs."
        ),
    )
    parser.add_argument(
        "--output-root",
        default=None,
        help=(
            "Root directory for docking outputs (default: auto-generated).  "
            "Each pair × engine gets a subdirectory under this root."
        ),
    )
    parser.add_argument(
        "-o", "--output",
        default=os.path.join("data", "output", "scores", "batch_results.csv"),
        help=(
            "Output results table "
            "(default: data/output/scores/batch_results.csv).  Contains "
            "one row per (pair, engine) with status (success/error) and "
            "paths to output directories.  Parent directories are created "
            "automatically."
        ),
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help=(
            "Only process the first N pairs (default: all).  Useful for "
            "testing the pipeline on a small subset before committing to "
            "a full batch run."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Don't actually run docking — just print what would happen.  "
            "Use this to verify pair resolution, PDB availability, and "
            "engine configuration before launching expensive compute."
        ),
    )

    args = parser.parse_args(argv)

    # --list-engines: show registered engines and exit
    if args.list_engines:
        engines = registry.list_engines()
        print("Registered docking engines:")
        for name in engines:
            plugin = registry.get(name)
            runner_ok = "✓" if plugin.runner else "✗"
            parser_ok = "✓" if plugin.parser else "✗"
            desc = plugin.description or "(no description)"
            print(f"  {name:15s}  runner={runner_ok}  parser={parser_ok}  {desc}")
        return

    if not args.pairs:
        parser.error("the following arguments are required: pairs")

    # Read pairs file
    sep = "\t" if args.pairs.endswith(".tsv") else ","
    pairs_df = pd.read_csv(args.pairs, sep=sep)

    results_df = batch_dock(
        pairs_df,
        engines=args.engines,
        pdb_dir=args.pdb_dir,
        output_root=args.output_root,
        limit=args.limit,
        dry_run=args.dry_run,
    )

    # Write results
    out_sep = "\t" if args.output.endswith(".tsv") else ","
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    results_df.to_csv(args.output, sep=out_sep, index=False)
    print(f"\nWrote {len(results_df)} results to {args.output}")

    # Summary
    if len(results_df) > 0:
        print("\n── Summary ──")
        print(results_df["status"].value_counts().to_string())


if __name__ == "__main__":
    main()
