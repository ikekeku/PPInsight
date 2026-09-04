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
import importlib.util
import os
import shutil
import sys
import tempfile
import time
from collections.abc import Callable
from pathlib import Path

import pandas as pd

from ppinsight import registry

_DEFAULT_BATCH_CORES = 1
_DEFAULT_BATCH_LIGHTDOCK_STEPS = 100
_DEFAULT_BATCH_LIGHTDOCK_SWARMS = 400
_DEFAULT_BATCH_LIGHTDOCK_GLOWWORMS = 200
_DEFAULT_BATCH_HADDOCK_SAMPLING = 10000
_DEFAULT_BATCH_HADDOCK_SELECT_TOP = 400
_DEFAULT_BATCH_HADDOCK_TOLERANCE = 5
_DEFAULT_BATCH_ROSETTA_N_RUNS = 5000
_DEFAULT_BATCH_ROSETTA_TOP_N = 20
_DEFAULT_BATCH_ROSETTA_CLUSTER_TOP_N = 200
_DEFAULT_BATCH_ROSETTA_RMSD_CUTOFF = 4.0

_SCREENING_LIGHTDOCK_STEPS = 50
_SCREENING_LIGHTDOCK_SWARMS = 50
_SCREENING_LIGHTDOCK_GLOWWORMS = 50
_SCREENING_HADDOCK_SAMPLING = 1000
_SCREENING_HADDOCK_SELECT_TOP = 100
_SCREENING_ROSETTA_N_RUNS = 100
_SCREENING_ROSETTA_TOP_N = 20
_SCREENING_ROSETTA_CLUSTER_TOP_N = 100

_RESULT_KEY_COLUMNS = ["proteinA", "proteinB", "engine"]

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


def _option_was_supplied(argv: list[str], *options: str) -> bool:
    """Return whether any option was explicitly supplied on the command line."""
    return any(
        argument == option or argument.startswith(f"{option}=")
        for argument in argv
        for option in options
    )


def _apply_screening_preset(args, argv: list[str]) -> None:
    """Apply low-cost screening values unless a specific flag overrides them."""
    if not args.screening:
        return

    if not _option_was_supplied(argv, "--lightdock-steps"):
        args.lightdock_steps = _SCREENING_LIGHTDOCK_STEPS
    if not _option_was_supplied(argv, "--lightdock-swarms"):
        args.lightdock_swarms = _SCREENING_LIGHTDOCK_SWARMS
    if not _option_was_supplied(argv, "--lightdock-glowworms"):
        args.lightdock_glowworms = _SCREENING_LIGHTDOCK_GLOWWORMS
    if not _option_was_supplied(
        argv, "--lightdock-anm", "--lightdock-no-anm"
    ):
        args.lightdock_anm = False
    if not _option_was_supplied(argv, "--lightdock-auto-clean-pdb"):
        args.lightdock_auto_clean_pdb = True

    if not _option_was_supplied(argv, "--haddock-sampling"):
        args.haddock_sampling = _SCREENING_HADDOCK_SAMPLING
    if not _option_was_supplied(argv, "--haddock-select-top"):
        args.haddock_select_top = _SCREENING_HADDOCK_SELECT_TOP
    if not _option_was_supplied(
        argv,
        "--haddock-skip-refinement",
        "--haddock-skip-flexref",
        "--haddock-skip-emref",
    ):
        args.haddock_skip_refinement = True

    if not _option_was_supplied(argv, "--rosetta-n-runs"):
        args.rosetta_n_runs = _SCREENING_ROSETTA_N_RUNS
    if not _option_was_supplied(argv, "--rosetta-top-n"):
        args.rosetta_top_n = _SCREENING_ROSETTA_TOP_N
    if not _option_was_supplied(argv, "--rosetta-cluster-top-n"):
        args.rosetta_cluster_top_n = _SCREENING_ROSETTA_CLUSTER_TOP_N
    if not _option_was_supplied(argv, "--rosetta-relax", "--rosetta-no-relax"):
        args.rosetta_relax = False


def _format_elapsed_time(seconds: float) -> str:
    """Format an elapsed duration as a compact wall-clock string."""
    total_seconds = max(0, int(round(seconds)))
    hours, remainder = divmod(total_seconds, 3600)
    minutes, remaining_seconds = divmod(remainder, 60)

    if hours:
        return f"{hours}h {minutes}m {remaining_seconds}s"
    if minutes:
        return f"{minutes}m {remaining_seconds}s"
    return f"{remaining_seconds}s"



def _result_row(
    protein_a: str,
    protein_b: str,
    label: str,
    family: str,
    engine: str,
    output_dir: str,
    status: str,
    *,
    error_type: str = "",
    error_message: str = "",
    log_path: str = "",
    preflight_warnings: str = "",
) -> dict[str, str]:
    """Build one stable-schema batch result row."""
    return {
        "proteinA": protein_a,
        "proteinB": protein_b,
        "label": label,
        "family": family,
        "engine": engine,
        "output_dir": output_dir,
        "status": status,
        "error_type": error_type,
        "error_message": error_message,
        "log_path": log_path,
        "preflight_warnings": preflight_warnings,
    }


def _log_path_for_run_dir(output_dir: str | None) -> str:
    """Return a known engine log path when one exists for *output_dir*."""
    if not output_dir:
        return ""
    run_dir = Path(output_dir)
    candidates = [run_dir / "log", run_dir / "lightdock.log", run_dir / "rosetta.log"]
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate)
    return ""


def _upsert_results(existing: pd.DataFrame, incoming: pd.DataFrame) -> pd.DataFrame:
    """Replace manifest rows by pair and engine, retaining one current row."""
    if incoming.empty:
        return existing.copy()
    incoming = incoming.drop_duplicates(_RESULT_KEY_COLUMNS, keep="last")
    if existing.empty:
        return incoming.copy().reset_index(drop=True)

    existing = existing.drop_duplicates(_RESULT_KEY_COLUMNS, keep="last")
    keys = incoming[_RESULT_KEY_COLUMNS].drop_duplicates()
    existing_index = existing.set_index(_RESULT_KEY_COLUMNS).index
    incoming_index = pd.MultiIndex.from_frame(keys)
    retained = existing.loc[~existing_index.isin(incoming_index)]
    return pd.concat([retained, incoming], ignore_index=True)


def _is_safe_output_dir(output_dir: str, output_root: str) -> bool:
    """Return whether a candidate run path is a child of the output root."""
    try:
        Path(output_dir).resolve().relative_to(Path(output_root).resolve())
    except ValueError:
        return False
    return True


def _remove_failed_run_dir(output_dir: str, output_root: str) -> bool:
    """Safely remove one failed run directory recorded under *output_root*."""
    if not output_dir or not _is_safe_output_dir(output_dir, output_root):
        return False
    run_dir = Path(output_dir)
    if not run_dir.is_dir():
        return False
    shutil.rmtree(run_dir)
    return True


def _pdb_coordinate_count(pdb_path: str) -> int:
    """Return the number of coordinate records in a PDB file."""
    with open(pdb_path, encoding="utf-8") as handle:
        return sum(
            line.startswith(("ATOM  ", "HETATM"))
            for line in handle
        )


def _pdb_hetatm_count(pdb_path: str) -> int:
    """Return the number of HETATM records in a PDB file."""
    with open(pdb_path, encoding="utf-8") as handle:
        return sum(line.startswith("HETATM") for line in handle)


def _preflight_engine(
    engine: str,
    receptor_pdb: str,
    ligand_pdb: str,
    engine_options: dict,
) -> tuple[str, str, list[str]]:
    """Validate one engine's dependencies and staged input compatibility."""
    warnings: list[str] = []

    for pdb_path in (receptor_pdb, ligand_pdb):
        if _pdb_coordinate_count(pdb_path) == 0:
            return (
                "pdb_invalid",
                f"No ATOM or HETATM records found in {pdb_path}.",
                warnings,
            )

    if engine == "lightdock":
        required = ("lightdock3_setup.py", "lightdock3.py")
        missing = [
            executable
            for executable in required
            if shutil.which(executable) is None
        ]
        if missing:
            return (
                "engine_unavailable",
                "LightDock executable(s) not found on PATH: " + ", ".join(missing),
                warnings,
            )
        hetero_atoms = _pdb_hetatm_count(receptor_pdb) + _pdb_hetatm_count(ligand_pdb)
        if hetero_atoms:
            message = (
                f"Found {hetero_atoms} HETATM record(s); LightDock DFIRE may "
                "reject unsupported residues."
            )
            if engine_options.get("auto_clean_pdb", False):
                warnings.append(message + " Auto-clean retry is enabled.")
            else:
                warnings.append(message + " Use --lightdock-auto-clean-pdb.")
        return "", "", warnings

    if engine == "haddock":
        if shutil.which("haddock3") is None:
            return (
                "engine_unavailable",
                "The haddock3 executable was not found on PATH.",
                warnings,
            )
        from ppinsight.pdb_to_haddock import copy_inputs

        with tempfile.TemporaryDirectory(
            prefix="ppinsight-haddock-preflight-"
        ) as temp_dir:
            staged_dir = Path(temp_dir) / "data"
            staged_dir.mkdir()
            try:
                staged_rec, staged_lig, _ = copy_inputs(
                    staged_dir,
                    receptor_pdb,
                    ligand_pdb,
                    auto_filter=bool(engine_options.get("auto_filter", True)),
                )
            except Exception as exc:
                return type(exc).__name__, str(exc), warnings

            staged_sets = []
            for pdb_path in (staged_rec, staged_lig):
                identifiers = set()
                with open(pdb_path, encoding="utf-8") as handle:
                    for line in handle:
                        if not line.startswith(("ATOM  ", "HETATM")):
                            continue
                        identifiers.add(line[72:76].strip() or line[21].strip())
                if not identifiers:
                    return (
                        "pdb_invalid",
                        f"No staged coordinates in {pdb_path}.",
                        warnings,
                    )
                staged_sets.append(identifiers)
            if staged_sets[0] & staged_sets[1]:
                return (
                    "haddock_chain_collision",
                    "HADDOCK staged partners share chain/seg identifiers: "
                    + ", ".join(sorted(staged_sets[0] & staged_sets[1])),
                    warnings,
                )
        return "", "", warnings

    if engine == "rosetta":
        if importlib.util.find_spec("pyrosetta") is None:
            return (
                "engine_unavailable",
                "PyRosetta is not installed in the active Python environment.",
                warnings,
            )
        return "", "", warnings

    return "unknown_engine", f"No preflight implementation for '{engine}'.", warnings


def preflight_batch(
    pairs_df: pd.DataFrame,
    engines: list[str],
    pdb_dir: str | None = None,
    limit: int | None = None,
    engine_kwargs: dict[str, dict] | None = None,
) -> pd.DataFrame:
    """Validate pairs and engine prerequisites without running docking jobs."""
    engine_kwargs = engine_kwargs or {}
    df = pairs_df.copy()
    df.columns = [column.strip() for column in df.columns]
    if limit:
        df = df.head(limit)

    results: list[dict[str, str]] = []
    for _, row in df.iterrows():
        protein_a = str(row["proteinA"]).strip()
        protein_b = str(row["proteinB"]).strip()
        label = str(row.get("label", "")).strip()
        family = str(row.get("family", "")).strip()
        receptor_pdb = _resolve_pdb_for_protein(protein_a, pdb_dir)
        ligand_pdb = _resolve_pdb_for_protein(protein_b, pdb_dir)

        for engine in engines:
            if receptor_pdb is None or ligand_pdb is None:
                missing = protein_a if receptor_pdb is None else protein_b
                results.append(_result_row(
                    protein_a, protein_b, label, family, engine, "",
                    "preflight_failed", error_type="pdb_missing",
                    error_message=f"PDB not found for {missing}.",
                ))
                continue

            error_type, error_message, warnings = _preflight_engine(
                engine,
                receptor_pdb,
                ligand_pdb,
                engine_kwargs.get(engine, {}),
            )
            status = "preflight_failed" if error_type else "preflight_ok"
            results.append(_result_row(
                protein_a,
                protein_b,
                label,
                family,
                engine,
                "",
                status,
                error_type=error_type,
                error_message=error_message,
                preflight_warnings="; ".join(warnings),
            ))
    return pd.DataFrame(results)



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
    engine_kwargs: dict[str, dict] | None = None,
    completed_runs: set[tuple[str, str, str]] | None = None,
    failed_run_dirs: dict[tuple[str, str, str], str] | None = None,
    clean_failed: bool = False,
    on_result: Callable[[dict], None] | None = None,
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
    engine_kwargs : dict[str, dict] | None
        Optional per-engine keyword arguments forwarded to each engine's
        registered runner.
    completed_runs : set[tuple[str, str, str]] | None
        Successful ``(proteinA, proteinB, engine)`` keys to skip. Intended for
        resumable batch execution.
    failed_run_dirs : dict[tuple[str, str, str], str] | None
        Failed-run directories from an existing results manifest.
    clean_failed : bool
        Remove a safely recorded failed-run directory before retrying it.
    on_result : callable | None
        Optional callback invoked for every newly recorded result row.

    Returns
    -------
    pd.DataFrame
        A results table with columns: ``proteinA``, ``proteinB``, ``label``,
        ``engine``, ``output_dir``, ``status``.
    """
    if output_root is None:
        from ppinsight.utils import _project_root
        output_root = os.path.join(_project_root(), "data", "output")
    output_root = os.path.abspath(os.path.expanduser(output_root))

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
    engine_kwargs = engine_kwargs or {}
    completed_runs = completed_runs or set()
    failed_run_dirs = failed_run_dirs or {}

    def record(result: dict) -> None:
        """Store one result and persist it through the optional callback."""
        results.append(result)
        if on_result is not None:
            on_result(result)

    for i, row in df.iterrows():
        pA = str(row["proteinA"]).strip()
        pB = str(row["proteinB"]).strip()
        label = str(row.get("label", "")).strip()
        family = str(row.get("family", "")).strip()

        print(f"\n[{i+1}/{total}] {pA} vs {pB} (label={label})")

        if dry_run:
            for eng in engines:
                print(f"  [dry-run] Would run {eng}")
                record(_result_row(
                    pA, pB, label, family, eng, "(dry-run)", "dry_run"
                ))
            continue

        # Resolve PDB files
        rec_pdb = _resolve_pdb_for_protein(pA, pdb_dir)
        lig_pdb = _resolve_pdb_for_protein(pB, pdb_dir)

        if not rec_pdb:
            print(f"  WARNING: No PDB found for {pA} — skipping")
            for eng in engines:
                record(_result_row(
                    pA, pB, label, family, eng, "", "pdb_missing_A",
                    error_type="pdb_missing",
                    error_message=f"PDB not found for {pA}.",
                ))
            continue
        if not lig_pdb:
            print(f"  WARNING: No PDB found for {pB} — skipping")
            for eng in engines:
                record(_result_row(
                    pA, pB, label, family, eng, "", "pdb_missing_B",
                    error_type="pdb_missing",
                    error_message=f"PDB not found for {pB}.",
                ))
            continue

        print(f"  Receptor: {rec_pdb}")
        print(f"  Ligand:   {lig_pdb}")

        for eng in engines:
            if (pA, pB, eng) in completed_runs:
                print(f"  [{eng}] already completed; skipping (--resume)")
                continue

            failed_run_dir = failed_run_dirs.get((pA, pB, eng), "")
            if clean_failed and failed_run_dir:
                if _remove_failed_run_dir(failed_run_dir, output_root):
                    print(f"  [{eng}] removed prior failed run: {failed_run_dir}")
                else:
                    print(
                        f"  [{eng}] did not remove prior failed run "
                        f"(missing or outside output root): {failed_run_dir}"
                    )

            try:
                plugin = registry.get(eng)
            except KeyError:
                print(f"  WARNING: Unknown engine '{eng}' — skipping")
                record(_result_row(
                    pA, pB, label, family, eng, "", "unknown_engine",
                    error_type="unknown_engine",
                    error_message=f"Unknown engine '{eng}'.",
                ))
                continue

            runner = plugin.runner
            if runner is None:
                print(f"  WARNING: Engine '{eng}' has no runner — skipping")
                record(_result_row(
                    pA, pB, label, family, eng, "", "no_runner",
                    error_type="no_runner",
                    error_message=f"Engine '{eng}' has no runner.",
                ))
                continue

            print(f"  Running {eng}...")
            t0 = time.time()
            run_kwargs = dict(engine_kwargs.get(eng, {}))
            try:
                out_dir = runner(rec_pdb, lig_pdb, output_root, label, **run_kwargs)
            except registry.EngineRunError as exc:
                elapsed = time.time() - t0
                print(f"  [{eng}] failed ({elapsed:.1f}s)")
                record(_result_row(
                    pA,
                    pB,
                    label,
                    family,
                    eng,
                    exc.output_dir or "",
                    "failed",
                    error_type=exc.error_type,
                    error_message=str(exc),
                    log_path=exc.log_path or _log_path_for_run_dir(exc.output_dir),
                ))
                continue
            except Exception as exc:
                elapsed = time.time() - t0
                print(f"  [{eng}] failed ({elapsed:.1f}s)")
                record(_result_row(
                    pA,
                    pB,
                    label,
                    family,
                    eng,
                    "",
                    "failed",
                    error_type=type(exc).__name__,
                    error_message=str(exc),
                ))
                continue
            elapsed = time.time() - t0

            status = "ok" if out_dir else "failed"
            print(f"  [{eng}] {status} ({elapsed:.1f}s)")
            if out_dir:
                record(_result_row(
                    pA, pB, label, family, eng, str(out_dir), "ok"
                ))
            else:
                record(_result_row(
                    pA,
                    pB,
                    label,
                    family,
                    eng,
                    "",
                    "failed",
                    error_type="runner_returned_none",
                    error_message=f"{eng} runner returned no output directory.",
                ))

    return pd.DataFrame(results)


def _engine_kwargs_from_args(args) -> dict[str, dict]:
    """Translate CLI flags into per-engine kwargs for registry runners."""
    haddock_skip_refinement = bool(args.haddock_skip_refinement)
    haddock_skip_flexref = bool(args.haddock_skip_flexref or haddock_skip_refinement)
    # emref depends on flexref outputs; skipping flexref implies skipping emref.
    haddock_skip_emref = bool(
        args.haddock_skip_emref
        or haddock_skip_flexref
        or haddock_skip_refinement
    )

    return {
        "lightdock": {
            "steps": args.lightdock_steps,
            "swarms": args.lightdock_swarms,
            "glowworms": args.lightdock_glowworms,
            "cores": args.cores,
            "anm": args.lightdock_anm,
            "scoring": args.lightdock_scoring,
            "auto_clean_pdb": args.lightdock_auto_clean_pdb,
        },
        "haddock": {
            "ncores": args.cores,
            "sampling": args.haddock_sampling,
            "select_top": args.haddock_select_top,
            "tolerance": args.haddock_tolerance,
            "skip_flexref": haddock_skip_flexref,
            "skip_emref": haddock_skip_emref,
        },
        "rosetta": {
            "n_runs": args.rosetta_n_runs,
            "top_n": args.rosetta_top_n,
            "relax": args.rosetta_relax,
            "cluster": not args.rosetta_no_cluster,
            "cluster_top_n": args.rosetta_cluster_top_n,
            "rmsd_cutoff": args.rosetta_rmsd_cutoff,
            "auto_filter": not args.rosetta_no_auto_filter,
            "pyrosetta_debug": args.rosetta_debug_pyrosetta,
            "save_top": args.rosetta_save_top,
        },
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    """CLI: batch-run docking for all pairs in a pairs file."""
    batch_started_at = time.perf_counter()
    parser = argparse.ArgumentParser(
        description=(
            "Run docking pipelines for every pair in a pairs file.  "
            "Produces a batch results table listing each run directory."
        ),
        epilog=(
            "Batch defaults when flags are omitted:\n"
            "  - LightDock: steps=100, swarms=400, glowworms=200, cores=1, "
            "ANM=enabled, scoring=LightDock default.\n"
            "  - HADDOCK: runs by default in batch with rigidbody sampling=10000, "
            "seletop select=400, tolerance=5, full refinement enabled, and "
            "ncores follows --cores.\n"
            "  - Rosetta: available through --engines rosetta (requires PyRosetta), "
            "with n_runs=5000, top_n=20, relax=on, clustering=on.\n"
            "Use --screening for a reduced, explicitly low-cost preset."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
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
        "--cores",
        type=int,
        default=_DEFAULT_BATCH_CORES,
        help=(
            "CPU cores to pass to supported engine runners (default: 1).  "
            "This maps to LightDock simulation cores and HADDOCK ncores. "
            "Rosetta trajectories run serially in the current pipeline."
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
            "Root directory for engine run directories and the default results "
            "table (default: data/output).  "
            "Batch creates engine-specific folders under this root (for "
            "example lightdock_runs/, haddock_runs/, rosetta_runs/), with one "
            "subdirectory per pair × engine run. Use -o to place the results "
            "table elsewhere."
        ),
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help=(
            "Output results table "
            "(default: <output-root>/scores/batch_results.csv).  Contains "
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
    parser.add_argument(
        "--resume",
        action="store_true",
        help=(
            "Reuse successful rows from the output results table. Failed or "
            "missing-input rows retry; each new result replaces its prior row."
        ),
    )
    parser.add_argument(
        "--clean-failed",
        action="store_true",
        help=(
            "With --resume, remove each safely recorded failed run directory "
            "under --output-root before retrying that pair and engine."
        ),
    )
    parser.add_argument(
        "--preflight",
        action="store_true",
        help=(
            "Validate PDB coordinates, engine availability, LightDock "
            "unsupported-residue risk, and HADDOCK staged chain/seg IDs "
            "without launching docking jobs."
        ),
    )
    parser.add_argument(
        "--screening",
        action="store_true",
        help=(
            "Use a reduced end-to-end preset for checking that pairs run: "
            "LightDock 50 steps/50 swarms/50 glowworms without ANM; HADDOCK "
            "1000 rigidbody models, select 100, no refinement; Rosetta 100 "
            "trajectories without FastRelax. Numeric and string engine flags "
            "override the corresponding preset value; some boolean preset "
            "flags (e.g. --lightdock-auto-clean-pdb, --haddock-skip-"
            "refinement) have no complementary negation flag."
        ),
    )

    lightdock_group = parser.add_argument_group("LightDock batch options")
    lightdock_group.add_argument(
        "--lightdock-steps",
        type=int,
        default=_DEFAULT_BATCH_LIGHTDOCK_STEPS,
        help=(
            "LightDock optimisation steps in batch mode (default: 100).  "
            "Use 10 for quick smoke tests; increase toward 100+ for broader "
            "sampling."
        ),
    )
    lightdock_group.add_argument(
        "--lightdock-swarms",
        type=int,
        default=_DEFAULT_BATCH_LIGHTDOCK_SWARMS,
        help=(
            "Override LightDock swarm count (default: 400)."
        ),
    )
    lightdock_group.add_argument(
        "--lightdock-glowworms",
        type=int,
        default=_DEFAULT_BATCH_LIGHTDOCK_GLOWWORMS,
        help=(
            "Override LightDock glowworms-per-swarm (default: 200)."
        ),
    )
    lightdock_group.add_argument(
        "--lightdock-scoring",
        default=None,
        help=(
            "LightDock scoring function (default: LightDock default scoring "
            "function)."
        ),
    )
    lightdock_group.add_argument(
        "--lightdock-anm",
        action="store_true",
        default=True,
        help=(
            "Enable ANM flexibility for LightDock (the production default)."
        ),
    )
    lightdock_group.add_argument(
        "--lightdock-no-anm",
        action="store_false",
        dest="lightdock_anm",
        help=(
            "Disable ANM flexibility for LightDock to reduce memory use."
        ),
    )
    lightdock_group.add_argument(
        "--lightdock-auto-clean-pdb",
        action="store_true",
        help=(
            "If LightDock scoring rejects unsupported non-protein residues, "
            "auto-generate protein-only copies under the run directory and "
            "retry that pair."
        ),
    )

    haddock_group = parser.add_argument_group("HADDOCK batch options")
    haddock_group.add_argument(
        "--haddock-sampling",
        type=int,
        default=_DEFAULT_BATCH_HADDOCK_SAMPLING,
        help=(
            "HADDOCK rigidbody sampling count in generated configs "
            "(default: 10000)."
        ),
    )
    haddock_group.add_argument(
        "--haddock-select-top",
        type=int,
        default=_DEFAULT_BATCH_HADDOCK_SELECT_TOP,
        help=(
            "HADDOCK seletop count in generated configs (default: 400)."
        ),
    )
    haddock_group.add_argument(
        "--haddock-tolerance",
        type=int,
        default=_DEFAULT_BATCH_HADDOCK_TOLERANCE,
        help=(
            "HADDOCK module output-fault tolerance percentage for rigidbody, "
            "flexref, and emref (default: 5)."
        ),
    )
    haddock_group.add_argument(
        "--haddock-skip-refinement",
        action="store_true",
        help=(
            "Skip HADDOCK refinement stages (flexref + emref) and continue "
            "from rigid-body outputs. Useful for brittle pairs and smoke tests."
        ),
    )
    haddock_group.add_argument(
        "--haddock-skip-flexref",
        action="store_true",
        help=(
            "Skip HADDOCK flexref stage. This also disables emref because "
            "emref depends on flexref outputs."
        ),
    )
    haddock_group.add_argument(
        "--haddock-skip-emref",
        action="store_true",
        help=(
            "Skip HADDOCK emref (water refinement) stage."
        ),
    )

    rosetta_group = parser.add_argument_group("Rosetta batch options")
    rosetta_group.add_argument(
        "--rosetta-n-runs",
        type=int,
        default=_DEFAULT_BATCH_ROSETTA_N_RUNS,
        help=(
            "Number of Rosetta docking trajectories per pair (default: 5000)."
        ),
    )
    rosetta_group.add_argument(
        "--rosetta-top-n",
        type=int,
        default=_DEFAULT_BATCH_ROSETTA_TOP_N,
        help=(
            "Rosetta top-N decoys used for final I_sc averaging "
            f"(default: {_DEFAULT_BATCH_ROSETTA_TOP_N})."
        ),
    )
    rosetta_group.add_argument(
        "--rosetta-relax",
        action="store_true",
        default=True,
        help=(
            "Enable Rosetta FastRelax preprocessing (the production default)."
        ),
    )
    rosetta_group.add_argument(
        "--rosetta-no-relax",
        action="store_false",
        dest="rosetta_relax",
        help=(
            "Skip Rosetta FastRelax preprocessing to shorten a screening run."
        ),
    )
    rosetta_group.add_argument(
        "--rosetta-no-cluster",
        action="store_true",
        help=(
            "Skip Rosetta decoy clustering in batch mode.  Clustering is "
            "enabled by default."
        ),
    )
    rosetta_group.add_argument(
        "--rosetta-cluster-top-n",
        type=int,
        default=_DEFAULT_BATCH_ROSETTA_CLUSTER_TOP_N,
        help=(
            "Top Rosetta decoys to include in clustering (default: 200)."
        ),
    )
    rosetta_group.add_argument(
        "--rosetta-rmsd-cutoff",
        type=float,
        default=_DEFAULT_BATCH_ROSETTA_RMSD_CUTOFF,
        help=(
            "Rosetta clustering Cα-RMSD cutoff in Å (default: 4.0)."
        ),
    )
    rosetta_group.add_argument(
        "--rosetta-no-auto-filter",
        action="store_true",
        help=(
            "Disable Rosetta DBREF accession-chain auto-filtering in batch "
            "mode."
        ),
    )
    rosetta_group.add_argument(
        "--rosetta-debug-pyrosetta",
        action="store_true",
        help=(
            "Enable verbose PyRosetta tracer output during Rosetta batch runs."
        ),
    )
    rosetta_group.add_argument(
        "--rosetta-save-top",
        type=int,
        default=0,
        help=(
            "Save top N Rosetta decoy structures per pair (default: 0)."
        ),
    )

    argv = list(sys.argv[1:] if argv is None else argv)
    args = parser.parse_args(argv)
    _apply_screening_preset(args, argv)

    if args.cores < 1:
        parser.error("--cores must be at least 1")
    if args.clean_failed and not args.resume:
        parser.error("--clean-failed requires --resume")

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

        # consrank isn't a registry EnginePlugin (it reranks a pool of poses
        # gathered from other engines rather than running/parsing one engine's
        # run directory), but its binary availability is checked the same
        # way, so surface it here for a single place to check tool readiness.
        from ppinsight.consrank import (
            consrank_binary_available,
            default_consrank_binary,
        )
        bin_ok = "✓" if consrank_binary_available() else "✗"
        print(
            f"  {'consrank':15s}  binary={bin_ok}       "
            f"Reference-free consensus ranking (Iter-CONSRANK) "
            f"[{default_consrank_binary()}]"
        )
        return

    if not args.pairs:
        parser.error("the following arguments are required: pairs")

    # Read pairs file
    sep = "\t" if args.pairs.endswith(".tsv") else ","
    pairs_df = pd.read_csv(args.pairs, sep=sep)

    from ppinsight.utils import _project_root

    output_root = args.output_root or os.path.join(_project_root(), "data", "output")
    output_root = os.path.abspath(os.path.expanduser(output_root))
    default_result_name = (
        "preflight_results.csv" if args.preflight else "batch_results.csv"
    )
    output_path = args.output or os.path.join(
        output_root, "scores", default_result_name
    )
    output_path = os.path.abspath(os.path.expanduser(output_path))
    out_sep = "\t" if output_path.endswith(".tsv") else ","

    if args.preflight:
        results_df = preflight_batch(
            pairs_df,
            engines=args.engines,
            pdb_dir=args.pdb_dir,
            limit=args.limit,
            engine_kwargs=_engine_kwargs_from_args(args),
        )
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        results_df.to_csv(output_path, sep=out_sep, index=False)
        print(f"\nWrote {len(results_df)} preflight results to {output_path}")
        print("\n── Preflight summary ──")
        if results_df.empty:
            print("No pairs to validate.")
        else:
            print(results_df["status"].value_counts().to_string())
        print(
            "Total preflight runtime: "
            f"{_format_elapsed_time(time.perf_counter() - batch_started_at)}"
        )
        return

    existing_results = pd.DataFrame()
    completed_runs: set[tuple[str, str, str]] = set()
    failed_run_dirs: dict[tuple[str, str, str], str] = {}

    if args.resume and os.path.isfile(output_path):
        existing_results = pd.read_csv(output_path, sep=out_sep)
        required_columns = {"proteinA", "proteinB", "engine", "status"}
        if not required_columns.issubset(existing_results.columns):
            parser.error(
                "--resume requires an existing PPInsight results table with "
                f"columns: {sorted(required_columns)}"
            )
        existing_results = _upsert_results(
            existing_results.iloc[0:0],
            existing_results,
        )
        completed_rows = existing_results[existing_results["status"] == "ok"]
        completed_runs = {
            (str(row.proteinA), str(row.proteinB), str(row.engine))
            for row in completed_rows.itertuples(index=False)
        }
        print(
            "Resuming "
            f"{len(completed_runs)} successful engine run(s) from {output_path}"
        )
        failed_rows = existing_results[existing_results["status"] == "failed"]
        failed_run_dirs = {
            (str(row.proteinA), str(row.proteinB), str(row.engine)):
            str(row.output_dir)
            for row in failed_rows.itertuples(index=False)
            if isinstance(row.output_dir, str) and row.output_dir
        }

    new_results: list[dict] = []
    _flush_counter = [0]
    _FLUSH_EVERY = 10

    def _flush_to_disk(force: bool = False) -> None:
        nonlocal existing_results
        _flush_counter[0] += 1
        if not force and _flush_counter[0] % _FLUSH_EVERY != 0:
            return
        if not new_results:
            return
        combined = _upsert_results(
            existing_results,
            pd.DataFrame(new_results),
        )
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        temporary_path = f"{output_path}.tmp"
        combined.to_csv(temporary_path, sep=out_sep, index=False)
        os.replace(temporary_path, output_path)
        existing_results = combined
        new_results.clear()

    def persist_result(result: dict) -> None:
        """Buffer a completed engine result and periodically flush to disk."""
        new_results.append(result)
        _flush_to_disk()

    results_df = batch_dock(
        pairs_df,
        engines=args.engines,
        pdb_dir=args.pdb_dir,
        output_root=output_root,
        limit=args.limit,
        dry_run=args.dry_run,
        engine_kwargs=_engine_kwargs_from_args(args),
        completed_runs=completed_runs,
        failed_run_dirs=failed_run_dirs,
        clean_failed=args.clean_failed,
        on_result=persist_result,
    )
    if new_results:
        _flush_to_disk(force=True)

    combined_results = _upsert_results(
        existing_results,
        results_df,
    )
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    combined_results.to_csv(output_path, sep=out_sep, index=False)
    print(f"\nWrote {len(combined_results)} results to {output_path}")

    # Summary
    print("\n── Summary ──")
    if len(results_df) > 0:
        print(results_df["status"].value_counts().to_string())
    else:
        print("No new engine runs.")
    print(
        "Total batch runtime: "
        f"{_format_elapsed_time(time.perf_counter() - batch_started_at)}"
    )


if __name__ == "__main__":
    main()
