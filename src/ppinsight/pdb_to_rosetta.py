"""
CLI helper for running PyRosetta protein-protein docking.

Usage matches the style of the other PPInsight pipeline scripts::

    pdb_to_rosetta 2UUY_rec 2UUY_lig
    pdb_to_rosetta 2UUY_rec 2UUY_lig --n-runs 5 --no-relax --save-top 3
"""

import argparse
import os
import sys
from contextlib import contextmanager
from pathlib import Path

# Shared path resolver
from ppinsight.utils import resolve_input_path


@contextmanager
def _pyrosetta_installer_env():
    """Ensure pyrosetta-installer shells out through this Python env."""
    env_bin_dir = os.path.dirname(sys.executable)
    original = {
        "PATH": os.environ.get("PATH"),
        "PYTHONNOUSERSITE": os.environ.get("PYTHONNOUSERSITE"),
        "PYTHONPATH": os.environ.get("PYTHONPATH"),
        "PIP_USER": os.environ.get("PIP_USER"),
    }

    os.environ["PATH"] = env_bin_dir + os.pathsep + (original["PATH"] or "")
    os.environ["PYTHONNOUSERSITE"] = "1"
    os.environ.pop("PYTHONPATH", None)
    os.environ.pop("PIP_USER", None)

    try:
        yield
    finally:
        for key, value in original.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def _ensure_pyrosetta():
    """Import PyRosetta, auto-installing via pyrosetta-installer if missing.

    PyRosetta is a large optional dependency (~500 MB).  The installer
    package (``pyrosetta-installer``) handles the platform-specific wheel
    download, but it may fail on unsupported platforms (e.g. ARM Linux).
    """
    try:
        import pyrosetta  # noqa: F401
    except ImportError:
        print("PyRosetta not found — installing via pyrosetta-installer …")
        try:
            import pyrosetta_installer
            with _pyrosetta_installer_env():
                pyrosetta_installer.install_pyrosetta(skip_if_installed=False)
            import pyrosetta  # noqa: F401
        except Exception as exc:
            print(
                f"ERROR: Could not auto-install PyRosetta: {exc}\n"
                "Install manually:  python -c "
                "\"import pyrosetta_installer; "
                "pyrosetta_installer.install_pyrosetta()\""
            )
            sys.exit(1)


def _make_output_dir(receptor_pdb, ligand_pdb,
                     base_root="data/output",
                     method="rosetta_runs"):
    """Create an organised output directory and return its path.

    Mirrors the LightDock/HADDOCK convention:
    ``<base_root>/<method>/<receptor_stem>_vs_<ligand_stem>/``
    """
    from ppinsight.pdb_to_lightdock import _project_root

    proj = _project_root()
    if not os.path.isabs(base_root):
        base_root = os.path.join(proj, base_root)

    rec_name = Path(receptor_pdb).stem
    lig_name = Path(ligand_pdb).stem
    run_base = os.path.join(base_root, method, f"{rec_name}_vs_{lig_name}")
    run_dir = run_base
    if os.path.exists(run_dir):
        idx = 1
        while True:
            candidate = f"{run_base}_{idx}"
            if not os.path.exists(candidate):
                run_dir = candidate
                break
            idx += 1
    os.makedirs(run_dir, exist_ok=True)
    return run_dir


# ── CLI ──────────────────────────────────────────────────────────

def main(argv=None):
    """CLI entrypoint for the Rosetta docking helper script."""
    parser = argparse.ArgumentParser(
        description="Run PyRosetta protein-protein docking for a receptor-ligand pair."
    )
    parser.add_argument("receptor", help="Path (or basename) of the receptor PDB file")
    parser.add_argument("ligand", help="Path (or basename) of the ligand PDB file")
    parser.add_argument(
        "--n-runs", type=int, default=10,
        help=(
            "Number of independent docking trajectories (default: 10).  "
            "Rosetta's stochastic search needs many runs to sample the "
            "energy landscape.  Use 10–100 for quick screening; the "
            "Rosetta docs recommend 10,000–100,000 for production global "
            "docking.  More runs = better coverage but linear cost increase."
        ),
    )
    parser.add_argument(
        "--top-n", type=int, default=20,
        help=(
            "Number of top-scoring decoys to average for the final reported "
            "score (default: 20).  Averaging smooths out stochastic noise.  "
            "Raise for high-n-runs benchmarks; lower (e.g. 5) for few-run "
            "quick screens."
        ),
    )
    parser.add_argument(
        "--no-relax", action="store_true",
        help=(
            "Skip FastRelax energy minimisation before docking.  Relaxation "
            "removes clashes and improves starting energies (recommended for "
            "crystal structures).  Skip only for speed during debugging or "
            "when structures are already relaxed.  Skipping may introduce "
            "artefactual high energies that degrade scoring."
        ),
    )
    parser.add_argument(
        "--save-top", type=int, default=0,
        help=(
            "Save the top N docked structures as PDB files (default: 0 = "
            "don't save).  Enable when you need to visualise or further "
            "analyse the best poses (e.g. in PyMOL or for DockQ evaluation)."
        ),
    )
    parser.add_argument(
        "--output-dir", default=None,
        help=(
            "Directory for output files (default: auto-generated under "
            "data/output/rosetta_runs/).  Repeated runs for the same pair "
            "auto-increment the folder name to avoid overwriting.  Set this "
            "explicitly to organise outputs when running many pairs or "
            "benchmarking."
        ),
    )
    parser.add_argument(
        "-q", "--quiet", action="store_true",
        help=(
            "Suppress progress messages.  Useful in batch/scripted runs "
            "where only the final score line matters."
        ),
    )
    parser.add_argument(
        "--input-dir", default=None,
        help=(
            "Directory to search when receptor/ligand are basenames instead "
            "of full paths (default: repo root).  Useful when PDB files "
            "live in a shared directory outside the project tree."
        ),
    )
    parser.add_argument(
        "--no-auto-filter",
        action="store_true",
        help=(
            "Disable PPInsight's accession-based DBREF chain filtering.  "
            "By default, accession-named mixed-complex PDBs are reduced to "
            "the chains mapped to that accession before Rosetta preparation.  "
            "Use this only when you intentionally want the full deposited "
            "complex or a non-accession partner definition."
        ),
    )
    parser.add_argument(
        "--no-cluster", action="store_true",
        help=(
            "Skip decoy clustering after docking.  Clustering groups decoys "
            "by Cα-RMSD to identify distinct binding modes — the best-scoring "
            "member of the largest cluster is the recommended prediction "
            "(standard Rosetta best practice).  Skip only when scipy or "
            "PyRosetta is unavailable, or for quick debugging."
        ),
    )
    parser.add_argument(
        "--cluster-top-n", type=int, default=200,
        help=(
            "Number of top-scoring decoys to cluster (default: 200).  Only "
            "the best decoys are worth clustering — low-scoring outliers add "
            "noise.  Raise for very large n-runs; lower for small test runs."
        ),
    )
    parser.add_argument(
        "--rmsd-cutoff", type=float, default=4.0,
        help=(
            "Cα-RMSD cutoff (Å) for hierarchical clustering (default: 4.0).  "
            "Decoys within this RMSD of each other are grouped together.  "
            "Tighter cutoffs (e.g. 2.0) produce more, smaller clusters — "
            "useful for high-resolution discrimination.  Looser cutoffs "
            "(e.g. 6.0) merge more poses — useful for coarse exploration."
        ),
    )

    args = parser.parse_args(argv)

    # ── resolve PDB paths ────────────────────────────────────────
    try:
        receptor = resolve_input_path(args.receptor, search_root=args.input_dir)
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        print("Hint: use --input-dir to specify the directory containing "
              "your PDB files.", file=sys.stderr)
        sys.exit(2)

    try:
        ligand = resolve_input_path(args.ligand, search_root=args.input_dir)
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        print("Hint: use --input-dir to specify the directory containing "
              "your PDB files.", file=sys.stderr)
        sys.exit(2)

    # ── ensure PyRosetta is available ────────────────────────────
    _ensure_pyrosetta()

    # Lazy import so the CLI arg-parsing stays fast even if PyRosetta
    # takes a moment to load.
    from ppinsight.docking import DockingPipeline

    # ── output directory ─────────────────────────────────────────
    output_dir = args.output_dir or _make_output_dir(receptor, ligand)
    verbose = not args.quiet

    if verbose:
        print(f"Receptor : {receptor}")
        print(f"Ligand   : {ligand}")
        print(f"Runs     : {args.n_runs}")
        print(f"Output   : {output_dir}")

    # ── run the pipeline ─────────────────────────────────────────
    pipeline = DockingPipeline(
        receptor, ligand,
        n_runs=args.n_runs,
        top_n=args.top_n,
        relax=not args.no_relax,
        verbose=verbose,
        auto_filter=not args.no_auto_filter,
    )
    result = pipeline.run()

    # ── optional outputs ─────────────────────────────────────────
    if args.save_top > 0:
        pipeline.save_top_structures(output_dir, top_n=args.save_top)

    # Always save per-decoy scores — needed for downstream clustering
    # and for 'ppinsight collect' to parse Rosetta results.
    csv_path = os.path.join(output_dir, "docking_scores.csv")
    pipeline.save_scores(csv_path)

    # ── optional clustering ──────────────────────────────────────
    # After docking, cluster decoys by Cα-RMSD to identify distinct
    # binding modes.  The best-scoring member of the largest cluster
    # is the recommended prediction (standard Rosetta best practice).
    if not args.no_cluster:
        if os.path.isfile(csv_path):
            try:
                from ppinsight.rosetta.analyze import cluster_and_rank
                clustered = cluster_and_rank(
                    csv_path,
                    output_dir,
                    top_n=args.cluster_top_n,
                    rmsd_cutoff=args.rmsd_cutoff,
                )
                if clustered is not None and verbose:
                    n_clusters = clustered["cluster"].nunique()
                    print(f"✓ Clustered decoys into {n_clusters} group(s)")
            except Exception as exc:
                if verbose:
                    print(f"⚠ Clustering skipped: {exc}")

    # ── summary ──────────────────────────────────────────────────
    if verbose:
        pipeline.print_summary()

    print(f"\nFinal docking score: {result['final_score']:.2f}")
    return result


if __name__ == "__main__":
    main()
