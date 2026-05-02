"""
Protein-Protein Docking Module

This module provides functionality for protein-protein docking using PyRosetta.

Main components:
- DockingPipeline: High-level interface for complete docking workflow
- prepare_structures: Structure preparation and relaxation
- run_docking: Docking protocol execution
- analyze_scores: Score analysis and ranking

Example usage:
    from ppinsight.rosetta import DockingPipeline

    pipeline = DockingPipeline("protein1.pdb", "protein2.pdb", n_runs=10)
    result = pipeline.run()
    print(f"Final docking score: {result['final_score']:.2f}")

Note:
    Importing this sub-package requires PyRosetta.  Install it with::

        python -c "import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()"
"""

try:
    # Eagerly import Rosetta sub-modules so that ``from ppinsight.rosetta
    # import DockingPipeline`` works.  If PyRosetta is missing, fall
    # through to the except block and set a placeholder.
    from .analyze import analyze_scores, cluster_and_rank, get_top_scores
    from .cluster import best_of_largest_cluster, cluster_decoys
    from .dock import run_docking, setup_docking_protocol
    from .pipeline import DockingPipeline
    from .prepare_structure import combine_proteins, prepare_structures, relax_structure
except ImportError:
    # PyRosetta is not installed — provide a helpful error on attribute access
    import warnings as _warnings
    _warnings.warn(
        "ppinsight.rosetta could not import PyRosetta. "
        "Rosetta docking features are unavailable. "
        "Install with: python -c "
        "\"import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()\"",
        ImportWarning,
        stacklevel=2,
    )
    DockingPipeline = None  # type: ignore[assignment,misc]

__all__ = [
    'DockingPipeline',
    'prepare_structures',
    'relax_structure',
    'combine_proteins',
    'run_docking',
    'setup_docking_protocol',
    'analyze_scores',
    'get_top_scores',
    'cluster_and_rank',
    'cluster_decoys',
    'best_of_largest_cluster',
]

__version__ = '0.1.0'
