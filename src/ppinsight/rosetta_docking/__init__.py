"""
Protein-Protein Docking Module

This module provides functionality for protein-protein docking using PyRosetta.

Main components:
- DockingPipeline: High-level interface for complete docking workflow
- prepare_structures: Structure preparation and relaxation
- run_docking: Docking protocol execution
- analyze_scores: Score analysis and ranking

Example usage:
    from ppinsight.docking import DockingPipeline
    
    pipeline = DockingPipeline("protein1.pdb", "protein2.pdb", n_runs=10)
    result = pipeline.run()
    print(f"Final docking score: {result['score']:.2f}")
"""

from .pipeline import DockingPipeline
from .prepare import prepare_structures, relax_structure, combine_proteins
from .dock import run_docking, setup_docking_protocol
from .analyze import analyze_scores, get_top_scores

__all__ = [
    'DockingPipeline',
    'prepare_structures',
    'relax_structure',
    'combine_proteins',
    'run_docking',
    'setup_docking_protocol',
    'analyze_scores',
    'get_top_scores'
]

__version__ = '0.1.0'
