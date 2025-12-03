"""
Docking execution module.

This module handles:
- Setting up docking protocols
- Running multiple docking simulations
- Managing docking movers and score functions
"""

import sys
from typing import List

try:
    import pyrosetta
    from pyrosetta import rosetta
except ImportError:
    print("ERROR: PyRosetta not found!")
    print("Please install PyRosetta: pip install pyrosetta-*.whl")
    sys.exit(1)


def setup_docking_protocol(low_res_cycles=50, high_res_cycles=100, 
                          randomize1=True, randomize2=True,
                          local_refine=False):
    """
    Setup a docking protocol mover.
    
    Args:
        low_res_cycles: Number of low-resolution (centroid) docking cycles
        high_res_cycles: Number of high-resolution refinement cycles
        randomize1: Randomize first rigid body orientation
        randomize2: Randomize second rigid body orientation
        local_refine: If True, perform local refinement only (no global search)
    
    Returns:
        Configured DockingProtocol mover
    """
    # Create docking protocol
    docking = rosetta.protocols.docking.DockingProtocol()
    
    # Set up docking parameters
    docking.set_low_res_protocol_only(False)
    docking.set_docking_local_refine(local_refine)
    
    # Set randomization (for global docking)
    if randomize1:
        docking.set_randomize1(True)
    if randomize2:
        docking.set_randomize2(True)
    
    # Set score functions
    scorefxn_low = rosetta.core.scoring.ScoreFunctionFactory.create_score_function("interchain_cen")
    scorefxn_high = pyrosetta.get_fa_scorefxn()
    
    docking.set_lowres_scorefxn(scorefxn_low)
    docking.set_highres_scorefxn(scorefxn_high)
    
    return docking


def run_single_docking(pose, docking_protocol=None, scorefxn=None):
    """
    Run a single docking simulation.
    
    Args:
        pose: Input Pose object (will be copied, not modified)
        docking_protocol: DockingProtocol mover. If None, creates default
        scorefxn: Score function for final scoring. If None, uses default
    
    Returns:
        Tuple of (docked_pose, score)
    """
    # Create a working copy
    work_pose = pyrosetta.Pose()
    work_pose.assign(pose)
    
    # Setup protocol if not provided
    if docking_protocol is None:
        docking_protocol = setup_docking_protocol()
    
    # Setup score function if not provided
    if scorefxn is None:
        scorefxn = pyrosetta.get_fa_scorefxn()
    
    # Run docking
    docking_protocol.apply(work_pose)
    
    # Get final score
    score = scorefxn(work_pose)
    
    return work_pose, score


def run_docking(pose, n_runs=10, save_all=True, verbose=False):
    """
    Run multiple docking simulations.
    
    This is the main docking function that runs the protocol multiple times
    and returns all results.
    
    Args:
        pose: Input complex Pose (from prepare_structures)
        n_runs: Number of independent docking runs (default: 10)
        save_all: If True, return all poses. If False, return only scores (default: True)
        verbose: If True, print progress (default: False)
    
    Returns:
        List of dictionaries with keys:
            - 'pose': Docked Pose object (if save_all=True)
            - 'score': Total score
            - 'run': Run number
    
    Example:
        >>> results = run_docking(complex_pose, n_runs=100)
        >>> best = min(results, key=lambda x: x['score'])
        >>> print(f"Best score: {best['score']:.2f}")
    """
    if verbose:
        print(f"Running {n_runs} docking simulations...")
    
    # Setup docking protocol once (reuse for efficiency)
    docking_protocol = setup_docking_protocol()
    scorefxn = pyrosetta.get_fa_scorefxn()
    
    results = []
    
    for i in range(n_runs):
        if verbose and (i + 1) % 10 == 0:
            print(f"  Completed {i + 1}/{n_runs} runs...")
        
        # Run docking
        docked_pose, score = run_single_docking(pose, docking_protocol, scorefxn)
        
        # Store results
        result = {
            'run': i + 1,
            'score': score
        }
        
        if save_all:
            result['pose'] = docked_pose
        
        results.append(result)
    
    if verbose:
        print(f"Docking complete. {n_runs} structures generated.")
    
    return results


def get_interface_score(pose, scorefxn=None):
    """
    Calculate the interface score (binding energy) for a docked complex.
    
    The interface score represents the energy change upon binding.
    
    Args:
        pose: Docked complex Pose
        scorefxn: Score function to use. If None, uses default
    
    Returns:
        Interface score (lower is better)
    """
    if scorefxn is None:
        scorefxn = pyrosetta.get_fa_scorefxn()
    
    # Create interface analyzer
    interface_analyzer = rosetta.protocols.analysis.InterfaceAnalyzerMover()
    interface_analyzer.set_scorefunction(scorefxn)
    
    # Apply to get interface metrics
    interface_analyzer.apply(pose)
    
    # Get interface score
    # This is stored in the pose's datacache
    interface_score = pose.scores['dG_separated']
    
    return interface_score


def save_docked_structure(pose, output_path):
    """
    Save a docked structure to PDB file.
    
    Args:
        pose: Docked Pose object
        output_path: Path for output PDB file
    """
    pose.dump_pdb(str(output_path))
