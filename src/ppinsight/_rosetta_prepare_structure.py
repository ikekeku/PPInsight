"""
Structure preparation module for protein docking.

This module handles:
- Loading PDB structures
- Structure relaxation/refinement
- Combining proteins with jumps for docking
"""

import sys
from pathlib import Path

try:
    import pyrosetta
    from pyrosetta import rosetta
except ImportError:
    print("ERROR: PyRosetta not found!")
    print("Please install PyRosetta: pip install pyrosetta-*.whl")
    sys.exit(1)


def initialize_pyrosetta(verbose=False):
    """
    Initialize PyRosetta with appropriate options.
    
    Args:
        verbose: If True, show PyRosetta output. Default False.
    
    Returns:
        True if successful
    """
    # Initialize with flags to handle disulfides properly
    init_flags = "-mute all -detect_disulf true -ignore_unrecognized_res true"
    
    if verbose:
        init_flags = "-detect_disulf true -ignore_unrecognized_res true"
    
    pyrosetta.init(init_flags)
    return True


def load_structure(pdb_path):
    """
    Load a protein structure from PDB file.
    
    Args:
        pdb_path: Path to PDB file (str or Path)
    
    Returns:
        PyRosetta Pose object
    
    Raises:
        FileNotFoundError: If PDB file doesn't exist
    """
    pdb_path = Path(pdb_path)
    
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")
    
    pose = pyrosetta.pose_from_pdb(str(pdb_path))
    
    # Fix disulfides if present
    try:
        rosetta.core.conformation.fix_disulfides(pose)
    except:
        pass  # No disulfides or already fixed
    
    return pose


def relax_structure(pose, scorefxn=None):
    """
    Relax a protein structure using FastRelax protocol.
    
    This performs energy minimization and side-chain repacking to 
    refine the structure before docking.
    
    Args:
        pose: PyRosetta Pose object
        scorefxn: Score function to use. If None, uses default ref2015
    
    Returns:
        Relaxed Pose object (modifies in place but also returns)
    """
    if scorefxn is None:
        scorefxn = pyrosetta.get_fa_scorefxn()
    
    # Setup FastRelax protocol
    fast_relax = rosetta.protocols.relax.FastRelax()
    fast_relax.set_scorefxn(scorefxn)
    
    # Apply relaxation
    fast_relax.apply(pose)
    
    return pose


def fix_structure_issues(pose):
    """
    Fix common structure issues before docking.
    
    Args:
        pose: PyRosetta Pose object
    
    Returns:
        Fixed Pose object
    """
    # Fix disulfide bonds
    try:
        rosetta.core.conformation.fix_disulfides(pose)
    except Exception as e:
        print(f"Warning: Could not auto-fix disulfides: {e}")
    
    # Remove any waters or ligands that might cause issues
    # (Optional - only if needed)
    
    return pose


def combine_proteins(pose1, pose2, jump_distance=15.0):
    """
    Combine two protein poses with a jump for docking.
    
    This creates a complex with two chains separated by the specified distance,
    with a proper fold tree for docking.
    
    Args:
        pose1: First protein Pose (will be chain A)
        pose2: Second protein Pose (will be chain B)
        jump_distance: Initial separation distance in Angstroms (default: 15.0)
    
    Returns:
        Combined Pose object with proper fold tree for docking
    """
    # Fix any structure issues first
    fix_structure_issues(pose1)
    fix_structure_issues(pose2)
    
    # Create combined pose
    combined_pose = pyrosetta.Pose()
    combined_pose.assign(pose1)
    
    # Append second protein with new chain
    rosetta.core.pose.append_pose_to_pose(
        combined_pose, 
        pose2, 
        new_chain=True
    )
    
    # Fix disulfides in combined pose
    try:
        rosetta.core.conformation.fix_disulfides(combined_pose)
    except:
        pass
    
    # Setup fold tree for docking
    chain1_end = pose1.total_residue()
    chain2_start = chain1_end + 1
    
    ft = rosetta.core.kinematics.FoldTree()
    ft.clear()
    
    # Add edges and jump
    ft.add_edge(1, chain1_end, -1)  # Chain 1 internal
    ft.add_edge(1, chain2_start, 1)  # Jump between chains
    ft.add_edge(chain2_start, combined_pose.total_residue(), -1)  # Chain 2 internal
    
    combined_pose.fold_tree(ft)
    
    # Set initial jump distance
    jump = combined_pose.jump(1)
    translation = rosetta.numeric.xyzVector_double_t(jump_distance, 0, 0)
    jump.set_translation(translation)
    combined_pose.set_jump(1, jump)
    
    return combined_pose


def prepare_structures(protein1_pdb, protein2_pdb, relax=True, jump_distance=15.0, verbose=False):
    """
    Complete structure preparation pipeline.
    
    This is the main function that:
    1. Initializes PyRosetta
    2. Loads both protein structures
    3. Optionally relaxes them
    4. Combines them for docking
    
    Args:
        protein1_pdb: Path to first protein PDB file
        protein2_pdb: Path to second protein PDB file
        relax: If True, relax structures before combining (default: True)
        jump_distance: Initial separation distance (default: 15.0 Å)
        verbose: If True, show detailed output (default: False)
    
    Returns:
        Combined Pose object ready for docking
    
    Example:
        >>> complex_pose = prepare_structures("protein1.pdb", "protein2.pdb")
        >>> print(f"Complex has {complex_pose.total_residue()} residues")
    """
    if verbose:
        print("Initializing PyRosetta...")
    initialize_pyrosetta(verbose=verbose)
    
    # Load structures
    if verbose:
        print(f"Loading {protein1_pdb}...")
    pose1 = load_structure(protein1_pdb)
    
    if verbose:
        print(f"Loading {protein2_pdb}...")
    pose2 = load_structure(protein2_pdb)
    
    if verbose:
        print(f"Protein 1: {pose1.total_residue()} residues")
        print(f"Protein 2: {pose2.total_residue()} residues")
    
    # Relax structures
    if relax:
        if verbose:
            print("Relaxing protein 1...")
        relax_structure(pose1)
        
        if verbose:
            print("Relaxing protein 2...")
        relax_structure(pose2)
    
    # Combine proteins
    if verbose:
        print(f"Combining proteins (jump distance: {jump_distance} Å)...")
    combined_pose = combine_proteins(pose1, pose2, jump_distance)
    
    if verbose:
        print(f"Combined complex: {combined_pose.total_residue()} residues")
    
    return combined_pose