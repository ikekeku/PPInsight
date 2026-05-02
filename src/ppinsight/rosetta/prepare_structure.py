"""
Structure preparation module for protein docking.

This module handles:
- Loading PDB structures
- Structure relaxation/refinement
- Combining proteins with jumps for docking
"""

import sys
import tempfile
from pathlib import Path

from ppinsight.utils import copy_pdb_selected_chains, dbref_chains_for_accession

try:
    import pyrosetta
    from pyrosetta import rosetta  # pylint: disable=no-member, import-error
except ImportError:
    print("ERROR: PyRosetta not found!")
    print("Please install PyRosetta: pip install pyrosetta-*.whl")
    sys.exit(1)


_CHAIN_ID_ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"


def initialize_pyrosetta(verbose=False):
    """
    Initialize PyRosetta with appropriate options.

    Args:
        verbose: If True, show PyRosetta output. Default False.

    Returns:
        True if successful
    """
    # Initialize with flags suitable for docking.
    # -detect_disulf false: prevents RuntimeError when docking perturbation
    #   separates chains that share a disulfide bond (the scoring
    #   function cannot find the partner after rigid-body moves).
    init_flags = "-mute all -detect_disulf false -ignore_unrecognized_res true"

    if verbose:
        init_flags = "-detect_disulf false -ignore_unrecognized_res true"

    pyrosetta.init(init_flags) # pylint: disable=no-member, import-error
    return True


def _pose_chain_ids(pose):
    """Return the PDB chain IDs for each conformation chain in *pose*."""
    pdb_info = pose.pdb_info()
    if pdb_info is None:
        raise RuntimeError("Pose is missing PDB chain information")

    return [
        pdb_info.chain(pose.chain_begin(chain_index))
        for chain_index in range(1, pose.num_chains() + 1)
    ]


def _assign_unique_chain_ids(pose):
    """Assign unique one-character chain IDs across all chains in *pose*."""
    if pose.num_chains() > len(_CHAIN_ID_ALPHABET):
        raise RuntimeError(
            "Rosetta docking supports at most "
            f"{len(_CHAIN_ID_ALPHABET)} uniquely addressable chains in this "
            "wrapper"
        )

    pdb_info = pose.pdb_info()
    if pdb_info is None:
        raise RuntimeError("Pose is missing PDB chain information")

    for chain_index in range(1, pose.num_chains() + 1):
        chain_id = _CHAIN_ID_ALPHABET[chain_index - 1]
        for residue_index in range(
            pose.chain_begin(chain_index),
            pose.chain_end(chain_index) + 1,
        ):
            pdb_info.chain(residue_index, chain_id)

    pdb_info.rebuild_pdb2pose()
    return _pose_chain_ids(pose)


def _load_partner_pose(
    pdb_path,
    verbose=False,
    partner_label="Partner",
    auto_filter=True,
):
    """Load one docking partner, filtering to accession-mapped protein chains."""
    pdb_path = Path(pdb_path)
    accession = pdb_path.stem.upper()
    allowed_chains = set()
    if auto_filter:
        allowed_chains = set(dbref_chains_for_accession(pdb_path, accession))

    filtered_path = None
    load_path = pdb_path
    if allowed_chains:
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as handle:
            filtered_path = Path(handle.name)

        copy_pdb_selected_chains(pdb_path, filtered_path, allowed_chains)
        load_path = filtered_path

        if verbose:
            print(
                f"Selected {partner_label.lower()} chains from DBREF for "
                f"{accession}: {', '.join(sorted(allowed_chains))}"
            )

    try:
        pose = pyrosetta.pose_from_pdb(str(load_path))
    finally:
        if filtered_path is not None:
            filtered_path.unlink(missing_ok=True)

    # Rosetta docking expects protein-only partners unless extra params are
    # supplied for ligands or other non-canonical residues.
    rosetta.core.pose.remove_nonprotein_residues(pose)

    if pose.total_residue() == 0:
        raise RuntimeError(
            f"{partner_label} contains no protein residues after Rosetta "
            f"preprocessing: {pdb_path}"
        )

    chain_ids = _pose_chain_ids(pose)
    if verbose:
        print(
            f"{partner_label}: {pose.total_residue()} protein residues across "
            f"chains {', '.join(chain_ids)}"
        )

    return pose, chain_ids


def load_structure(pdb_path, auto_filter=True):
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

    pose, _ = _load_partner_pose(pdb_path, auto_filter=auto_filter)

    # Disulfide detection is handled by the -detect_disulf init flag.
    # No additional fix_disulfides call is needed.

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
    fast_relax = rosetta.protocols.relax.FastRelax()  # pylint: disable=no-member
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
    # Disulfide detection is handled by the -detect_disulf init flag
    # passed during pyrosetta.init(). Nothing else to fix here.

    return pose


def combine_proteins(
    pose1,
    pose2,
    jump_distance=15.0,
    partner1_chain_ids=None,
    partner2_chain_ids=None,
):
    """
    Combine two protein poses with a jump for docking.

    This preserves internal chain breaks within each partner and then lets
    Rosetta build the docking fold tree for the two partner groups.

    Args:
        pose1: First protein Pose (will be chain A)
        pose2: Second protein Pose (will be chain B)
        jump_distance: Initial separation distance in Angstroms (default: 15.0)

    Returns:
        Combined Pose object with proper fold tree for docking
    """
    fix_structure_issues(pose1)
    fix_structure_issues(pose2)

    if partner1_chain_ids is None:
        partner1_chain_ids = _pose_chain_ids(pose1)
    if partner2_chain_ids is None:
        partner2_chain_ids = _pose_chain_ids(pose2)

    combined_pose = pyrosetta.Pose()
    combined_pose.assign(pose1)

    # Keep the ligand's internal chain topology intact instead of forcing a
    # polymer bond across chain termini.
    combined_pose.append_pose_by_jump(pose2, 1)

    combined_chain_ids = _assign_unique_chain_ids(combined_pose)
    expected_chain_count = len(partner1_chain_ids) + len(partner2_chain_ids)
    if len(combined_chain_ids) != expected_chain_count:
        raise RuntimeError(
            "Combined Rosetta pose has an unexpected chain count: "
            f"expected {expected_chain_count}, found {len(combined_chain_ids)}"
        )

    partner_string = (
        "".join(combined_chain_ids[: len(partner1_chain_ids)])
        + "_"
        + "".join(combined_chain_ids[len(partner1_chain_ids) :])
    )

    movable_jumps = rosetta.utility.vector1_int()  # pylint: disable=no-member
    movable_jumps.append(1)
    rosetta.protocols.docking.setup_foldtree(  # pylint: disable=no-member
        combined_pose,
        partner_string,
        movable_jumps,
    )
    rosetta.core.pose.add_comment(  # pylint: disable=no-member
        combined_pose,
        "ppinsight_partners",
        partner_string,
    )

    # Disulfide detection already handled by -detect_disulf init flag.

    jump = combined_pose.jump(1)
    translation = rosetta.numeric.xyzVector_double_t(jump_distance, 0, 0)  # pylint: disable=no-member
    jump.set_translation(translation)
    combined_pose.set_jump(1, jump)

    return combined_pose


def prepare_structures(
    protein1_pdb, protein2_pdb,
    relax=True, jump_distance=15.0, verbose=False,
    auto_filter=True,
):
    """
    Complete structure preparation pipeline.

    Args:
        protein1_pdb: Path to first protein PDB file
        protein2_pdb: Path to second protein PDB file
        relax: If True, relax structures before combining (default: True)
        jump_distance: Initial separation distance (default: 15.0 Å)
        verbose: If True, show detailed output (default: False)

    Returns:
        Combined Pose object ready for docking
    """
    if verbose:
        print("Initializing PyRosetta...")
    initialize_pyrosetta(verbose=verbose)

    if verbose:
        print(f"Loading {protein1_pdb}...")
    pose1, pose1_chain_ids = _load_partner_pose(
        protein1_pdb,
        verbose=verbose,
        partner_label="Receptor",
        auto_filter=auto_filter,
    )

    if verbose:
        print(f"Loading {protein2_pdb}...")
    pose2, pose2_chain_ids = _load_partner_pose(
        protein2_pdb,
        verbose=verbose,
        partner_label="Ligand",
        auto_filter=auto_filter,
    )

    if verbose:
        print(f"Protein 1: {pose1.total_residue()} residues")
        print(f"Protein 2: {pose2.total_residue()} residues")

    if relax:
        if verbose:
            print("Relaxing protein 1...")
        relax_structure(pose1)

        if verbose:
            print("Relaxing protein 2...")
        relax_structure(pose2)

    if verbose:
        print(f"Combining proteins (jump distance: {jump_distance} Å)...")
    combined_pose = combine_proteins(
        pose1,
        pose2,
        jump_distance,
        partner1_chain_ids=pose1_chain_ids,
        partner2_chain_ids=pose2_chain_ids,
    )

    if verbose:
        print(f"Combined complex: {combined_pose.total_residue()} residues")

    return combined_pose
