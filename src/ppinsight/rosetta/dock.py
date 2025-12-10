"""
Docking execution module.

This module handles:
- Setting up docking protocols
- Running multiple docking simulations
- Managing docking movers and score functions
"""

import sys
# pylint: disable=no-member, import-error
from pyrosetta.rosetta.protocols.moves import SequenceMover
from pyrosetta.rosetta.protocols.docking import FaDockingSlideIntoContact
from pyrosetta.rosetta.protocols.minimization_packing import MinMover, PackRotamersMover
from pyrosetta.rosetta.core.pack.task import TaskFactory
# pylint: enable=no-member, import-error

try:
    import pyrosetta
    from pyrosetta import rosetta
except ImportError:
    print("ERROR: PyRosetta not found!")
    print("Please install PyRosetta: pip install pyrosetta-*.whl")
    sys.exit(1)


def setup_docking_protocol():
    """
    Setup a complete docking protocol using DockMCMProtocol.

    Returns:
        Configured docking protocol
    """
    # Use the full DockMCMProtocol which includes low-res and high-res docking
    docking = rosetta.protocols.docking.DockMCMProtocol()  # pylint: disable=no-member

    # Set score functions
    scorefxn = pyrosetta.get_fa_scorefxn()
    docking.set_scorefxn(scorefxn)

    return docking


def setup_simple_docking():
    """
    Setup a simple but complete docking workflow.

    Returns:
        Configured SequenceMover with full docking pipeline
    """
    # Create a sequence of movers for docking

    # Score function
    scorefxn = pyrosetta.get_fa_scorefxn()

    # Create movemap for docking
    movemap = rosetta.core.kinematics.MoveMap()
    movemap.set_jump(1, True)  # Allow rigid body movement
    movemap.set_bb(False)      # Don't move backbone
    movemap.set_chi(True)      # Allow side chain movement

    # Create movers
    # 1. Slide proteins into contact
    slide_into_contact = FaDockingSlideIntoContact(1)  # jump number = 1

    # 2. Minimize
    min_mover = MinMover()
    min_mover.movemap(movemap)
    min_mover.score_function(scorefxn)

    # Combine into sequence
    sequence = SequenceMover()
    sequence.add_mover(slide_into_contact)
    sequence.add_mover(min_mover)

    return sequence


def run_single_docking(pose, docking_protocol=None, scorefxn=None, randomize=True):
    """
    Run a single docking simulation.

    Args:
        pose: Input Pose object (will be copied, not modified)
        docking_protocol: Docking mover. If None, creates default
        scorefxn: Score function for final scoring. If None, uses default
        randomize: If True, randomize initial orientation

    Returns:
        Tuple of (docked_pose, score)
    """
    # Create a working copy
    work_pose = pyrosetta.Pose()
    work_pose.assign(pose)

    # Setup score function if not provided
    if scorefxn is None:
        scorefxn = pyrosetta.get_fa_scorefxn()

    # Randomize initial position if requested
    if randomize:
        # Randomize both translation and rotation
        rigid_body_perturb = rosetta.protocols.rigid.RigidBodyPerturbMover(1, 8.0, 8.0)
        rigid_body_perturb.apply(work_pose)

    # Setup protocol if not provided
    if docking_protocol is None:
        docking_protocol = setup_simple_docking()

    # Run docking
    try:
        docking_protocol.apply(work_pose)
    except RuntimeError as e:
        print(f"Docking error: {e}")


    # Get final score
    score = scorefxn(work_pose)

    return work_pose, score


def setup_full_docking_protocol():
    """
    Setup a complete docking protocol with all stages.

    Returns:
        Configured SequenceMover with complete docking pipeline
    """
    # from pyrosetta.rosetta.protocols.moves import SequenceMover
    # from pyrosetta.rosetta.protocols.docking import FaDockingSlideIntoContact
    # from pyrosetta.rosetta.protocols.minimization_packing import MinMover, PackRotamersMover
    # from pyrosetta.rosetta.core.pack.task import TaskFactory

    # scorefxn = pyrosetta.get_fa_scorefxn()
    scorefxn = pyrosetta.create_score_function("ref2015")
    # Reduce disulfide weight
    scorefxn.set_weight(rosetta.core.scoring.ScoreType.dslf_fa13, 0.5)

    # Create movemap
    movemap = rosetta.core.kinematics.MoveMap()
    movemap.set_jump(1, True)
    movemap.set_bb(False)
    movemap.set_chi(True)

    # Stage 1: Randomize orientation
    rigid_body_perturb = rosetta.protocols.rigid.RigidBodyPerturbMover(1, 8.0, 8.0)

    # Stage 2: Slide into contact
    slide_into_contact = FaDockingSlideIntoContact(1)

    # Stage 3: Pack rotamers at interface
    task_factory = TaskFactory()
    pack_mover = PackRotamersMover()
    pack_mover.task_factory(task_factory)
    pack_mover.score_function(scorefxn)

    # Stage 4: Minimize
    min_mover = MinMover()
    min_mover.movemap(movemap)
    min_mover.score_function(scorefxn)

    # Combine into sequence
    sequence = SequenceMover()
    sequence.add_mover(rigid_body_perturb)
    sequence.add_mover(slide_into_contact)
    sequence.add_mover(pack_mover)
    sequence.add_mover(min_mover)

    return sequence


def run_docking(pose, n_runs=10, save_all=True, verbose=False, use_full_protocol=False):
    """
    Run multiple docking simulations.

    This is the main docking function that runs the protocol multiple times
    and returns all results.

    Args:
        pose: Input complex Pose (from prepare_structures)
        n_runs: Number of independent docking runs (default: 10)
        save_all: If True, return all poses. If False, return only scores (default: True)
        verbose: If True, print progress (default: False)
        use_full_protocol: If True, use full protocol with packing (slower but better)

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
    if use_full_protocol:
        docking_protocol = setup_full_docking_protocol()
    else:
        docking_protocol = setup_simple_docking()

    scorefxn = pyrosetta.get_fa_scorefxn()

    results = []

    for i in range(n_runs):
        if verbose and (i + 1) % 10 == 0:
            print(f"  Completed {i + 1}/{n_runs} runs...")

        # Run docking with randomization
        docked_pose, score = run_single_docking(
            pose,
            docking_protocol,
            scorefxn,
            randomize=True  # Always randomize for global docking
        )

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

    # Get interface score from pose datacache
    try:
        interface_score = pose.scores['dG_separated']
    except KeyError:
        # Fallback: calculate manually
        interface_score = scorefxn(pose)

    return interface_score


def save_docked_structure(pose, output_path):
    """
    Save a docked structure to PDB file.

    Args:
        pose: Docked Pose object
        output_path: Path for output PDB file
    """
    pose.dump_pdb(str(output_path))
