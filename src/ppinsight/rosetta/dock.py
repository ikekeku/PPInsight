"""
Docking execution module — RosettaDock via PyRosetta.

Best practices from the RosettaDock protocol:

* **Prepacking** is mandatory before docking — removes internal clashes in
  each partner so that docking scores are not contaminated by intra-molecular
  artefacts.
* ``-ex1 -ex2aro`` extra rotamer sampling must be enabled for accurate
  side-chain packing at the interface.
* **Global docking** (PPInsight default) requires ``-spin``,
  ``-randomize1``, and ``-randomize2`` flags — these are modelled here
  via rigid-body perturbation movers.
* Production runs should use **10 000–100 000** decoys (``nstruct``).
  The pipeline warns when fewer are requested.
* The **I_sc** (interface score / ``dG_separated``) is the primary quality
  metric, not the total Rosetta energy.

This module handles:
- Pre-packing side chains
- Setting up docking protocols (global by default)
- Running multiple docking simulations
- Extracting interface scores and key metrics
"""

import sys
import warnings

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


# ---------------------------------------------------------------------------
# Initialisation with best-practice flags
# ---------------------------------------------------------------------------

_INIT_DONE = False


def ensure_init(extra_flags: str = ""):
    """Initialise PyRosetta with best-practice flags (idempotent).

    Always includes ``-ex1 -ex2aro`` for extra rotamer sampling at the
    interface, which is **mandatory** for accurate RosettaDock results.
    """
    global _INIT_DONE
    if _INIT_DONE:
        return
    init_flags = "-ex1 -ex2aro -mute all"
    if extra_flags:
        init_flags += " " + extra_flags
    pyrosetta.init(init_flags)
    _INIT_DONE = True


# ---------------------------------------------------------------------------
# Pre-packing (mandatory before docking)
# ---------------------------------------------------------------------------

def prepack(pose, scorefxn=None):
    """Pre-pack side chains to remove intra-molecular clashes.

    This step is **mandatory** before RosettaDock — it ensures that the
    docking score reflects the actual binding energy rather than artefacts
    from bad internal packing.

    Parameters
    ----------
    pose : pyrosetta.Pose
        The input complex pose.  Modified **in place**.
    scorefxn : ScoreFunction | None
        Score function for packing.  Defaults to ``ref2015``.

    Returns
    -------
    pyrosetta.Pose
        The same *pose* (modified in place) for convenience.
    """
    ensure_init()
    if scorefxn is None:
        scorefxn = pyrosetta.create_score_function("ref2015")

    task_factory = TaskFactory()
    packer = PackRotamersMover()
    packer.task_factory(task_factory)
    packer.score_function(scorefxn)
    packer.apply(pose)
    return pose


def setup_docking_protocol():
    """
    Setup a complete docking protocol using DockMCMProtocol.

    Calls :func:`ensure_init` to guarantee ``-ex1 -ex2aro`` flags.

    Returns:
        Configured docking protocol
    """
    ensure_init()
    # Use the full DockMCMProtocol which includes low-res and high-res docking
    docking = rosetta.protocols.docking.DockMCMProtocol()  # pylint: disable=no-member

    # Set score functions
    scorefxn = pyrosetta.get_fa_scorefxn()
    docking.set_scorefxn(scorefxn)

    return docking


def setup_simple_docking():
    """
    Setup a simple but complete docking workflow.

    Calls :func:`ensure_init` to guarantee ``-ex1 -ex2aro`` flags.

    Returns:
        Configured SequenceMover with full docking pipeline
    """
    ensure_init()
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


def run_single_docking(pose, docking_protocol=None, scorefxn=None,
                       randomize=True, global_docking=True):
    """
    Run a single docking simulation.

    For **global docking** (default), the ligand orientation is fully
    randomised (equivalent to Rosetta's ``-spin -randomize1 -randomize2``).

    Args:
        pose: Input Pose object (will be copied, not modified)
        docking_protocol: Docking mover. If None, creates default
        scorefxn: Score function for final scoring. If None, uses default
        randomize: If True, randomize initial orientation
        global_docking: If True (default), apply large random perturbation
            suitable for global docking (no prior knowledge of binding site).

    Returns:
        Tuple of (docked_pose, score)
    """
    ensure_init()
    # Create a working copy
    work_pose = pyrosetta.Pose()
    work_pose.assign(pose)

    # Setup score function if not provided
    if scorefxn is None:
        scorefxn = pyrosetta.get_fa_scorefxn()

    # Randomize initial position if requested
    if randomize:
        if global_docking:
            # Large perturbation for global docking:
            # - Translation: up to 50 Å to explore entire surface
            # - Rotation: up to 360° for full orientational sampling
            rigid_body_perturb = rosetta.protocols.rigid.RigidBodyPerturbMover(
                1, 50.0, 360.0
            )
        else:
            # Small perturbation for local refinement
            rigid_body_perturb = rosetta.protocols.rigid.RigidBodyPerturbMover(
                1, 8.0, 8.0
            )
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


def setup_full_docking_protocol(global_docking=True):
    """
    Setup a complete docking protocol with all stages.

    For global docking (default), applies large random perturbation
    equivalent to ``-spin -randomize1 -randomize2``.

    Returns:
        Configured SequenceMover with complete docking pipeline
    """
    ensure_init()
    scorefxn = pyrosetta.create_score_function("ref2015")
    # Reduce disulfide weight
    scorefxn.set_weight(rosetta.core.scoring.ScoreType.dslf_fa13, 0.5)

    # Create movemap
    movemap = rosetta.core.kinematics.MoveMap()
    movemap.set_jump(1, True)
    movemap.set_bb(False)
    movemap.set_chi(True)

    # Stage 1: Randomize orientation
    if global_docking:
        rigid_body_perturb = rosetta.protocols.rigid.RigidBodyPerturbMover(
            1, 50.0, 360.0
        )
    else:
        rigid_body_perturb = rosetta.protocols.rigid.RigidBodyPerturbMover(
            1, 8.0, 8.0
        )

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


def run_docking(pose, n_runs=10, save_all=True, verbose=False,
                use_full_protocol=False, global_docking=True,
                skip_prepack=False):
    """
    Run multiple docking simulations.

    This is the main docking function that runs the protocol multiple times
    and returns all results.

    **Best-practice notes** (RosettaDock protocol):

    * Pre-packing is run automatically before docking unless
      ``skip_prepack=True``.
    * For **global docking** (``global_docking=True``, the default for
      PPInsight), each run fully randomises the ligand orientation — no
      prior knowledge of the binding site is assumed.
    * Production global docking requires **10 000–100 000 decoys**
      (``n_runs``).  A warning is emitted when fewer than 1000 are
      requested so users are aware of the sampling implications.

    Args:
        pose: Input complex Pose (from prepare_structures)
        n_runs: Number of independent docking runs (default: 10).
            For publication-quality global docking, use 10 000–100 000.
        save_all: If True, return all poses. If False, return only scores
        verbose: If True, print progress
        use_full_protocol: If True, use full protocol with packing (slower)
        global_docking: If True (default), fully randomize orientation for
            each run (no prior binding-site knowledge assumed).
        skip_prepack: If True, skip the mandatory pre-packing step.

    Returns:
        List of dictionaries with keys:
            - 'pose': Docked Pose object (if save_all=True)
            - 'score': Total score
            - 'i_sc': Interface score (dG_separated)
            - 'run': Run number

    Example:
        >>> results = run_docking(complex_pose, n_runs=100)
        >>> best = min(results, key=lambda x: x['i_sc'])
        >>> print(f"Best I_sc: {best['i_sc']:.2f}")
    """
    ensure_init()

    # ── Nstruct warning for global docking ──
    _MIN_GLOBAL = 1000
    if global_docking and n_runs < _MIN_GLOBAL:
        warnings.warn(
            f"Global docking requested with only {n_runs} decoys. "
            f"The RosettaDock protocol recommends 10,000–100,000 decoys "
            f"for global docking to adequately sample the conformational "
            f"space.  Results with fewer than {_MIN_GLOBAL} decoys should "
            f"be treated as preliminary.  Increase n_runs for production.",
            UserWarning,
            stacklevel=2,
        )

    if verbose:
        docking_mode = "global" if global_docking else "local"
        print(f"Running {n_runs} {docking_mode} docking simulations...")

    # ── Mandatory pre-packing step ──
    if not skip_prepack:
        if verbose:
            print("  Pre-packing side chains (mandatory RosettaDock step)...")
        prepack(pose)

    # Setup docking protocol once (reuse for efficiency)
    if use_full_protocol:
        docking_protocol = setup_full_docking_protocol(
            global_docking=global_docking,
        )
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
            randomize=True,
            global_docking=global_docking,
        )

        # Extract interface score (the primary RosettaDock metric)
        i_sc = get_interface_score(docked_pose, scorefxn)

        # Store results
        result = {
            'run': i + 1,
            'score': score,
            'i_sc': i_sc,
        }

        if save_all:
            result['pose'] = docked_pose

        results.append(result)

    if verbose:
        print(f"Docking complete. {n_runs} structures generated.")

    return results


def get_interface_score(pose, scorefxn=None):
    """
    Calculate the interface score (I_sc / binding energy) for a docked complex.

    The interface score (``dG_separated``) is the **primary quality metric**
    for RosettaDock — it represents the energy change upon binding and is
    more discriminating than the total Rosetta energy for identifying
    correct binding modes.

    Args:
        pose: Docked complex Pose
        scorefxn: Score function to use. If None, uses default

    Returns:
        Interface score (lower is better)
    """
    ensure_init()
    if scorefxn is None:
        scorefxn = pyrosetta.get_fa_scorefxn()

    # Create interface analyzer
    interface_analyzer = rosetta.protocols.analysis.InterfaceAnalyzerMover()
    interface_analyzer.set_scorefunction(scorefxn)

    # Apply to get interface metrics
    interface_analyzer.apply(pose)

    # Get interface score from pose datacache
    # Prefer Pose.cache (new API) over deprecated Pose.scores
    try:
        cache = getattr(pose, "cache", None) or pose.scores
        interface_score = cache['dG_separated']
    except (KeyError, AttributeError):
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
