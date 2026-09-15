"""
Docking execution module — RosettaDock via PyRosetta.

Best practices from the RosettaDock protocol:

* **Prepacking** is mandatory before docking — removes internal clashes in
  each partner so that docking scores are not contaminated by intra-molecular
  artefacts. The packer task is restricted to repacking: PyRosetta's default
  task allows design at every position, which would silently replace the
  input sequence.
* ``-ex1 -ex2aro`` extra rotamer sampling must be enabled for accurate
  side-chain packing at the interface.
* **Global docking** (PPInsight default) randomises both partners'
  orientations and spins the ligand before each trajectory (Rosetta's
  ``-randomize1 -randomize2 -spin``), then runs the standard two-stage
  ``DockingProtocol``: centroid low-resolution Monte Carlo search followed by
  full-atom high-resolution refinement with interface repacking.
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

try:
    import pyrosetta
    from pyrosetta import rosetta
    from pyrosetta.rosetta.core.pack.task import TaskFactory, operation
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
except ImportError:
    print("ERROR: PyRosetta not found!")
    print("Please install PyRosetta: pip install pyrosetta-*.whl")
    sys.exit(1)


DOCKING_JUMP = 1

# Rosetta's -dock_pert defaults for local docking: 3 degrees, 8 Angstrom.
LOCAL_PERTURB_ROT_DEG = 3.0
LOCAL_PERTURB_TRANS_ANG = 8.0

# DockingProtocol's low-resolution filter rejects trajectories that end
# without interchain contact or with clashes; the Rosetta app re-runs such
# jobs (FAIL_RETRY). Cap the retries so a pathological input can't loop.
MAX_DOCKING_ATTEMPTS = 10


class DockingTrajectoryError(RuntimeError):
    """Raised when every docking attempt was rejected by Rosetta's filters."""


# ---------------------------------------------------------------------------
# Initialisation with best-practice flags
# ---------------------------------------------------------------------------

_INIT_DONE = False


def ensure_init(extra_flags: str = "", *, mute: bool = True):
    """Initialise PyRosetta with best-practice flags (idempotent).

    Always includes ``-ex1 -ex2aro`` for extra rotamer sampling at the
    interface, which is **mandatory** for accurate RosettaDock results.
    ``prepare_structure.initialize_pyrosetta`` passes the same flags, so the
    pipeline gets them whichever entry point initialises first.
    """
    global _INIT_DONE
    if _INIT_DONE:
        return

    is_initialized = getattr(pyrosetta, "is_initialized", None)
    if callable(is_initialized) and is_initialized():
        _INIT_DONE = True
        return

    init_flags = "-ex1 -ex2aro"
    if mute:
        init_flags += " -mute all"
    if extra_flags:
        init_flags += " " + extra_flags
    pyrosetta.init(init_flags)
    _INIT_DONE = True


# ---------------------------------------------------------------------------
# Pre-packing (mandatory before docking)
# ---------------------------------------------------------------------------

def repack_only_task_factory():
    """Task factory that repacks side chains without changing the sequence.

    A bare ``TaskFactory()`` yields a task where every residue is
    designable, so a packer built from it rewrites the protein's sequence.
    ``RestrictToRepacking`` turns that off; ``InitializeFromCommandline``
    picks up ``-ex1 -ex2aro``; ``IncludeCurrent`` keeps the input rotamers
    in the search.
    """
    task_factory = TaskFactory()
    task_factory.push_back(operation.InitializeFromCommandline())
    task_factory.push_back(operation.IncludeCurrent())
    task_factory.push_back(operation.RestrictToRepacking())
    return task_factory


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

    packer = PackRotamersMover()
    packer.task_factory(repack_only_task_factory())
    packer.score_function(scorefxn)
    packer.apply(pose)
    return pose


# ---------------------------------------------------------------------------
# Docking protocol
# ---------------------------------------------------------------------------

def setup_docking_protocol(global_docking=True):
    """Build the standard RosettaDock ``DockingProtocol`` for jump 1.

    Global docking runs both stages: a centroid low-resolution Monte Carlo
    rigid-body search, then full-atom high-resolution refinement
    (``DockMCMProtocol``) with interface repacking and side-chain recovery.
    Local docking (``global_docking=False``) skips the low-resolution stage
    and only refines around the starting orientation.

    The fold tree is not rebuilt here (``autofoldtree=False``):
    ``prepare_structure.combine_proteins`` already set up the docking fold
    tree for the receptor/ligand partner split.

    Calls :func:`ensure_init` to guarantee ``-ex1 -ex2aro`` flags.
    """
    ensure_init()
    protocol = rosetta.protocols.docking.DockingProtocol(  # pylint: disable=no-member
        DOCKING_JUMP,
        False,                 # low_res_protocol_only
        not global_docking,    # docking_local_refine
        False,                 # autofoldtree
    )
    movable_jumps = rosetta.utility.vector1_int()  # pylint: disable=no-member
    movable_jumps.append(DOCKING_JUMP)
    protocol.set_movable_jumps(movable_jumps)
    return protocol


def randomize_partners(pose, jump=DOCKING_JUMP):
    """Fully randomise the starting orientation for global docking.

    Equivalent to Rosetta's ``-randomize1 -randomize2 -spin``: each partner
    is rotated uniformly at random about its own centroid, then the
    downstream partner is spun about the axis joining the two centroids.
    Modifies *pose* in place.
    """
    rigid = rosetta.protocols.rigid  # pylint: disable=no-member
    rigid.RigidBodyRandomizeMover(pose, jump, rigid.partner_upstream).apply(pose)
    rigid.RigidBodyRandomizeMover(pose, jump, rigid.partner_downstream).apply(pose)
    rigid.RigidBodySpinMover(jump).apply(pose)


def run_single_docking(pose, docking_protocol=None, scorefxn=None,
                       randomize=True, global_docking=True):
    """
    Run a single docking simulation.

    For **global docking** (default), the starting orientation is fully
    randomised (see :func:`randomize_partners`); for local docking a small
    ``-dock_pert``-style perturbation is applied instead.

    Args:
        pose: Input Pose object (will be copied, not modified)
        docking_protocol: Docking mover. If None, creates default
        scorefxn: Score function for final scoring. If None, uses default
        randomize: If True, perturb the initial orientation
        global_docking: If True (default), fully randomise the orientation
            and run the two-stage protocol (no prior knowledge of the
            binding site).

    Returns:
        Tuple of (docked_pose, score)

    Raises:
        DockingTrajectoryError: if Rosetta's docking filters rejected every
            one of ``MAX_DOCKING_ATTEMPTS`` trajectories.
    """
    ensure_init()
    if scorefxn is None:
        scorefxn = pyrosetta.get_fa_scorefxn()
    if docking_protocol is None:
        docking_protocol = setup_docking_protocol(global_docking=global_docking)

    ms = rosetta.protocols.moves.MoverStatus  # pylint: disable=no-member
    for _attempt in range(MAX_DOCKING_ATTEMPTS):
        work_pose = pyrosetta.Pose()
        work_pose.assign(pose)

        if randomize:
            if global_docking:
                randomize_partners(work_pose)
            else:
                rosetta.protocols.rigid.RigidBodyPerturbMover(  # pylint: disable=no-member
                    DOCKING_JUMP, LOCAL_PERTURB_ROT_DEG, LOCAL_PERTURB_TRANS_ANG,
                ).apply(work_pose)

        docking_protocol.apply(work_pose)
        status = docking_protocol.get_last_move_status()
        if status == ms.MS_SUCCESS and work_pose.is_fullatom():
            return work_pose, scorefxn(work_pose)
        if status in (ms.FAIL_DO_NOT_RETRY, ms.FAIL_BAD_INPUT):
            raise DockingTrajectoryError(
                f"DockingProtocol refused the input pose (status {status})."
            )

    raise DockingTrajectoryError(
        f"All {MAX_DOCKING_ATTEMPTS} docking trajectories were rejected by "
        "Rosetta's low-resolution filters (no interchain contact or clashes)."
    )


def run_docking(pose, n_runs=10, save_all=True, verbose=False,
                global_docking=True, skip_prepack=False):
    """
    Run multiple docking simulations.

    This is the main docking function that runs the protocol multiple times
    and returns all results.

    **Best-practice notes** (RosettaDock protocol):

    * Pre-packing is run automatically before docking unless
      ``skip_prepack=True``.
    * For **global docking** (``global_docking=True``, the default for
      PPInsight), each run fully randomises the starting orientation and
      runs the two-stage low-/high-resolution protocol — no prior knowledge
      of the binding site is assumed.
    * Production global docking requires **10 000–100 000 decoys**
      (``n_runs``).  A warning is emitted when fewer than 1000 are
      requested so users are aware of the sampling implications.

    Args:
        pose: Input complex Pose (from prepare_structures)
        n_runs: Number of independent docking runs (default: 10).
            For publication-quality global docking, use 10 000–100 000.
        save_all: If True, return all poses. If False, return only scores
        verbose: If True, print progress
        global_docking: If True (default), fully randomize orientation for
            each run (no prior binding-site knowledge assumed).
        skip_prepack: If True, skip the mandatory pre-packing step.

    Returns:
        List of dictionaries with keys:
            - 'pose': Docked Pose object (if save_all=True)
            - 'score': Legacy alias of the Rosetta total score
            - 'total_score': Rosetta total score
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
    docking_protocol = setup_docking_protocol(global_docking=global_docking)
    scorefxn = pyrosetta.get_fa_scorefxn()

    results = []

    for i in range(n_runs):
        if verbose and (i + 1) % 10 == 0:
            print(f"  Completed {i + 1}/{n_runs} runs...")

        docked_pose, total_score = run_single_docking(
            pose,
            docking_protocol,
            scorefxn,
            randomize=True,
            global_docking=global_docking,
        )

        # Extract interface score (the primary RosettaDock metric)
        i_sc = get_interface_score(docked_pose, scorefxn)

        result = {
            'run': i + 1,
            'description': f"decoy_{i + 1}",
            'score': total_score,
            'total_score': total_score,
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
