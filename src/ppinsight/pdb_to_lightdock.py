"""
pdb_to_lightdock.py
-------------------
Runs the **full** LightDock best-practice protocol for a receptor-ligand pair:

    1. ``lightdock3_setup.py`` — prepare structures, compute swarm positions
    2. ``lightdock3.py``       — run the GSO simulation
    3. ``lgd_generate_conformations.py`` — generate PDB models (all swarms)
    4. ``lgd_cluster_bsas.py`` — cluster models per swarm (BSAS, 4 Å cutoff)
    5. ``lgd_rank.py``         — global ranking across all swarms

Steps 3-5 run by default (the official LightDock best practice).  Pass
``--skip-postprocess`` to stop after the simulation (step 2).

ANM (Anisotropic Network Model) flexibility is enabled by default (``-anm``).
Pass ``--no-anm`` to disable it for rigid docking.

The scoring function defaults to DFIRE but can be changed via ``--scoring``.

Input PDBs: passed on the command line, e.g.::

    python pdb_to_lightdock.py 2UUY_rec.pdb 2UUY_lig.pdb

Output folders::

    data/output/lightdock_runs/<receptor>_vs_<ligand>/
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
import warnings

from ppinsight.utils import (  # noqa: F401 — re-exported
    _project_root,
    resolve_input_path,
)


def run_command(cmd, cwd=None):
    """
    Print and execute a shell command.
    - cmd: list of strings, e.g. ["lightdock3.py", "setup.json", "100"]
    - cwd: directory in which to run the command
    """
    print(">>", " ".join(cmd))
    try:
        subprocess.run(cmd, cwd=cwd, check=True)
    except FileNotFoundError:
        exe = cmd[0]
        print(
            f"\nERROR: '{exe}' not found on $PATH.",
            file=sys.stderr,
        )
        print(
            "Hint: LightDock CLI tools (lightdock3_setup.py, lightdock3.py, "
            "lgd_generate_conformations.py, lgd_cluster_bsas.py, lgd_rank.py) "
            "must be installed and on your PATH.\n"
            "Install: pip install lightdock  (included with ppinsight)\n"
            "Verify:  which lightdock3_setup.py",
            file=sys.stderr,
        )
        sys.exit(2)


def make_output_dir(receptor_pdb, ligand_pdb,
                    base_root="data/output",
                    method="lightdock_runs"):
    """
    Create an organized output directory structure.

    Final structure::

        data/output/
            lightdock_runs/
                <receptor_name>_vs_<ligand_name>/

    If a run folder already exists, a numeric suffix (``_1``, ``_2``, …)
    is appended to avoid clobbering previous results.
    """
    if not os.path.isabs(base_root):
        base_root = os.path.join(_project_root(), base_root)
    os.makedirs(base_root, exist_ok=True)

    method_dir = os.path.join(base_root, method)
    os.makedirs(method_dir, exist_ok=True)

    rec_name = os.path.splitext(os.path.basename(receptor_pdb))[0]
    lig_name = os.path.splitext(os.path.basename(ligand_pdb))[0]
    run_folder_name = f"{rec_name}_vs_{lig_name}"
    run_dir = os.path.join(method_dir, run_folder_name)
    if os.path.exists(run_dir):
        idx = 1
        while True:
            candidate = os.path.join(method_dir, f"{run_folder_name}_{idx}")
            if not os.path.exists(candidate):
                run_dir = candidate
                break
            idx += 1
    os.makedirs(run_dir, exist_ok=True)
    return run_dir


# ---------------------------------------------------------------------------
# Stage 1 — setup
# ---------------------------------------------------------------------------

def _run_lightdock_setup(working_dir, rec_basename, lig_basename,
                         swarms, glowworms, *, anm=True):
    """Run the LightDock setup step (creates setup.json and initial files).

    Parameters
    ----------
    anm : bool
        Enable ANM flexibility (default True — LightDock best practice).
    """
    cmd = [
        "lightdock3_setup.py",
        rec_basename,
        lig_basename,
        "--noxt",
        "--noh",
        "--now",
    ]
    if anm:
        cmd.append("-anm")
    if swarms:
        cmd += ["-s", str(swarms)]
    if glowworms:
        cmd += ["-g", str(glowworms)]
    run_command(cmd, cwd=working_dir)


# ---------------------------------------------------------------------------
# Stage 2 — simulation
# ---------------------------------------------------------------------------

def _run_lightdock_simulation(working_dir, steps, cores, swarm_list, scoring=None):
    """Run the LightDock simulation step (lightdock3.py).

    Parameters
    ----------
    scoring : str | None
        Scoring function name (e.g. ``"dfire"``, ``"fastdfire"``,
        ``"cpydock"``).  ``None`` uses LightDock's default (fastdfire).
    """
    cmd = ["lightdock3.py", "setup.json", str(steps), "-c", str(cores)]
    if scoring:
        cmd += ["-s", scoring]
    if swarm_list:
        cmd += ["-l"] + list(map(str, swarm_list))
    run_command(cmd, cwd=working_dir)


# ---------------------------------------------------------------------------
# Stage 3 — generate PDB conformations (ALL swarms)
# ---------------------------------------------------------------------------

def _run_lightdock_generation(working_dir, rec_basename, lig_basename,
                              steps, glowworms):
    """Generate PDB models from GSO results in **every** swarm directory.

    The official LightDock protocol requires PDB models to be generated
    before clustering can run.  This iterates over all ``swarm_*`` dirs.
    """
    num_models = glowworms if glowworms else 200
    swarm_dirs = sorted(glob.glob(os.path.join(working_dir, "swarm_*")))
    if not swarm_dirs:
        warnings.warn(
            "No swarm_* directories found — skipping generation.",
            stacklevel=2,
        )
        return
    for swarm_path in swarm_dirs:
        swarm_name = os.path.basename(swarm_path)
        out_file = os.path.join(swarm_name, f"gso_{steps}.out")
        if not os.path.isfile(os.path.join(working_dir, out_file)):
            continue  # skip swarms that weren't simulated
        cmd = [
            "lgd_generate_conformations.py",
            rec_basename,
            lig_basename,
            out_file,
            str(num_models),
        ]
        run_command(cmd, cwd=working_dir)


# ---------------------------------------------------------------------------
# Stage 4 — BSAS clustering (per swarm)
# ---------------------------------------------------------------------------

def _run_lightdock_clustering(working_dir, steps):
    """Run BSAS clustering in every swarm directory.

    Executes ``lgd_cluster_bsas.py gso_<steps>.out`` inside each swarm
    directory.  This produces ``cluster.repr`` files used by the ranking
    stage and score parser.
    """
    swarm_dirs = sorted(glob.glob(os.path.join(working_dir, "swarm_*")))
    if not swarm_dirs:
        warnings.warn(
            "No swarm_* directories found — skipping clustering.",
            stacklevel=2,
        )
        return
    for swarm_path in swarm_dirs:
        gso_file = os.path.join(swarm_path, f"gso_{steps}.out")
        if not os.path.isfile(gso_file):
            continue
        cmd = ["lgd_cluster_bsas.py", f"gso_{steps}.out"]
        run_command(cmd, cwd=swarm_path)


# ---------------------------------------------------------------------------
# Stage 5 — global ranking across all swarms
# ---------------------------------------------------------------------------

def _run_lightdock_ranking(working_dir, steps):
    """Produce a global ranking of cluster representatives across swarms.

    Runs ``lgd_rank.py`` which reads all ``cluster.repr`` files and writes
    a ``rank_by_scoring.list`` (or ``rank_by_luciferin.list``) at the
    simulation root.

    ``lgd_rank.py`` requires two positional arguments: ``num_swarms`` and
    ``steps``.
    """
    num_swarms = len(glob.glob(os.path.join(working_dir, "swarm_*")))
    if num_swarms == 0:
        warnings.warn("No swarm_* directories found — skipping ranking.", stacklevel=2)
        return
    cmd = ["lgd_rank.py", str(num_swarms), str(steps), "-c", "1"]
    run_command(cmd, cwd=working_dir)


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def _execute_lightdock_stages(receptor_pdb, ligand_pdb, working_dir, opts):
    """Execute the full LightDock pipeline.

    By default all five stages (setup → simulate → generate → cluster →
    rank) are executed.  Pass ``skip_postprocess=True`` in *opts* to stop
    after the simulation (step 2).
    """
    rec_base = os.path.basename(receptor_pdb)
    lig_base = os.path.basename(ligand_pdb)
    steps = opts.get("steps", 100)

    # Stage 1 — setup
    _run_lightdock_setup(
        working_dir, rec_base, lig_base,
        opts.get("swarms"), opts.get("glowworms"),
        anm=opts.get("anm", True),
    )

    # Stage 2 — simulation
    _run_lightdock_simulation(
        working_dir, steps,
        opts.get("cores", 1),
        opts.get("swarm_list"),
        scoring=opts.get("scoring"),
    )

    # Stages 3–5: post-processing (generate + cluster + rank)
    if opts.get("skip_postprocess", False):
        print("\nSkipping post-processing (generate/cluster/rank).")
        print("Pass without --skip-postprocess for the full protocol.\n")
    else:
        # Stage 3 — generate PDB models in ALL swarms
        _run_lightdock_generation(
            working_dir, rec_base, lig_base, steps, opts.get("glowworms"),
        )
        # Stage 4 — BSAS clustering per swarm
        _run_lightdock_clustering(working_dir, steps)
        # Stage 5 — global ranking
        _run_lightdock_ranking(working_dir, steps)

    print("\nDone.")
    print("Results saved in:", working_dir)


def lightdock_pipeline(receptor_pdb, ligand_pdb, working_dir, opts=None, **kwargs):
    """Run the full LightDock workflow using an opts dict or kwargs.

    Supported keys
    ~~~~~~~~~~~~~~
    swarms, glowworms, steps, swarm_list, cores, anm (bool),
    scoring (str), skip_postprocess (bool)

    Backward-compatible: ``generate_models=True`` is accepted but ignored
    (generation is now always on unless ``skip_postprocess`` is set).
    """
    opts = dict(opts or {})
    opts.update(kwargs)

    # Backward compat: old callers may pass generate_models=True/False.
    # Map False → skip_postprocess=True, True → leave as default.
    if "generate_models" in opts and not opts["generate_models"]:
        opts.setdefault("skip_postprocess", True)

    # Copy input PDBs into the working directory (preserve filename)
    shutil.copy(receptor_pdb, working_dir)
    shutil.copy(ligand_pdb, working_dir)

    _execute_lightdock_stages(receptor_pdb, ligand_pdb, working_dir, opts)


def _cleanup_previous_outputs(workdir):
    """Remove previous LightDock outputs from a run directory.

    Returns a tuple (file_count, folder_count) with the number of removed items.
    """
    items = glob.glob(os.path.join(workdir, "lightdock*"))
    items += glob.glob(os.path.join(workdir, "swarm_*"))
    items += [os.path.join(workdir, "setup.json"), os.path.join(workdir, "init")]

    file_count = folder_count = 0
    for p in items:
        if os.path.exists(p):
            if os.path.isdir(p):
                shutil.rmtree(p, ignore_errors=True)
                folder_count += 1
            elif os.path.isfile(p):
                os.remove(p)
                file_count += 1
    return file_count, folder_count


def main(argv=None):
    """CLI entrypoint for the LightDock helper script.

    Accepts receptor and ligand PDB paths and optional LightDock parameters.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Run the full LightDock global-docking protocol: "
            "setup → simulate → generate → cluster → rank."
        ),
    )
    parser.add_argument("receptor", help="Path to receptor PDB file")
    parser.add_argument("ligand", help="Path to ligand PDB file")
    parser.add_argument(
        "--swarms",
        type=int,
        default=None,
        help=(
            "Number of swarms — independent search starting points spread "
            "over the protein surface.  More swarms explore more of the "
            "surface but increase runtime linearly.  Use the default for "
            "standard-size proteins; raise to 400+ for large complexes "
            "where the binding site is unknown.  (LightDock default if not set.)"
        ),
    )
    parser.add_argument(
        "--glowworms",
        type=int,
        default=None,
        help=(
            "Number of glowworms per swarm — agents that optimise poses "
            "within each swarm.  More glowworms give finer local sampling "
            "but increase per-swarm cost.  The default (200) works well for "
            "most cases; raise for very flexible interfaces.  "
            "(LightDock default if not set.)"
        ),
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=100,
        help=(
            "Number of GSO optimisation steps per swarm (default: 100).  "
            "More steps let glowworms converge further.  100 is the "
            "LightDock tutorial recommendation; 50 may suffice for quick "
            "screening, while 200+ is useful for difficult targets."
        ),
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help=(
            "Number of CPU cores for the simulation step (default: 1).  "
            "Each swarm runs independently, so this scales well.  Set to "
            "the number of physical cores on your machine for fastest runs."
        ),
    )
    parser.add_argument(
        "--no-anm",
        action="store_true",
        help=(
            "Disable ANM (Anisotropic Network Model) backbone flexibility.  "
            "By default ANM is enabled (LightDock best practice).  "
            "Use --no-anm only when docking rigid bodies (e.g. antibody–"
            "antigen where backbone movement is minimal) or for faster "
            "debugging runs.  Disabling ANM reduces per-step cost but "
            "can miss correct poses that require backbone rearrangement."
        ),
    )
    parser.add_argument(
        "--scoring",
        type=str,
        default=None,
        help=(
            "LightDock scoring function (e.g. 'dfire', 'fastdfire', "
            "'cpydock').  Default: DFIRE (knowledge-based, general-purpose).  "
            "Use 'cpydock' for antibody–antigen or when electrostatics "
            "dominate; use 'fastdfire' for quicker screening at slight "
            "accuracy cost.  The pipeline auto-detects the score direction "
            "for whichever function you choose."
        ),
    )
    parser.add_argument(
        "--skip-postprocess",
        action="store_true",
        help=(
            "Skip post-processing (model generation, BSAS clustering, and "
            "global ranking).  Use this when you only need raw GSO output "
            "for custom analysis, or for debugging the simulation step in "
            "isolation.  Without post-processing, 'collect' cannot parse "
            "cluster-level results."
        ),
    )
    parser.add_argument(
        "--input-dir",
        default=None,
        help=(
            "Directory to search when receptor/ligand are basenames instead "
            "of full paths (default: repo root).  Useful when PDB files "
            "live in a shared directory outside the project tree."
        ),
    )

    args = parser.parse_args(argv)

    receptor = args.receptor
    ligand = args.ligand
    try:
        receptor = resolve_input_path(receptor, search_root=args.input_dir)
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        print("Hint: use --input-dir to specify the directory containing "
              "your PDB files.", file=sys.stderr)
        sys.exit(2)
    try:
        ligand = resolve_input_path(ligand, search_root=args.input_dir)
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        print("Hint: use --input-dir to specify the directory containing "
              "your PDB files.", file=sys.stderr)
        sys.exit(2)

    workdir = make_output_dir(receptor, ligand, method="lightdock_runs")

    print(f"\nCleaning previous LightDock outputs in: {workdir} (if any)")
    file_count, folder_count = _cleanup_previous_outputs(workdir)
    print(f"\nRemoved {file_count} files and {folder_count} folders.")

    print(
        "\nRunning LightDock global-docking pipeline for:\n"
        f" Receptor: {receptor}\n"
        f" Ligand:   {ligand}\n"
        f" Output dir: {workdir}\n",
    )

    opts = {
        "swarms": args.swarms,
        "glowworms": args.glowworms,
        "steps": args.steps,
        "swarm_list": None,
        "cores": args.cores,
        "anm": not args.no_anm,
        "scoring": args.scoring,
        "skip_postprocess": args.skip_postprocess,
    }

    lightdock_pipeline(receptor, ligand, working_dir=workdir, opts=opts)

    summary_lines = [
        "LightDock global docking completed.",
        f"Run folder: {workdir}",
        "Parameters used:",
        "- Swarms: " + (str(opts["swarms"]) if opts["swarms"] else "default"),
        "- Glowworms: "
        + (str(opts["glowworms"]) if opts["glowworms"] else "default"),
        "- Steps: " + str(opts["steps"]),
        "- CPU cores: " + str(opts["cores"]),
        "- ANM flexibility: " + ("enabled" if opts["anm"] else "disabled"),
        "- Scoring function: " + (opts["scoring"] or "default (DFIRE)"),
        "- Post-processing: "
        + ("skipped" if opts["skip_postprocess"]
           else "full (generate+cluster+rank)"),
        f"Results directory: {workdir}",
    ]
    print("\n".join(summary_lines))


if __name__ == '__main__':
    main()
