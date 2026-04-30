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
import re
import shutil
import subprocess
import sys
import warnings

from ppinsight.utils import _project_root, resolve_input_path  # noqa: F401

_LIGHTDOCK_ALLOWED_PDB_RECORDS = {
    "ATOM",
    "MODEL",
    "TER",
    "ENDMDL",
    "END",
}
_LIGHTDOCK_UNSUPPORTED_RESIDUE_RE = re.compile(
    r"\[NotSupportedInScoringError\]\s+Residue\s+(?P<residue>\S+)\s+"
    r"or atom\s+(?P<atom>\S+)\s+not supported"
)


class LightDockSimulationError(RuntimeError):
    """Raised when LightDock simulation fails after setup completed."""

    def __init__(
        self,
        message: str,
        *,
        cleanable: bool = False,
        unsupported_residue: str | None = None,
        simulation_output: str = "",
    ):
        super().__init__(message)
        self.cleanable = cleanable
        self.unsupported_residue = unsupported_residue
        self.simulation_output = simulation_output


def clean_pdb_for_lightdock(
    input_pdb: str,
    output_pdb: str | None = None,
) -> dict[str, int | str]:
    """Write a protein-only PDB copy that avoids common LightDock failures.

    LightDock's DFIRE scoring accepts standard amino-acid coordinates but can
    fail on hetero records such as sulfate ions. This helper keeps only core
    coordinate records needed for protein-only docking.
    """
    if output_pdb is None:
        stem, ext = os.path.splitext(input_pdb)
        output_pdb = f"{stem}_lightdock_clean{ext or '.pdb'}"

    kept_records = 0
    removed_hetatm = 0
    removed_other = 0
    saw_end = False

    with open(input_pdb, encoding="utf-8") as src, open(
        output_pdb, "w", encoding="utf-8"
    ) as dst:
        for line in src:
            record = line[:6].strip()
            if record in _LIGHTDOCK_ALLOWED_PDB_RECORDS:
                dst.write(line)
                kept_records += 1
                if record.startswith("END"):
                    saw_end = True
                continue

            if record == "HETATM":
                removed_hetatm += 1
            elif line.strip():
                removed_other += 1

        if not saw_end:
            dst.write("END\n")

    return {
        "input_path": input_pdb,
        "output_path": output_pdb,
        "kept_records": kept_records,
        "removed_hetatm": removed_hetatm,
        "removed_other": removed_other,
    }


def _stage_cleaned_lightdock_inputs(
    receptor_pdb: str,
    ligand_pdb: str,
    working_dir: str,
):
    """Create cleaned protein-only copies under the run directory."""
    clean_dir = os.path.join(working_dir, "cleaned_inputs")
    os.makedirs(clean_dir, exist_ok=True)

    receptor_clean = os.path.join(clean_dir, os.path.basename(receptor_pdb))
    ligand_clean = os.path.join(clean_dir, os.path.basename(ligand_pdb))

    receptor_report = clean_pdb_for_lightdock(receptor_pdb, receptor_clean)
    ligand_report = clean_pdb_for_lightdock(ligand_pdb, ligand_clean)

    return receptor_clean, ligand_clean, {
        "receptor": receptor_report,
        "ligand": ligand_report,
    }


def _lightdock_simulation_outputs(working_dir: str, steps: int) -> list[str]:
    """Return all swarm simulation outputs for the requested step count."""
    pattern = os.path.join(working_dir, "swarm_*", f"gso_{steps}.out")
    return sorted(glob.glob(pattern))


def _unsupported_lightdock_residue(simulation_output: str | None):
    """Return the unsupported residue marker reported by LightDock, if any."""
    if not simulation_output:
        return None
    match = _LIGHTDOCK_UNSUPPORTED_RESIDUE_RE.search(simulation_output)
    if not match:
        return None
    return match.group("residue"), match.group("atom")


def _raise_if_lightdock_simulation_failed(
    working_dir: str,
    steps: int,
    simulation_output,
):
    """Fail fast when LightDock produced no simulation outputs."""
    if simulation_output is None:
        return
    if _lightdock_simulation_outputs(working_dir, steps):
        return

    unsupported = _unsupported_lightdock_residue(simulation_output)
    if unsupported:
        residue, atom = unsupported
        raise LightDockSimulationError(
            "LightDock simulation failed because the scoring function "
            f"rejected unsupported residue {residue} (atom {atom}).",
            cleanable=True,
            unsupported_residue=residue,
            simulation_output=simulation_output,
        )

    raise LightDockSimulationError(
        "LightDock simulation failed because no swarm outputs were produced.",
        simulation_output=simulation_output,
    )


def _is_interactive_session() -> bool:
    """Return True when the command is attached to an interactive terminal."""
    return sys.stdin.isatty() and sys.stdout.isatty()


def _should_clean_and_retry_interactively(
    exc: LightDockSimulationError,
    working_dir: str,
) -> bool:
    """Prompt the user to clean and retry when interactive."""
    residue = exc.unsupported_residue or "an unsupported residue"
    print(
        "\nLightDock rejected a non-protein residue during scoring: "
        f"{residue}."
    )
    print(
        "PPInsight can create protein-only copies of the input PDBs under:\n"
        f"  {os.path.join(working_dir, 'cleaned_inputs')}"
    )
    print("Original input files will not be modified.")

    try:
        reply = input("Clean the PDBs and retry LightDock? [y/N]: ")
    except EOFError:
        return False
    return reply.strip().lower() in {"y", "yes"}


def _print_cleaned_input_report(cleaned_inputs):
    """Print a short summary of generated cleaned PDB copies."""
    print("\nCreated LightDock-ready protein-only copies:")
    for role in ("receptor", "ligand"):
        report = cleaned_inputs[role]
        print(
            f"  {role.capitalize()}: {report['output_path']} "
            f"(removed {report['removed_hetatm']} HETATM records, "
            f"{report['removed_other']} other records)"
        )


def _print_lightdock_failure_and_exit(exc: Exception):
    """Emit a user-facing LightDock error with actionable guidance."""
    print(f"\nERROR: {exc}", file=sys.stderr)
    if isinstance(exc, LightDockSimulationError) and exc.cleanable:
        print(
            "Hint: rerun with --auto-clean-pdb to create protein-only copies "
            "and retry automatically, or clean the PDB manually.",
            file=sys.stderr,
        )
    sys.exit(1)


def _run_lightdock_attempt(receptor: str, ligand: str, workdir: str, opts):
    """Run one LightDock attempt after clearing prior run artifacts."""
    print(f"\nCleaning previous LightDock outputs in: {workdir} (if any)")
    file_count, folder_count = _cleanup_previous_outputs(workdir)
    print(f"\nRemoved {file_count} files and {folder_count} folders.")

    print(
        "\nRunning LightDock global-docking pipeline for:\n"
        f" Receptor: {receptor}\n"
        f" Ligand:   {ligand}\n"
        f" Output dir: {workdir}\n",
    )
    lightdock_pipeline(receptor, ligand, working_dir=workdir, opts=opts)


def _print_lightdock_success_summary(
    workdir: str,
    opts,
    *,
    cleaned_retry_used: bool = False,
):
    """Print the final successful-run summary for the CLI."""
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
    ]
    if cleaned_retry_used:
        summary_lines.append(
            "- Input cleanup: retried with protein-only copies under cleaned_inputs/"
        )
    summary_lines.append(f"Results directory: {workdir}")
    print("\n".join(summary_lines))


def run_command(cmd, cwd=None):
    """
    Print and execute a shell command.
    - cmd: list of strings, e.g. ["lightdock3.py", "setup.json", "100"]
    - cwd: directory in which to run the command
    """
    print(">>", " ".join(cmd))
    try:
        proc = subprocess.Popen(
            cmd,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
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

    output_lines: list[str] = []
    assert proc.stdout is not None
    for line in proc.stdout:
        print(line, end="")
        output_lines.append(line)

    return_code = proc.wait()
    output = "".join(output_lines)
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, cmd, output=output)
    return output


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
    try:
        return run_command(cmd, cwd=working_dir)
    except subprocess.CalledProcessError as exc:
        unsupported = _unsupported_lightdock_residue(exc.output)
        if unsupported:
            residue, atom = unsupported
            raise LightDockSimulationError(
                "LightDock simulation failed because the scoring function "
                f"rejected unsupported residue {residue} (atom {atom}).",
                cleanable=True,
                unsupported_residue=residue,
                simulation_output=exc.output or "",
            ) from exc
        raise


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
    simulation_output = _run_lightdock_simulation(
        working_dir, steps,
        opts.get("cores", 1),
        opts.get("swarm_list"),
        scoring=opts.get("scoring"),
    )
    _raise_if_lightdock_simulation_failed(working_dir, steps, simulation_output)

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
    clean_group = parser.add_mutually_exclusive_group()
    clean_group.add_argument(
        "--auto-clean-pdb",
        action="store_true",
        help=(
            "If LightDock fails because DFIRE rejects unsupported non-protein "
            "residues (for example ions or other HETATM records), create "
            "protein-only copies under the run directory and retry "
            "automatically.  Use this for scripted runs where interactive "
            "prompts would block.  Original input files are not modified."
        ),
    )
    clean_group.add_argument(
        "--no-clean-pdb",
        action="store_true",
        help=(
            "Disable the interactive clean-and-retry prompt when LightDock "
            "rejects unsupported residues.  Use this when you want the "
            "command to fail immediately so you can inspect the input files "
            "manually."
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

    cleaned_retry_used = False

    try:
        _run_lightdock_attempt(receptor, ligand, workdir, opts)
    except LightDockSimulationError as exc:
        if not exc.cleanable:
            _print_lightdock_failure_and_exit(exc)

        should_retry = False
        if args.auto_clean_pdb:
            should_retry = True
        elif not args.no_clean_pdb and _is_interactive_session():
            should_retry = _should_clean_and_retry_interactively(exc, workdir)

        if not should_retry:
            _print_lightdock_failure_and_exit(exc)

        (
            cleaned_receptor,
            cleaned_ligand,
            cleaned_inputs,
        ) = _stage_cleaned_lightdock_inputs(
            receptor,
            ligand,
            workdir,
        )
        _print_cleaned_input_report(cleaned_inputs)
        print("\nRetrying LightDock with cleaned protein-only copies...")

        try:
            _run_lightdock_attempt(cleaned_receptor, cleaned_ligand, workdir, opts)
        except (LightDockSimulationError, subprocess.CalledProcessError) as retry_exc:
            _print_lightdock_failure_and_exit(retry_exc)
        cleaned_retry_used = True
    except subprocess.CalledProcessError as exc:
        _print_lightdock_failure_and_exit(exc)

    _print_lightdock_success_summary(
        workdir,
        opts,
        cleaned_retry_used=cleaned_retry_used,
    )


if __name__ == '__main__':
    main()
