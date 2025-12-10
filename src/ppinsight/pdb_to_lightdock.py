"""
pdb_to_lightdock.py
-------------------
Runs LightDock for a receptor-ligand pair.

Input PDBs: you pass them on the command line, e.g.
    python pdb_to_lightdock.py examples/ppinsight_data/input_files/2UUY_rec.pdb \
                               examples/ppinsight_data/input_files/2UUY_lig.pdb

Output folders:
    examples/ppinsight_data/output_files/lightdock_runs/<receptor>_vs_<ligand>/
"""

# typing.Optional is not used in this small module
import argparse  # for CLI flags
import glob     # to find files using wildcards
import os        # for paths and directories
import subprocess  # to run LightDock command-line tools
import sys       # to read command-line arguments
import shutil    # to copy files


def _project_root():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))


def resolve_input_path(path):
    """Resolve an input path by returning it if it exists or searching the repo for the basename.

    Accepts short names like '2UUY_rec' or '2UUY_rec.pdb' and returns the absolute path
    of the first matching file found under the repository root.
    """
    if not path:
        raise FileNotFoundError("Empty input path")

    p = os.path.expanduser(path)
    p = os.path.abspath(p)
    if os.path.exists(p):
        return p

    # Search repo for basename (with and without .pdb)
    base = os.path.basename(path)
    candidates = [base]
    if not base.lower().endswith('.pdb'):
        candidates.append(base + '.pdb')

    proj = _project_root()
    for c in candidates:
        pattern = os.path.join(proj, '**', c)
        matches = glob.glob(pattern, recursive=True)
        if matches:
            found = os.path.abspath(matches[0])
            print(f"Resolved '{path}' -> '{found}'")
            return found

    raise FileNotFoundError(f"Could not find input file '{path}' (searched repo {proj}).")



def run_command(cmd, cwd=None):
    """
    Print and execute a shell command.
    - cmd: list of strings, e.g. ["lightdock3.py", "setup.json", "100"]
    - cwd: directory in which to run the command
    """
    print(">>", " ".join(cmd))          # show the command so the user sees what is happening
    subprocess.run(cmd, cwd=cwd, check=True)  # run the command; stop if it fails (check=True)


def make_output_dir(receptor_pdb, ligand_pdb,
                    base_root="examples/ppinsight_data/output_files",
                    method="lightdock_runs"):
    """
    Create an organized output directory structure.

    Final structure:
        examples/ppinsight_data/output_files/
            lightdock_runs/
                <receptor_name>_vs_<ligand_name>/

    - receptor_pdb, ligand_pdb: paths the user passed in
    - base_root: top-level output folder shared by all methods
    - method: subfolder for this tool (LightDock)

    If folders already exist, they are reused (not deleted).
    """

    # Resolve base_root relative to the project root to avoid creating
    # nested output paths when the current working directory changes
    # (e.g. running from a different folder on shared filesystems).
    if not os.path.isabs(base_root):
        base_root = os.path.join(_project_root(), base_root)
    # Make sure the base root folder exists:
    # /.../PPInsight/examples/ppinsight_data/output_files/
    os.makedirs(base_root, exist_ok=True)

    # Inside that, make sure the method folder exists:
    # e.g. examples/ppinsight_data/output_files/lightdock_runs/
    method_dir = os.path.join(base_root, method)
    os.makedirs(method_dir, exist_ok=True)

    # Get just the filenames without path and .pdb extension
    # e.g. "2UUY_rec.pdb" -> "2UUY_rec"
    rec_name = os.path.splitext(os.path.basename(receptor_pdb))[0]
    lig_name = os.path.splitext(os.path.basename(ligand_pdb))[0]

    # Final run folder name: "<receptor>_vs_<ligand>"
    # e.g. "2UUY_rec_vs_2UUY_lig"
    run_folder_name = f"{rec_name}_vs_{lig_name}"

    # Full path to this specific run
    run_dir = os.path.join(method_dir, run_folder_name)

    # Create the run directory if it does not exist (do not delete if it already exists)
    os.makedirs(run_dir, exist_ok=True)

    # Return the full path so the rest of the script can use it
    return run_dir


def _run_lightdock_setup(working_dir, rec_basename, lig_basename, swarms, glowworms):
    """Run the LightDock setup step (creates setup.json and initial files)."""
    cmd = [
        "lightdock3_setup.py",
        rec_basename,
        lig_basename,
        "--noxt",
        "--noh",
        "--now",
        "-anm",
    ]
    if swarms:
        cmd += ["-s", str(swarms)]
    if glowworms:
        cmd += ["-g", str(glowworms)]
    run_command(cmd, cwd=working_dir)


def _run_lightdock_simulation(working_dir, steps, cores, swarm_list):
    """Run the LightDock simulation step (lightdock3.py)."""
    cmd = ["lightdock3.py", "setup.json", str(steps), "-c", str(cores)]
    if swarm_list:
        cmd += ["-l"] + list(map(str, swarm_list))
    run_command(cmd, cwd=working_dir)


def _run_lightdock_generation(working_dir, rec_basename, lig_basename, steps, glowworms):
    """Run the model generation step (lgd_generate_conformations.py)."""
    swarm = 0
    out_file = os.path.join(f"swarm_{swarm}", f"gso_{steps}.out")
    num_models = glowworms if glowworms else 200
    cmd = [
        "lgd_generate_conformations.py",
        rec_basename,
        lig_basename,
        out_file,
        str(num_models),
    ]
    run_command(cmd, cwd=working_dir)


def _execute_lightdock_stages(receptor_pdb, ligand_pdb, working_dir, opts):
    """Execute LightDock setup, simulation and optional generation.

    This helper holds the commands so the public pipeline function stays small
    (keeps local variable counts low for linters).
    """
    # Run the three LightDock stages
    _run_lightdock_setup(
        working_dir,
        os.path.basename(receptor_pdb),
        os.path.basename(ligand_pdb),
        opts.get("swarms"),
        opts.get("glowworms"),
    )
    _run_lightdock_simulation(
        working_dir,
        opts.get("steps", 100),
        opts.get("cores", 1),
        opts.get("swarm_list"),
    )
    if opts.get("generate_models", False):
        _run_lightdock_generation(
            working_dir,
            os.path.basename(receptor_pdb),
            os.path.basename(ligand_pdb),
            opts.get("steps", 100),
            opts.get("glowworms"),
        )
    else:
        print("\nSkipping model generation; pass --generate to enable it.")

    print("\nDone.")
    print("Docked models saved in:", os.path.join(working_dir, "swarm_0"))


def lightdock_pipeline(receptor_pdb, ligand_pdb, working_dir, opts=None, **kwargs):
    """Run the full LightDock workflow using an opts dict or kwargs.

    Supported keys: swarms, glowworms, steps, swarm_list, cores, generate_models
    """
    opts = dict(opts or {})
    opts.update(kwargs)

    # Copy input PDBs into the working directory (preserve filename)
    shutil.copy(receptor_pdb, working_dir)
    shutil.copy(ligand_pdb, working_dir)

    # Delegate to the detailed runner to keep this function small and clear.
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
        description="Run LightDock for a receptor-ligand pair."
    )
    parser.add_argument("receptor", help="Path to receptor PDB file")
    parser.add_argument("ligand", help="Path to ligand PDB file")
    parser.add_argument(
        "--swarms",
        type=int,
        default=None,
        help=("Number of swarms (LightDock default if not set)"),
    )
    parser.add_argument(
        "--glowworms",
        type=int,
        default=None,
        help=("Number of glowworms (LightDock default if not set)"),
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=100,
        help=("Number of LightDock steps (default: 100)"),
    )
    parser.add_argument(
        "--swarm-list",
        type=str,
        default="0",
        help=("Comma-separated list of swarm indices to run (e.g. '0,1')"),
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help=("Number of CPU cores to use"),
    )
    parser.add_argument(
        "--generate",
        action="store_true",
        help=("If set, run lgd_generate_conformations.py at the end"),
    )

    args = parser.parse_args(argv)

    receptor = args.receptor
    ligand = args.ligand
    # Allow users to pass short basenames (e.g. '2UUY_rec' or '2UUY_rec.pdb').
    # If the provided path doesn't exist, search the repository for a matching file.
    try:
        receptor = resolve_input_path(receptor)
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        sys.exit(2)

    try:
        ligand = resolve_input_path(ligand)
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        sys.exit(2)
    # options are consumed via the opts dict below; avoid extra local variables
    # to keep pylint happy (pack into opts directly)

    # Create a descriptive output directory for this receptor-ligand pair.
    # This will be:
    #   examples/ppinsight_data/output_files/lightdock_runs/<rec_vs_lig>/
    workdir = make_output_dir(receptor, ligand, method="lightdock_runs")

    print(f"\nCleaning previous LightDock outputs in: {workdir} (if any)")
    file_count, folder_count = _cleanup_previous_outputs(workdir)
    print(f"\nRemoved {file_count} files and {folder_count} folders.")

    # Run the LightDock pipeline with mostly default settings
    print(
        "\nRunning LightDock pipeline for:\n"
        f" Receptor: {receptor}\n"
        f" Ligand:   {ligand}\n"
        f" Output dir: {workdir}\n",
    )
    # Pack runtime options into a single dict to keep the local variable count low
    opts = {
        "swarms": args.swarms,
        "glowworms": args.glowworms,
        "steps": args.steps,
        "swarm_list": (
            [int(x) for x in args.swarm_list.split(",") if x.strip()]
            if args.swarm_list
            else None
        ),
        "cores": args.cores,
        "generate_models": args.generate,
    }

    lightdock_pipeline(receptor, ligand, working_dir=workdir, opts=opts)

    # -----------------------------
    # Summary
    # -----------------------------
    summary_lines = [
        "LightDock docking completed.",
        f"Run folder: {workdir}",
        "Parameters used:",
        "- Swarms: " + (str(opts["swarms"]) if opts["swarms"] else "default"),
        "- Glowworms: "
        + (str(opts["glowworms"]) if opts["glowworms"] else "default"),
        "- Steps: " + str(opts["steps"]),
        "- CPU cores: " + str(opts["cores"]),
        f"Docked models (if generated) are in: {os.path.join(workdir, 'swarm_0')}",
    ]
    print("\n".join(summary_lines))


if __name__ == '__main__':
    main()
