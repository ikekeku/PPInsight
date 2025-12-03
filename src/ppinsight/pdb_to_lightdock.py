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

import argparse  # for CLI flags
import glob     # to find files using wildcards
import os        # for paths and directories
import subprocess  # to run LightDock command-line tools
import sys       # to read command-line arguments
import shutil    # to copy files


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

    # Make sure the base root folder exists:
    # examples/ppinsight_data/output_files/
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


def lightdock_pipeline(receptor_pdb, ligand_pdb,
                       working_dir,
                       swarms=None, glowworms=None,
                       steps=100, swarm_list=None, cores=1,
                       generate_models=False):
    """
    Runs the full LightDock workflow:

    1) Setup  (lightdock3_setup.py -> creates setup.json and initial files)
    2) Simulation (lightdock3.py)
    3) Model generation (lgd_generate_conformations.py)

    Parameters:
      receptor_pdb, ligand_pdb : input PDB file paths
      working_dir              : where all LightDock files will be written
      swarms, glowworms        : docking parameters (None → LightDock defaults)
      steps                    : number of LightDock steps (tutorial uses 100)
      swarm_list               : which swarms to run (e.g. [0])
      cores                    : number of CPU cores for the simulation
    """

    # Build paths to the copies of the PDB files INSIDE the working directory.
    # Example:
    #   working_dir = examples/ppinsight_data/output_files/lightdock_runs/2UUY_rec_vs_2UUY_lig/
    #   rec = .../2UUY_rec.pdb (inside that directory)
    #   lig = .../2UUY_lig.pdb
    rec = os.path.join(working_dir, os.path.basename(receptor_pdb))
    lig = os.path.join(working_dir, os.path.basename(ligand_pdb))

    # Copy the original PDB files into the working directory
    shutil.copy(receptor_pdb, rec)
    shutil.copy(ligand_pdb, lig)

    # -----------------------------
    # 1) LightDock setup
    # -----------------------------
    # lightdock3_setup.py prepares LightDock input files.
    # Flags:
    #   --noxt : remove terminal OXT atoms
    #   --noh  : remove hydrogens
    #   --now  : remove water molecules
    # These are recommended preprocessing steps.
    # When we run the command with cwd=working_dir, pass filenames relative
    # to that working directory (not the full path). Using the basenames
    # avoids the subprocess looking for a path that is relative to the
    # working directory (which would otherwise duplicate the run_dir prefix
    # and make the files look missing).
    cmd = ["lightdock3_setup.py", os.path.basename(rec), os.path.basename(lig), "--noxt", "--noh", "--now", "-anm"]

    # If user passed a custom number of swarms, add "-s <value>"
    if swarms:
        cmd += ["-s", str(swarms)]

    # If user passed a custom number of glowworms, add "-g <value>"
    if glowworms:
        cmd += ["-g", str(glowworms)]

    # Run the setup command inside the working directory
    run_command(cmd, cwd=working_dir)

    # -----------------------------
    # 2) Run LightDock simulation
    # -----------------------------
    # "setup.json" is created by the previous step in working_dir.
    # Here we tell LightDock how many steps to run.
    cmd = ["lightdock3.py", "setup.json", str(steps), "-c", str(cores)]

    # If a list of swarms is provided, pass it using "-l"
    # e.g., swarm_list = [0] means "only run swarm 0"
    if swarm_list:
        cmd += ["-l"] + list(map(str, swarm_list))

    # Run the simulation
    run_command(cmd, cwd=working_dir)

    # -----------------------------
    # 3) Generate docked models
    # -----------------------------
    # After the simulation, LightDock produces files like:
    #   swarm_0/gso_100.out (if steps=100)
    #
    # Here we choose:
    #   - swarm 0 (simple default)
    #   - last step (gso_<steps>.out)
    swarm = 0
    # Use a path relative to working_dir when invoking the generator so
    # the subprocess finds the swarm output file correctly.
    out_file = os.path.join(f"swarm_{swarm}", f"gso_{steps}.out")

    # Number of models to generate:
    # if glowworms is set, use that; otherwise use 200 as a standard default.
    num_models = glowworms if glowworms else 200

    # Command to generate the PDB conformations from the .out file
    cmd = [
        "lgd_generate_conformations.py",
        os.path.basename(rec),          # receptor PDB filename inside working_dir
        os.path.basename(lig),          # ligand PDB filename inside working_dir
        out_file,     # relative path to LightDock output with swarm positions
        str(num_models)  # how many models to generate
    ]

    # Run model generation only if requested by the caller
    if generate_models:
        run_command(cmd, cwd=working_dir)
    else:
        print("\nSkipping model generation (lgd_generate_conformations.py). Pass --generate to enable it.")

    # Print where the final docked models are located
    print("\nDone.")
    print("Docked models saved in:", os.path.join(working_dir, f"swarm_{swarm}"))


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        description="Run LightDock for a receptor-ligand pair."
    )
    parser.add_argument("receptor", help="Path to receptor PDB file")
    parser.add_argument("ligand", help="Path to ligand PDB file")
    parser.add_argument("--swarms", type=int, default=None, help="Number of swarms (LightDock default if not set)")
    parser.add_argument("--glowworms", type=int, default=None, help="Number of glowworms (LightDock default if not set)")
    parser.add_argument("--steps", type=int, default=100, help="Number of LightDock steps (default: 100)")
    parser.add_argument("--swarm-list", type=str, default="0", help="Comma-separated list of swarm indices to run (e.g. '0,1')")
    parser.add_argument("--cores", type=int, default=1, help="Number of CPU cores to use")
    parser.add_argument("--generate", action="store_true", help="If set, run lgd_generate_conformations.py at the end")

    args = parser.parse_args()

    receptor = args.receptor
    ligand = args.ligand
    generate_models = args.generate
    swarms = args.swarms
    glowworms = args.glowworms
    steps = args.steps
    cores = args.cores
    swarm_list = None
    if args.swarm_list:
        swarm_list = [int(x) for x in args.swarm_list.split(",") if x.strip()]

    # Create a descriptive output directory for this receptor-ligand pair.
    # This will be:
    #   examples/ppinsight_data/output_files/lightdock_runs/<rec_vs_lig>/
    workdir = make_output_dir(receptor, ligand, method="lightdock_runs")
    
    print(f"\nCleaning previous LightDock outputs in: {workdir} (if any)")
    items = glob.glob(os.path.join(workdir, "lightdock*")) # find all files/folders starting with "lightdock"
    items += glob.glob(os.path.join(workdir, "swarm_*")) # find all files/folders starting with "swarm_"
    items += [os.path.join(workdir, "setup.json"), os.path.join(workdir, "init")] # add specific files/folders to delete

    file_count = folder_count = 0
    for p in items: # iterate over all found items
        if os.path.exists(p):
            if os.path.isdir(p): # if it's a directory, remove it and its contents
                shutil.rmtree(p, ignore_errors=True) # ignore errors if directory doesn't exist
                folder_count += 1
            elif os.path.isfile(p): # if it's a file, remove it
                os.remove(p)
                folder_count += 1
    print(f"\nRemoved {file_count} files and {folder_count} folders.")

    # Run the LightDock pipeline with mostly default settings    
    print(f"\nRunning LightDock pipeline for:\n Receptor: {receptor}\n Ligand:   {ligand}\n Output dir: {workdir}\n")
    lightdock_pipeline(
        receptor,
        ligand,
        working_dir=workdir,
        swarms=swarms,
        glowworms=glowworms,
        steps=steps,
        swarm_list=swarm_list,
        cores=cores,
        generate_models=generate_models,
    )

    # -----------------------------
    # Summary
    # -----------------------------
    summary = f"""
    LightDock docking completed.
    Run folder: {workdir}
    Parameters used:
    - Swarms: {swarms if swarms else 'default'}
    - Glowworms: {glowworms if glowworms else 'default'}
    - Steps: {steps}
    - CPU cores: {cores}
    Docked models (if generated) are in: {os.path.join(workdir, 'swarm_0')}
    """
    print(summary.strip())
