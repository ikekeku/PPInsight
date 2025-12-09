"""
pdb_to_haddock.py
-----------------
Prepare a HADDOCK3 run directory and a tutorial-style config (.cfg).

This script DOES NOT run HADDOCK by default. It only stages input files
and writes a config file the user can run manually or with --run.

It mirrors the structure of the notebook version and has a small CLI.
"""
from pathlib import Path
from typing import Optional
import argparse
import shutil
import os
import subprocess
import platform
import shlex
import glob
import time
import sys


def _project_root():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))


def resolve_input_path(path: str) -> str:
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
            info(f"Resolved '{path}' -> '{found}'")
            return found
    raise FileNotFoundError(f"Could not find input file '{path}' (searched repo {proj}).")


def info(msg: str):
    """Simple logger for informational messages.

    Kept minimal so tests can run with predictable output.
    """
    print(f"[INFO] {msg}")


def run_command(cmd, cwd=None):
    """Print and run a command.

    This is a thin wrapper around subprocess.run so tests can monkeypatch
    `run_command` to avoid executing external binaries.
    """
    print(">>", " ".join(map(str, cmd)))
    subprocess.run(list(map(str, cmd)), cwd=cwd, check=True)


BASE_ROOT = Path("examples/ppinsight_data/output_files")
METHOD = "haddock_runs"
CONTAINER_IMAGE = "ghcr.io/haddocking/haddock3:latest"


def make_run_dir(rec, lig, runname, base_root=BASE_ROOT, method=METHOD):
    """Create a run directory and its data subdirectory.

    If directories already exist, they are reused. Return a tuple of
    (run_dir, data_dir).
    """
    base_root = Path(base_root)
    base_root.mkdir(parents=True, exist_ok=True)

    method_dir = base_root / method
    method_dir.mkdir(exist_ok=True)
    # Use the explicit runname when provided so the on-disk directory
    # matches the "run_dir" value written into the HADDOCK .cfg file.
    # Fall back to the receptor_vs_ligand pattern if runname is falsy.
    if runname:
        run_dir = method_dir / runname
    else:
        run_dir = method_dir / f"{Path(rec).stem}_vs_{Path(lig).stem}"
    run_dir.mkdir(exist_ok=True)

    data_dir = run_dir / "data"
    data_dir.mkdir(exist_ok=True)

    info(f"Created run folder: {run_dir}")
    return run_dir, data_dir


def copy_inputs(data_dir, rec, lig, ambig=None):
    """Copy receptor, ligand (and optionally ambig) into the run data dir.

    Accepts either full paths or short basenames; callers should resolve
    inputs using `resolve_input_path` if they want repository searching.
    """
    # Validate inputs and give helpful errors if files are missing
    rec_path = Path(rec)
    lig_path = Path(lig)
    if not rec_path.exists():
        raise FileNotFoundError(f"Receptor file not found: {rec}")
    if not lig_path.exists():
        raise FileNotFoundError(f"Ligand file not found: {lig}")

    rec_dst = data_dir / rec_path.name
    lig_dst = data_dir / lig_path.name
    shutil.copy(str(rec_path), str(rec_dst))
    shutil.copy(str(lig_path), str(lig_dst))

    ambig_dst = None
    if ambig:
        ambig_path = Path(ambig)
        if not ambig_path.exists():
            raise FileNotFoundError(f"Ambiguous restraints file not found: {ambig}")
        ambig_dst = data_dir / ambig_path.name
        shutil.copy(str(ambig_path), str(ambig_dst))

    info("Copied input structures into data/ folder")
    return rec_dst, lig_dst, ambig_dst


def write_cfg(cfg_path: Path, runname: str, mode: str, ncores: int, rec_rel: str, lig_rel: str, ambig_rel: str):
    # If no ambiguous restraints are provided, enable ab-initio sampling
    # for the rigid-body stage using `cmrest = true`.
    abinitio_block = "" if ambig_rel else "cmrest = true\n"

    text = f"""# ====================================================================
# Protein-protein docking example (auto-generated)

# directory in which the scoring will be done
run_dir = "{runname}"

# execution mode
mode = "{mode}"
ncores = {ncores}

# molecules to be docked
molecules =  [
    "{rec_rel}",
    "{lig_rel}"
    ]

# ====================================================================
# Parameters for each stage are defined below, prefer full paths
# ====================================================================
[topoaa]
autohis = false
[topoaa.mol1]
nhisd = 0
nhise = 1
hise_1 = 75
[topoaa.mol2]
nhisd = 1
hisd_1 = 76
nhise = 1
hise_1 = 15

[rigidbody]
tolerance = 20
ambig_fname = "{ambig_rel}"
sampling = 20

{abinitio_block}
[caprieval]
reference_fname = ""

[seletop]
select = 5

[flexref]
tolerance = 20
ambig_fname = "{ambig_rel}"

[emref]
tolerance = 20
ambig_fname = "{ambig_rel}"

[clustfcc]
min_population = 1

[seletopclusts]
top_models = 4
# ====================================================================
"""
    cfg_path.write_text(text.strip() + "\n")
    info(f"Wrote config file: {cfg_path}")


def _choose_runname(base_dir: Path, base_name: str, overwrite: bool) -> str:
    """Choose a runname when none is provided.

    If overwrite is False and a directory with base_name exists under base_dir,
    append a numeric suffix (_1, _2, ...) to avoid clobbering previous runs.
    """
    candidate = base_name
    method_dir = base_dir
    i = 1
    while (method_dir / candidate).exists():
        if overwrite:
            # caller will remove it explicitly
            break
        candidate = f"{base_name}_{i}"
        i += 1
    return candidate


def haddock_pipeline(rec, lig, runname: str = None, mode: str = "local", ncores: int = 4,
                    ambig: str = None, base_root=BASE_ROOT, method=METHOD,
                    run_haddock: bool = False, haddock_cmd: str = "haddock3",
                    overwrite: bool = False,
                    container: str = "auto", container_image: str = CONTAINER_IMAGE,
                    workspace_root: Optional[str] = None):
    """Prepare a HADDOCK run directory, copy inputs, write cfg, and optionally run HADDOCK.

    - rec/lig may be full paths or short basenames; they will be resolved against
      the repository if not found on the filesystem.
    - If runname is None, a sensible default based on the receptor/ligand names
      is chosen; if that directory already exists, a numeric suffix is appended
      unless overwrite=True.
    """
    # Resolve inputs so short basenames work (search the repo)
    rec = resolve_input_path(rec)
    lig = resolve_input_path(lig)

    # Determine runname: if not provided, base it on receptor and ligand
    base_name = f"{Path(rec).stem}_vs_{Path(lig).stem}"
    method_root = Path(base_root) / method
    method_root.mkdir(parents=True, exist_ok=True)
    chosen_runname = runname if runname else _choose_runname(method_root, base_name, overwrite)

    run_dir = method_root / chosen_runname
    if run_dir.exists() and overwrite:
        # remove existing directory tree to allow a fresh run
        shutil.rmtree(run_dir)

    # 1) Create organized output folder
    run_dir, data_dir = make_run_dir(rec, lig, chosen_runname, base_root=base_root, method=method)

    # 2) Copy input files (they are already resolved)
    rec_dst, lig_dst, ambig_dst = copy_inputs(data_dir, rec, lig, ambig)

    # data-relative paths inside .cfg
    rec_rel = f"data/{rec_dst.name}"
    lig_rel = f"data/{lig_dst.name}"
    ambig_rel = f"data/{ambig_dst.name}" if ambig_dst else ""

    # 3) Write config
    cfg_path = run_dir / f"{chosen_runname}.cfg"
    # If there are existing .cfg files in the run directory, report and remove
    existing_cfgs = list(run_dir.glob("*.cfg"))
    if existing_cfgs:
        info(f"Found existing cfg files in {run_dir}:")
        for p in existing_cfgs:
            info(f" - {p.name}")
        # Remove them so the new cfg is the authoritative one for this run
        for p in existing_cfgs:
            try:
                p.unlink()
                info(f"Removed old cfg: {p.name}")
            except Exception:
                info(f"Could not remove existing cfg: {p.name}; leaving it in place")
    write_cfg(cfg_path, chosen_runname, mode, ncores, rec_rel, lig_rel, ambig_rel)

    # Optionally run haddock3 with the config file
    if run_haddock:
        # Before running, try to detect a common failure: an included CNS
        # binary may be built for a different CPU architecture (e.g. arm64
        # binary on an x86_64 machine). If so, fail early with a clear message.
        try:
            import haddock as _haddock_mod  # type: ignore
            haddock_dir = Path(_haddock_mod.__file__).resolve().parent
            cns_path = haddock_dir / "bin" / "cns"
            if cns_path.exists():
                # Use the `file` utility to inspect the binary's arch.
                try:
                    out = subprocess.run(["file", str(cns_path)], capture_output=True, text=True, check=True)
                    arch_info = out.stdout.lower()
                    host_arch = platform.machine().lower()
                    # Simple checks for common mismatches
                    if ("arm64" in arch_info or "aarch64" in arch_info) and ("arm" not in host_arch and "aarch" not in host_arch):
                        raise RuntimeError(
                            f"The HADDOCK cns binary at {cns_path} is built for ARM (arm64) but the host CPU is {host_arch}.\n"
                            "This will cause a 'Bad CPU type in executable' error.\n"
                            "Possible fixes: install a HADDOCK build for your architecture, run on a matching machine, or use a compatible conda package."
                        )
                    if ("intel" in arch_info or "x86_64" in arch_info or "mach-o 64-bit executable x86_64" in arch_info) and ("x86_64" not in host_arch and "amd64" not in host_arch):
                        raise RuntimeError(
                            f"The HADDOCK cns binary at {cns_path} appears to be x86_64 but the host CPU is {host_arch}.\n"
                            "This architecture mismatch will prevent execution.\n"
                        )
                except subprocess.CalledProcessError:
                    # If `file` is not available or fails, continue and let haddock report
                    pass
        except Exception as exc:
            # If we detected a clear incompatibility, surface a helpful error
            if isinstance(exc, RuntimeError):
                raise
            # Otherwise ignore and attempt to run haddock (it may still work)

        # When running with cwd=run_dir, pass the filename relative to that dir
        # Allow haddock_cmd to be a wrapper or an alternate command (e.g. 'echo').
        # If the user requested containerized execution or it can be auto-selected,
        # prefer running haddock3 inside a container to avoid host GLIBC/CNS ABI issues.

        def _detect_container_runtime():
            """Detect available container runtimes.

            Returns 'docker', 'apptainer' (or 'singularity'), or None.
            """
            if shutil.which("docker"):
                return "docker"
            # Apptainer may be installed under 'apptainer' or legacy 'singularity'
            if shutil.which("apptainer"):
                return "apptainer"
            if shutil.which("singularity"):
                return "apptainer"
            return None

        def _run_in_container(runtime: str, image: str, host_workspace: str, run_dir_path: Path, cfg_name: str):
            """Construct and run a container command that executes haddock3 inside the image.

            This function uses run_command to preserve the testability.
            """
            host_workspace = os.path.abspath(host_workspace)
            # Compute run_dir relative to the workspace root so we can `cd` correctly inside container
            rel_run = os.path.relpath(str(run_dir_path), host_workspace)
            # Ensure paths are safe for commands
            rel_run_posix = rel_run.replace(os.path.sep, "/")
            if runtime == "docker":
                # Mount the workspace into /workspace in the container
                workdir = f"/workspace/{os.path.dirname(rel_run_posix)}" if os.path.dirname(rel_run_posix) else "/workspace"
                cfg_base = os.path.basename(cfg_name)
                cmd = [
                    "docker", "run", "--rm",
                    "-v", f"{host_workspace}:/workspace",
                    "-w", workdir,
                    image,
                    "bash", "-lc", f"haddock3 {shlex.quote(cfg_base)}"
                ]
            else:
                # Apptainer can execute docker images directly via the docker:// prefix
                image_spec = image
                if not image_spec.startswith("docker://") and not image_spec.startswith("shub://"):
                    image_spec = f"docker://{image_spec}"
                # apptainer exec --bind host:/workspace --pwd /workspace docker://image bash -lc 'cd ... && haddock3 cfg'
                cmd = [
                    "apptainer", "exec",
                    "--bind", f"{host_workspace}:/workspace",
                    "--pwd", "/workspace",
                    image_spec,
                    "bash", "-lc", f"cd /workspace/{rel_run_posix} && haddock3 {shlex.quote(os.path.basename(cfg_name))}"
                ]
            run_command(cmd)

        # Decide container use
        chosen_container = None
        if container == "auto":
            chosen_container = _detect_container_runtime()
        elif container in ("docker", "apptainer", "singularity"):
            chosen_container = container if container != "singularity" else "apptainer"

        if chosen_container:
            # Determine a sensible host workspace root for binding; fall back to project root
            host_ws = workspace_root if workspace_root else _project_root()
            try:
                _run_in_container(chosen_container, container_image, host_ws, run_dir, cfg_path.name)
            except subprocess.CalledProcessError:
                # If the container run failed, surface the error to the caller (same behavior as local run)
                raise
        else:
            # No container runtime found; run locally
            cmd = shlex.split(haddock_cmd) + [cfg_path.name]
            try:
                run_command(cmd, cwd=run_dir)
            except subprocess.CalledProcessError as e:
                # Inspect stderr/msg to see if it's a Bad CPU type issue and suggest fixes
                msg = str(e)
                if "Bad CPU type" in msg or "bad cpu" in msg.lower():
                    raise RuntimeError(
                        "HADDOCK (CNS) failed to execute due to CPU architecture mismatch.\n"
                        "Check that the 'cns' binary in the haddock package matches your machine architecture.\n"
                        "If you installed haddock via conda, try installing a variant built for your platform or use a different channel."
                    ) from e
                raise

    return run_dir, cfg_path


def main(argv=None):
    parser = argparse.ArgumentParser(description="Prepare HADDOCK3 config and input files (tutorial style)")
    parser.add_argument("receptor", help="Path to receptor PDB file")
    parser.add_argument("ligand", help="Path to ligand PDB file")
    parser.add_argument("--runname", default=None, help="Config run_dir name (default: auto-increment to avoid clobbering)")
    parser.add_argument("--mode", default="local", help="HADDOCK execution mode")
    parser.add_argument("--ncores", type=int, default=4, help="CPU cores")
    parser.add_argument("--ambig", default=None, help="Path to ambiguous restraints (.tbl)")
    parser.add_argument("--run", action="store_true", help="If set, run haddock3 with the generated config")
    parser.add_argument("--haddock-cmd", default="haddock3", help="Command to invoke HADDOCK (default: haddock3). Can be a wrapper or 'echo' for testing.")
    parser.add_argument("--overwrite", action="store_true", help="If set, overwrite an existing run directory with the same name")
    parser.add_argument("--container", default="auto", help="Container runtime to use: 'auto', 'docker', 'apptainer' (default: auto)")
    parser.add_argument("--container-image", default=CONTAINER_IMAGE, help="Container image to use when running HADDOCK inside a container")

    args = parser.parse_args(argv)

    # Resolve receptor/ligand paths (accept basenames or full paths)
    try:
        receptor_path = resolve_input_path(args.receptor)
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        sys.exit(2)

    try:
        ligand_path = resolve_input_path(args.ligand)
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        sys.exit(2)

    # Run the pipeline
    run_dir, cfg_path = haddock_pipeline(
        receptor_path,
        ligand_path,
        runname=args.runname,
        mode=args.mode,
        ncores=args.ncores,
        ambig=args.ambig,
        run_haddock=args.run,
        haddock_cmd=args.haddock_cmd,
        overwrite=args.overwrite,
        container=args.container,
        container_image=args.container_image,
        workspace_root=_project_root(),
    )

    # User-facing summary specific to HADDOCK
    summary_lines = [
        "HADDOCK staging completed.",
        f"Run name: {run_dir.name}",
        f"Run directory: {run_dir}",
        f"Config file: {cfg_path}",
        f"Mode: {args.mode}",
        f"CPU cores: {args.ncores}",
        f"Receptor (staged): {run_dir / 'data' / Path(receptor_path).name}",
        f"Ligand   (staged): {run_dir / 'data' / Path(ligand_path).name}",
        f"Ambig restraints: {args.ambig if args.ambig else 'None'}",
        f"HADDOCK command to be executed: {args.haddock_cmd if args.run else '(not running)'}",
    ]
    print('\n'.join(summary_lines))


if __name__ == '__main__':
    main()
