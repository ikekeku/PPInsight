"""
pdb_to_haddock.py
-----------------
Prepare a HADDOCK3 run directory and a tutorial-style config (.cfg).

This script DOES NOT run HADDOCK by default. It only stages input files
and writes a config file the user can run manually or with --run.

It mirrors the structure of the notebook version and has a small CLI.
"""
from pathlib import Path
import argparse
import shutil
import os
import subprocess
import platform
import shlex


def info(msg: str):
    print(f"[INFO] {msg}")


def run_command(cmd, cwd=None):
    """Print and run a command. Kept small so tests can monkeypatch it."""
    print(">>", " ".join(map(str, cmd)))
    subprocess.run(list(map(str, cmd)), cwd=cwd, check=True)


BASE_ROOT = Path("examples/ppinsight_data/output_files")
METHOD = "haddock_runs"


def make_run_dir(rec, lig, runname, base_root=BASE_ROOT, method=METHOD):
    base_root = Path(base_root)
    base_root.mkdir(parents=True, exist_ok=True)

    method_dir = base_root / method
    method_dir.mkdir(exist_ok=True)

    run_dir = method_dir / f"{Path(rec).stem}_vs_{Path(lig).stem}"
    run_dir.mkdir(exist_ok=True)

    data_dir = run_dir / "data"
    data_dir.mkdir(exist_ok=True)

    info(f"Created run folder: {run_dir}")
    return run_dir, data_dir


def copy_inputs(data_dir, rec, lig, ambig=None):
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


def haddock_pipeline(rec, lig, runname="run1", mode="local", ncores=4, ambig=None, base_root=BASE_ROOT, method=METHOD, run_haddock=False, haddock_cmd: str = "haddock3"):
    # 1) Create organized output folder
    run_dir, data_dir = make_run_dir(rec, lig, runname, base_root=base_root, method=method)

    # 2) Copy input files
    rec_dst, lig_dst, ambig_dst = copy_inputs(data_dir, rec, lig, ambig)

    # data-relative paths inside .cfg
    rec_rel = f"data/{rec_dst.name}"
    lig_rel = f"data/{lig_dst.name}"
    ambig_rel = f"data/{ambig_dst.name}" if ambig_dst else ""

    # 3) Write config
    cfg_path = run_dir / f"{runname}.cfg"
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
    write_cfg(cfg_path, runname, mode, ncores, rec_rel, lig_rel, ambig_rel)

    # Optionally run haddock3 with the config file
    if run_haddock:
        # Before running, try to detect a common failure: an included CNS
        # binary may be built for a different CPU architecture (e.g. arm64
        # binary on an x86_64 machine). If so, fail early with a clear message.
        try:
            import haddock as _haddock_mod
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

    info("\nDone.")
    info(f"Files staged in: {data_dir}")
    info(f"Config file: {cfg_path}")
    return run_dir, cfg_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare HADDOCK3 config and input files (tutorial style)")
    parser.add_argument("receptor", help="Path to receptor PDB file")
    parser.add_argument("ligand", help="Path to ligand PDB file")
    parser.add_argument("--runname", default="run1", help="Config run_dir name")
    parser.add_argument("--mode", default="local", help="HADDOCK execution mode")
    parser.add_argument("--ncores", type=int, default=4, help="CPU cores")
    parser.add_argument("--ambig", default=None, help="Path to ambiguous restraints (.tbl)")
    parser.add_argument("--run", action="store_true", help="If set, run haddock3 with the generated config")
    parser.add_argument("--haddock-cmd", default="haddock3", help="Command to invoke HADDOCK (default: haddock3). Can be a wrapper or 'echo' for testing.")

    args = parser.parse_args()

    haddock_pipeline(args.receptor, args.ligand, runname=args.runname, mode=args.mode, ncores=args.ncores, ambig=args.ambig, run_haddock=args.run, haddock_cmd=args.haddock_cmd)
