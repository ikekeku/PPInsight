Running HADDOCK on Hyak (or other Linux x86_64 clusters)
======================================================

This short note explains how to prepare and run HADDOCK on an HPC system such as Hyak (which provides Linux x86_64 compute nodes). It also documents how to check and avoid the "Bad CPU type in executable" error we encountered when a macOS-conda-installed HADDOCK contained an arm64 CNS binary.

Summary
-------
- The issue: locally-installed HADDOCK on macOS may include a CNS binary built for arm64 while the host is x86_64 (or vice versa). Running such a binary produces "Bad CPU type in executable" or similar errors.
- The solution: run HADDOCK on a system whose CPU architecture matches the HADDOCK/CNS binaries (Linux x86_64 on Hyak), or use a container image built for the target architecture.

Recommended options
-------------------
1) Run on Hyak (recommended when you have access)
   - Use a Hyak login node to create a conda env and install HADDOCK, or (better) use the cluster's provided software modules / container images.
   - Submit a job script that runs `haddock3` inside the compute node where Linux x86_64 is available.

2) Use a container (Docker / Singularity)
   - Pull an x86_64-compatible image containing haddock and cns. This avoids architecture mismatch on mixed-architecture hosts.

3) Rebuild or install the correct variant of HADDOCK for your architecture
   - If you must run locally, install a conda package built for your CPU (e.g. on macOS M1/M2 use arm64 builds; on x86_64 use Intel builds).

Quick checks to run before execution
-----------------------------------
- Inspect the CNS binary inside the haddock installation:

```bash
file $(python -c "import haddock; import pathlib; print(pathlib.Path(haddock.__file__).resolve().parent / 'bin' / 'cns')")
```

- Check the machine architecture:

```bash
uname -m
```

If the CNS binary is "Mach-O 64-bit executable arm64" but `uname -m` reports "x86_64`, you'll get the Bad CPU error.

Hyak example job script (SLURM)
--------------------------------
Below is a small SLURM job script you can adapt for Hyak. It assumes haddock is available via a conda environment or module on the compute node.

```bash
#!/bin/bash
#SBATCH --job-name=haddock_job
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=haddock-%j.out

# Load modules or activate your conda env
# module load miniconda
# conda activate haddock-env

# Optional: check CNS binary arch before running
CNS_PATH=$(python -c "import haddock, pathlib; print(pathlib.Path(haddock.__file__).resolve().parent / 'bin' / 'cns')")
if [ -f "$CNS_PATH" ]; then
  echo "CNS binary:" $(file "$CNS_PATH")
fi

# Run haddock from the run directory. Use relative cfg filename.
cd /path/to/your/haddock/run_dir
haddock3 run1.cfg
```

Notes
-----
- On Hyak, prefer running inside the compute node (via sbatch) rather than on login nodes.
- If using Singularity, bind your run directory into the container and run haddock3 inside it.
- If you see architecture mismatch, prefer container or matching conda package rather than trying to mix architectures.

Contact
-------
If you'd like, I can prepare a ready-to-use Singularity definition or Dockerfile that packages an x86_64 haddock + cns so you can run it anywhere with consistent architecture.