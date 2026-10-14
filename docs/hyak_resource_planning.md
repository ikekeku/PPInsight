# Hyak Resource Planning

This guide helps estimate the Hyak resources needed for PPInsight docking jobs.
It is based on one production-scale, bound-redocking run of VEGFA with VEGFR1
and VEGFR2. Treat the recorded times as planning data, not guarantees: protein
size, input quality, queue conditions, and engine versions all affect runtime.

## What a job requests

The supplied VEGFA-VEGFR Slurm script requests one node, 24 CPU cores, 32 GB
of memory, and up to three days of wall time.

- LightDock and HADDOCK receive the requested core count.
- The current Rosetta batch runner executes trajectories serially, so it does
  not use all 24 requested cores during its phase.
- A job's wall time is the maximum time Slurm permits it to run. Slurm ends a
  job that reaches this limit, even if the current docking stage is nearly
  complete.

## Measured production run

The first VEGFA-VEGFR1 pair completed at the production settings used by the
Slurm script:

| Engine | Settings | Measured time | Raw output size |
|---|---|---:|---:|
| LightDock | 400 swarms, 200 glowworms, 100 steps, ANM off | about 11 h 25 m | about 11 GB |
| HADDOCK | 10,000 rigid-body models, 400 selected, full refinement | about 2 h 50 m | about 1.6 GB |
| Rosetta | 1,000 trajectories, relaxation and clustering enabled | about 4 h 23 m | about 376 MB |
| Total | Engines run sequentially | about 18 h 37 m | about 13 GB |

The job's measured peak resident memory was about 25 GB. A 32 GB request was
enough for this run, but it leaves limited headroom. Request 40 GB or more when
the partition allows it if the inputs are larger or you want more protection
from memory-related failure.

## Choosing wall time

For two pairs at comparable production settings, the measured first-pair time
implies at least 37 hours of docking before collection, DockQ calculations, or
unexpected retries. A 30-hour allocation is therefore too short for this
specific all-engine, two-pair workflow.

Suggested starting points:

| Goal | Suggested wall time | Notes |
|---|---:|---|
| One pair, all three engines | 24-36 hours | Based on the measured VEGFR1 run. |
| Two comparable pairs, all three engines | 72 hours | Allows for the measured runtime, post-processing, and variation between pairs. |
| Pipeline smoke test | 2-6 hours | Use PPInsight's `--screening` preset; this only verifies that the workflow works and is not production sampling. |

Use `sbatch --test-only` before submitting when a maintenance window is near.
Slurm can start a job before maintenance only when the requested wall time fits
in the available window.

## Storage and logs

LightDock writes a PDB model for every generated pose. At the production
settings above, this is 80,000 poses per pair before clustering. The resulting
output and clustering messages are expected to be large.

- Reserve at least 20 GB of scratch space per comparable pair for raw docking
  output; more is prudent for larger proteins or repeated runs.
- The observed Slurm standard-output log was about 109 MB and contained about
  1.6 million lines, mostly LightDock clustering messages. These messages are
  useful for debugging but are noisy for routine monitoring.
- Use `tail -n 50 <job>.out` or `grep` for engine milestones instead of opening
  the complete log interactively.
- A future logging improvement can write verbose LightDock output into a
  dedicated engine log file while retaining short progress messages in the
  Slurm log. Do not discard the detailed log without another retained copy.

## Restarting after interruption

Raw engine output that has already been written remains on disk after a Slurm
time limit or job failure. The batch results manifest records which pair-engine
runs completed successfully, allowing `ppinsight batch --resume` to skip those
runs on the next submission.

Current PPInsight writes the manifest after each completed engine run. Deploy
the current `batch_dock.py` to Hyak before relying on this behavior; older
versions wrote only after every ten results, which was too infrequent for a
six-run batch.

## HADDOCK runtime version

The conda environment pins its HADDOCK package version, but the automatic
container path can instead run the separate image
`ghcr.io/haddocking/haddock3:latest`. A container image is its own packaged
environment, so the conda pin does not control it.

For reproducible production results, record the exact runtime version used in
the job log and prefer either:

1. the tested `haddock3` executable in the activated Hyak environment, or
2. a container image pinned to a specific verified version or digest.

Avoid `latest` for benchmark runs that you expect to reproduce later.

## Interpreting partial runs

If docking completes but collection, quality evaluation, or plotting is
interrupted, the earlier output directories remain usable. Re-run only the
incomplete downstream command after confirming its inputs exist. A failed
Slurm job status means the shell script did not finish all of its steps; it
does not mean every earlier result was deleted.