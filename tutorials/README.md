# PPInsight End-to-End Tutorial

This tutorial has two parts:

1. **Part 1, the real pipeline:** fetch → dock → collect → compare.
2. **Part 2, prebuilt scores:** explore `tutorials/all_scores.tsv` without
   running any docking engine.

---

## Prerequisites

```bash
# Install ppinsight (from the repo root)
pip install -e .

# Optional: DockQ for quality evaluation
pip install -e ".[quality]"
```

Part 2 needs no docking engines: the repo includes example outputs under
`examples/`. For Part 1, install the engines you plan to use.

---

## What's in this directory

| File / Folder | Purpose |
|---|---|
| `README.md` | This walkthrough |
| `all_scores.tsv` | Pre-built unified scores assembled from real repo example outputs |
| `walkthrough_scores.tsv` | Raw LightDock + HADDOCK scores from the original walkthrough |
| `demo_pairs.tsv` | Pairs file with interaction labels (for ROC) |
| `build_demo_scores.py` | Script that assembles `all_scores.tsv` from repo example data |
| `sample_output/` | Example plot PNGs produced by the commands below |

---

## Part 1 — Real Pipeline Workflow

This section follows the actual PPInsight product flow:

`fetch -> dock (single-engine or batch) -> collect -> compare`

### Step 0 — Understand the handoff points

PPInsight keeps each stage separate on purpose:

- `ppinsight fetch` downloads sequences and available PDB structures.
- `ppinsight batch` or the single-engine commands run docking.
- `ppinsight collect` converts engine outputs into a unified scores file.
- `ppinsight compare` is the visualization and summary entry point.

Practical detail: **`batch` does not write the unified `scores.tsv` by
itself.** It writes a batch results table plus engine run directories,
and `ppinsight collect` is still the next step.

### Step 0.5 — Run quick preflight checks

```bash
# Which engines are currently available in this environment?
ppinsight batch --list-engines
```

### Step 1 — Fetch protein data

```bash
ppinsight fetch P69905 P68871 \
    --fasta tutorials/fetch_demo.fasta \
    --csv tutorials/fetch_demo.csv \
    --pdb-dir data/input
```

This writes the FASTA and CSV files plus `data/input/P69905.pdb` and
`data/input/P68871.pdb`. The files are named by accession, so you can use
the same accessions in `proteinA` and `proteinB`. Each file holds the first
PDB entry UniProt lists for that protein; check which one with
`grep -m1 '^HEADER' data/input/P69905.pdb`.

The docking steps below use the bundled `2UUY` example instead.

### Step 2 — Prepare identifiers for docking inputs

For single-engine runs, you can call the docking command directly with
two PDB stems (receptor and ligand). A pairs file is only required when
you use `ppinsight batch`.

Minimal pairs-file example (for optional batch mode):

```tsv
proteinA	proteinB	label
P69905	P68871	interaction
```

These names must match the PDB filename stems in the directory you pass
to `--pdb-dir`.

The dry-run below returns to the shipped `tutorials/demo_pairs.tsv`
example because those prepared `2UUY` structures are already bundled in
the repo.

### Step 3 — Start with a single engine run

Before launching batch mode, validate one engine and one pair so path
resolution and runtime dependencies are confirmed.

```bash
ppinsight lightdock 2UUY_rec 2UUY_lig \
    --input-dir examples/ppinsight_data/input_files \
    --steps 10 --skip-postprocess
```

When this works, move on to batch mode.

### Step 3.5 — Batch docking (optional once single-run is validated)

Preview a batch against the bundled example PDBs:

```bash
ppinsight batch tutorials/demo_pairs.tsv \
    --engines lightdock \
    --pdb-dir examples/ppinsight_data/input_files \
    --dry-run \
    -o data/output/scores/batch_results.csv
```

Remove `--dry-run` to run it. Batch writes a results table and engine run
directories, not a scores file; `collect` is still the next step. For batch
defaults and flags, see the README's batch section.

### Step 4 — Collect scores from run directories

`collect` works on actual engine run directories. For a single known
pair, pass `--pair`. When you are annotating rows from a parsed pairs
file, use `--label-file` (legacy alias: `--pairs`) as needed.

Try it on the bundled LightDock output:

```bash
ppinsight collect examples/lightdock/simulation \
    --pair 2UUY_rec:2UUY_lig \
    -o data/output/scores/lightdock_example.tsv
```

For your own pair-level run directories the pattern is:

```bash
ppinsight collect data/output/lightdock_runs/2UUY_rec_vs_2UUY_lig/ \
    --pair 2UUY_rec:2UUY_lig \
    -o data/output/scores/lightdock_2uuy.tsv
```

If more than one engine has scored the same pair, pass those pair-level
run directories together:

```bash
ppinsight collect \
    data/output/lightdock_runs/2UUY_rec_vs_2UUY_lig/ \
    data/output/haddock_runs/2UUY_rec_vs_2UUY_lig/ \
    --pair 2UUY_rec:2UUY_lig \
    -o data/output/scores/scores.tsv
```

### Step 5 — Compare

Once you have a unified scores file, `ppinsight compare` becomes the
single entry point for summaries, classification, and plots.

```bash
ppinsight compare data/output/scores/scores.tsv --metric dockq --table
ppinsight compare data/output/scores/scores.tsv --plot-type quality_bar --pair proteinA:proteinB
ppinsight compare data/output/scores/scores.tsv --guide
```

## Part 2 — Explore The Prebuilt Tutorial Scores

The rest of this walkthrough uses the shipped example outputs and the
prebuilt `tutorials/all_scores.tsv` file. This is the fastest way to
learn `compare` without running docking locally.

## Step 6 — Understand the shipped tutorial data

The three protein pairs represented in the shipped tutorial assets:

| Pair | Receptor | Ligand | Engine(s) |
|------|----------|--------|-----------|
| 2UUY | 2UUY\_rec | 2UUY\_lig | LightDock, HADDOCK |
| e2aP / hpr | e2aP\_1F3G | hpr\_ensemble | HADDOCK |
| COL\_D / IMM\_D | COL\_D | IMM\_D | Rosetta |

The raw docking outputs live in `examples/`:
- `examples/haddock3/run1-test/` — HADDOCK caprieval stages
- `examples/lightdock/simulation/` — LightDock swarm outputs
- `examples/rosetta/Protein-Protein-Docking/output_files/` — Rosetta score files

---

## Step 7 — Rebuild `all_scores.tsv` from the shipped example outputs

If you want to rebuild `all_scores.tsv` from scratch (or see how
collection works), run the build script:

```bash
cd tutorials/
python build_demo_scores.py
```

Or use `ppinsight collect` directly on individual engine outputs:

```bash
# Collect LightDock scores
ppinsight collect examples/lightdock/simulation \
    --pair 2UUY_rec:2UUY_lig -o lightdock_scores.tsv

# Collect HADDOCK scores
ppinsight collect examples/haddock3/run1-test \
    --pair e2aP_1F3G:hpr_ensemble -o haddock_scores.tsv
```

---

## Step 8 — Explore what's in the scores file

```bash
# What metrics are available?
ppinsight compare tutorials/all_scores.tsv --list-metrics

# What pairs are available?
ppinsight compare tutorials/all_scores.tsv --list-pairs

# Data-aware guide: which plot types and flags work with this file?
ppinsight compare tutorials/all_scores.tsv --guide
```

---

## Step 9 — Tabular summary

Every `ppinsight compare` command prints a tabular summary before any
plot.  To see **only** the table (no plot), use `--table`:

```bash
ppinsight compare tutorials/all_scores.tsv --metric dockq --table
```

This prints per-engine statistics (count, mean, median, min, max) for
the chosen metric across all pairs.

---

## Step 10 — Violin plot (default)

```bash
ppinsight compare tutorials/all_scores.tsv --metric dockq
```

This shows the full score distribution for each engine, faceted by
protein pair.

---

## Step 11 — Strip plot

```bash
ppinsight compare tutorials/all_scores.tsv --metric dockq \
    --plot-type strip -o tutorials/sample_output/strip_dockq.png
```

Every pose is drawn as a point, with a box showing the median and
interquartile range.  The number under each engine name is how many poses
that engine contributed.

## Step 12 — CAPRI quality classification

```bash
# Print CAPRI quality tier counts
ppinsight compare tutorials/all_scores.tsv --capri-quality

# Quality bar chart: % of poses in each CAPRI tier, one protein pair at a time
ppinsight compare tutorials/all_scores.tsv --plot-type quality_bar \
    --pair 2UUY_rec:2UUY_lig -o tutorials/sample_output/quality_bar.png
```

The quality bar is **pose-level** and covers one protein pair.  Tiers come
from each pose's DockQ, with the thresholds listed in the legend.  The
**n = X** label under each engine name is the number of classified poses.


## Step 13 — CDF

```bash
ppinsight compare tutorials/all_scores.tsv --metric dockq \
    --plot-type cdf -o tutorials/sample_output/cdf_dockq.png
```

For DockQ, vertical CAPRI threshold lines are drawn automatically.


## Step 14 — Filter by pair

```bash
ppinsight compare tutorials/all_scores.tsv --metric dockq \
    --pair 2UUY_rec:2UUY_lig -o tutorials/sample_output/violin_2uuy.png
```

---

## Step 15 — What the shipped tutorial data does **not** support

Run

```bash
ppinsight compare tutorials/all_scores.tsv --guide
```

to see what the shipped tutorial file supports:

- `violin`
- `strip`
- `quality_bar`
- `cdf`
- tabular summaries and ranking tables

It does **not** support:

- `roc` — because `tutorials/all_scores.tsv` has no `label` column
- `scatter` — because the models do not share the same pairs
- `funnel` and `topn` — because the file has no `pose_id` column and no
  DockQ rows measured against a native structure (`quality_dockq`)

The `dockq` metric used above is HADDOCK's own caprieval value, so plots
label it "DockQ (HADDOCK caprieval)".  It is measured against HADDOCK's own
best-scoring model, not against a native complex; see the FAQ entry on
`dockq` vs `quality_dockq`.

Scatter becomes valid only after a real aligned benchmark.
Reference command once you have aligned multi-model coverage:

```bash
ppinsight compare data/output/scores/aligned_scores.tsv --metric dockq \
    --plot-type scatter --models haddock lightdock \
    -o tutorials/sample_output/scatter_haddock_vs_lightdock.png
```

The funnel and top-N plots become valid once you run `ppinsight quality`
on a collected scores file for one target with a native structure (use
`--plot-type topn` in place of `funnel` for the top-N plot):

```bash
ppinsight quality data/output/scores/scores.tsv data/input/native.pdb \
    -o data/output/scores/scores_with_dockq.tsv
ppinsight compare data/output/scores/scores_with_dockq.tsv \
    --plot-type funnel -o tutorials/sample_output/funnel.png
```

---

## Step 16 — Compare per-model score files

If you collected scores into separate files per engine:

```bash
ppinsight compare lightdock_scores.tsv haddock_scores.tsv \
    --metric dockq --names LightDock HADDOCK
```

---

## Step 17 — Rank engines per pair

`--rank` prints a per-pair ranking table showing which engine scored
best on each protein pair:

```bash
ppinsight compare tutorials/all_scores.tsv --metric dockq --rank
```

---

## Step 18 — Classify interactions

If your collected scores file includes `label` values (`interaction` /
`non-interaction`), `--classify` produces a confusion-matrix summary. Use an
engine's own score; DockQ is refused because non-binding pairs have no native
structure:

```bash
ppinsight compare data/output/scores/scores.tsv --metric score --classify
```

This is not available on the shipped `tutorials/all_scores.tsv` because
that file has no `label` column.

> Override the decision threshold with `--threshold 0.5`.

---

## Step 19 — Check the run record (provenance)

When `ppinsight collect` writes a scores file it also writes a companion
**run record** — a JSON sidecar documenting where the numbers came from
(source directories, engine parameters, timestamps).

```bash
cat tutorials/all_scores.tsv.provenance.json | python -m json.tool
```

This is useful for auditing results or reproducing a collection run.

---

## Sample output

The `sample_output/` folder contains example PNGs produced by these
commands.  Use them as a reference for what to expect.

---

See `data/README.md` for the `data/` directory layout.

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `--metric X not found` | Run `--list-metrics` to see available metrics |
| `--plot-type scatter` errors | Both `--models` must have scored the same pairs on the same metric |
| `--plot-type funnel` says no pose_id matched | Run `ppinsight quality` on the collected scores file, not on a folder of copied or renamed poses |
| `ModuleNotFoundError: DockQ` | Install with `pip install -e ".[quality]"` |
| `lightdock` not found | Verify the CLI tools are on `$PATH`, for example: `which lightdock3_setup.py` |
