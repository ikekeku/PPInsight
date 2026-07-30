# PPInsight End-to-End Tutorial

This tutorial is split into two honest paths:

1. **Real pipeline map** — fetch → dock → collect → compare.
2. **Prebuilt tutorial-data walkthrough** — explore the precomputed
    tutorial scores file in `tutorials/all_scores.tsv`.

By the end you will have:

1. A clear map of the real PPInsight pipeline and its intermediate files
2. A rebuilt tutorial scores file from repo example outputs
3. Verified compare commands that actually work on the shipped tutorial data
4. CAPRI quality summaries from the shipped example outputs

---

## Prerequisites

```bash
# Install ppinsight (from the repo root)
pip install -e .

# Optional: DockQ for quality evaluation
pip install -e ".[quality]"
```

For the prebuilt tutorial-data walkthrough, you do **not** need HADDOCK,
LightDock, or Rosetta installed because the repo already includes example
outputs under `examples/`.

For an end-to-end run, install the engine runtimes you plan to
use. In this environment, `lightdock3_setup.py` and `lightdock3.py` are
available, but HADDOCK is not, so the examples below distinguish
verified commands from the intended full workflow.

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

Verified here with real UniProt accessions:

```bash
ppinsight fetch P69905 P68871 \
    --fasta tutorials/fetch_demo.fasta \
    --csv tutorials/fetch_demo.csv \
    --pdb-dir data/input
```

That command was run in this environment successfully.  It wrote:

- `tutorials/fetch_demo.fasta`
- `tutorials/fetch_demo.csv`
- `data/input/P69905.pdb`
- `data/input/P68871.pdb`

This step is only to show the fetch-to-pairs handoff. The shipped
docking walkthrough switches back to the bundled `2UUY` example below.

The fetch step now writes pair-ready accession aliases, so the same
accessions can be used directly in `proteinA` and `proteinB`.

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

Verified here as a dry-run against the shipped example input PDBs:

```bash
ppinsight batch tutorials/demo_pairs.tsv \
    --engines lightdock \
    --pdb-dir examples/ppinsight_data/input_files \
    --dry-run \
    -o data/output/scores/batch_results.csv
```

In this environment the dry-run completed successfully for the pair
`2UUY_rec` vs `2UUY_lig`.

That step writes a batch results table plus engine-specific run
directories. It does not create the unified scores file; `collect` is
still the next step.

For a real run, remove `--dry-run`.  If you have all engine runtimes
installed, the intended form is:

```bash
ppinsight batch data/input/pairs/pairs.csv \
    --engines lightdock haddock rosetta \
    --pdb-dir data/input/ \
    -o data/output/scores/batch_results.csv
```

Batch-mode defaults used when you do not pass extra engine flags:

| Engine | Defaults in `ppinsight batch` |
|---|---|
| LightDock | `steps=10`, `cores=1`, `ANM=enabled`, swarms auto, glowworms auto, scoring uses the LightDock default |
| HADDOCK | `ncores` follows `--cores` (batch stages HADDOCK runs by default) |
| Rosetta | Enabled by default in batch mode (requires PyRosetta); defaults are `n_runs=10`, `top_n=20`, `relax=enabled`, `cluster=enabled`, `cluster_top_n=200`, `rmsd_cutoff=4.0`, `auto_filter=enabled`. |

Useful batch controls:

- `--cores N` sets CPU cores for supported runners.
- `--lightdock-no-anm` runs LightDock rigid-body from the start.
- `--lightdock-auto-clean-pdb` retries unsupported-residue failures with protein-only cleaned copies.

On known LightDock ANM setup atom-mismatch failures, PPInsight prompts for
confirmation before retrying that pair with ANM disabled.

### Step 4 — Collect scores from run directories

`collect` works on actual engine run directories. For a single known
pair, pass `--pair`. When you are annotating rows from a parsed pairs
file, use `--label-file` (legacy alias: `--pairs`) as needed.

Verified here with the shipped LightDock example output:

```bash
ppinsight collect examples/lightdock/simulation \
    --pair 2UUY_rec:2UUY_lig \
    -o /tmp/ppinsight_lightdock_collect.tsv
```

That command completed successfully in this environment and wrote 200
rows plus a provenance sidecar.

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
ppinsight compare data/output/scores/scores.tsv --plot-type quality_bar
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

## Step 11 — Ridge plot

```bash
ppinsight compare tutorials/all_scores.tsv --metric dockq \
    --plot-type ridge -o tutorials/sample_output/ridge_dockq.png
```

## Step 12 — CAPRI quality classification

```bash
# Print CAPRI quality tier counts
ppinsight compare tutorials/all_scores.tsv --capri-quality

# Quality bar chart (stacked: incorrect / acceptable / medium / high)
ppinsight compare tutorials/all_scores.tsv --plot-type quality_bar \
    -o tutorials/sample_output/quality_bar.png
```

The quality bar is **pose-level**.  The **N = X** label beside each bar
is the number of classified poses, not the number of protein pairs.


## Step 13 — CDF

```bash
ppinsight compare tutorials/all_scores.tsv --metric dockq \
    --plot-type cdf -o tutorials/sample_output/cdf_dockq.png
```

For `dockq`, vertical CAPRI threshold lines are drawn automatically.


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
- `ridge`
- `quality_bar`
- `cdf`
- tabular summaries and ranking tables

It does **not** support:

- `roc` — because `tutorials/all_scores.tsv` has no `label` column
- `scatter` — because the models do not share the same pairs
- `difference` — for the same alignment reason

Scatter and difference become valid only after a real aligned benchmark.
Reference commands once you have aligned multi-model coverage:

```bash
ppinsight compare data/output/scores/aligned_scores.tsv --metric dockq \
    --plot-type scatter --models haddock lightdock \
    -o tutorials/sample_output/scatter_haddock_vs_lightdock.png

ppinsight compare data/output/scores/aligned_scores.tsv --metric dockq \
    --plot-type difference --models haddock lightdock \
    -o tutorials/sample_output/difference_haddock_vs_lightdock.png
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
`non-interaction`), `--classify` produces a confusion-matrix summary:

```bash
ppinsight compare data/output/scores/scores.tsv --metric dockq --classify
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

## Condensed pipeline reference

If you have your own prepared PDB files and all required engine runtimes
installed, the end-to-end PPInsight sequence is:

```bash
# 1. Fetch or prepare structures
ppinsight fetch P69905 P68871 --pdb-dir data/input

# 2. Create a flat pairs file
ppinsight parse data/input/pairs/my_table.tsv \
    -o data/input/pairs/pairs.csv --stats

# 3. Run one engine for one pair (recommended first)
ppinsight lightdock PAIR_A PAIR_B --input-dir data/input/

# 3b. Optional: run many pairs/engines with batch mode
ppinsight batch data/input/pairs/pairs.csv \
    --engines lightdock haddock rosetta \
    --pdb-dir data/input/ \
    -o data/output/scores/batch_results.csv

# 4. Collect from pair-level run directories (single-engine or batch)
ppinsight collect data/output/lightdock_runs/PAIR_A_vs_PAIR_B/ \
    --pair PAIR_A:PAIR_B \
    -o data/output/scores/scores.tsv

# 5. Compare
ppinsight compare data/output/scores/scores.tsv --metric dockq
ppinsight compare data/output/scores/scores.tsv --capri-quality
ppinsight compare data/output/scores/scores.tsv --metric dockq --rank

# 6. (Optional) Evaluate against a native structure
ppinsight quality data/output/scores/scores.tsv native.pdb \
    -o data/output/scores/scores_quality.tsv
ppinsight quality data/output/lightdock_runs/PAIR_A_vs_PAIR_B/ native.pdb \
    --engine lightdock \
    -o quality.tsv
```

> See `data/README.md` for details on the `data/` directory layout.

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `--metric X not found` | Run `--list-metrics` to see available metrics |
| `--plot-type difference` or `scatter` errors | Both models must share the same pairs on the same metric |
| Scatter errors about alignment | Both `--models` must have scored the same pairs |
| `ModuleNotFoundError: DockQ` | Install with `pip install -e ".[quality]"` |
| `lightdock` not found | Verify the CLI tools are on `$PATH`, for example: `which lightdock3_setup.py` |
