# PPInsight End-to-End Tutorial

This tutorial walks through the **complete PPInsight workflow** using
real docking data from three engines: **HADDOCK**, **LightDock**, and
**Rosetta**.  By the end you will have:

1. A unified scores file combining all three engines
2. Tabular summaries comparing engine performance
3. Plots (violin, box, bar, heatmap, quality bar) saved as PNGs
4. CAPRI quality classification of docked models

---

## Prerequisites

```bash
# Install ppinsight (from the repo root)
pip install -e .

# Optional: DockQ for quality evaluation
pip install -e ".[quality]"
```

> You do **not** need HADDOCK, LightDock, or Rosetta installed.
> This tutorial uses pre-computed docking outputs that ship with the repo
> under `examples/`.

---

## What's in this directory

| File / Folder | Purpose |
|---|---|
| `README.md` | This walkthrough |
| `all_scores.tsv` | Pre-built unified scores (11,830 rows, 3 engines, 3 pairs, 10 metrics) |
| `walkthrough_scores.tsv` | Raw LightDock + HADDOCK scores from the original walkthrough |
| `demo_pairs.tsv` | Pairs file with interaction labels (for ROC) |
| `build_demo_scores.py` | Script that assembles `all_scores.tsv` from repo example data |
| `sample_output/` | Example plot PNGs produced by the commands below |

---

## Step 0 — Understand the data

The three protein pairs used in this tutorial:

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

## Step 1 — Collect scores

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

## Step 2 — Explore what's in the scores file

```bash
# What metrics are available?
ppinsight compare tutorials/all_scores.tsv --list-metrics

# What pairs are available?
ppinsight compare tutorials/all_scores.tsv --list-pairs

# Data-aware guide: which plot types and flags work with this file?
ppinsight compare tutorials/all_scores.tsv --guide
```

---

## Step 3 — Tabular summary

Every `ppinsight compare` command prints a tabular summary before any
plot.  To see **only** the table (no plot), use `--table`:

```bash
ppinsight compare tutorials/all_scores.tsv --metric dockq --table
```

This prints per-engine statistics (count, mean, median, min, max) for
the chosen metric across all pairs.

---

## Step 4 — Violin plot (default)

```bash
ppinsight compare tutorials/all_scores.tsv --metric dockq
```

This shows the full score distribution for each engine, faceted by
protein pair.

---

## Step 5 — Box plot

```bash
ppinsight compare tutorials/all_scores.tsv --metric dockq \
    --plot-type box -o tutorials/sample_output/box_dockq.png
```

---

## Step 6 — Bar chart (engine comparison)

```bash
ppinsight compare tutorials/all_scores.tsv --metric dockq \
    --plot-type bar -o tutorials/sample_output/bar_dockq.png
```

The bar chart requires ≥ 2 models or ≥ 2 pairs.  If the data doesn't
meet this prerequisite, the CLI will error with a hint.

---

## Step 7 — Heatmap

```bash
ppinsight compare tutorials/all_scores.tsv --metric dockq \
    --plot-type heatmap -o tutorials/sample_output/heatmap_dockq.png
```

Requires ≥ 2 pairs with `proteinA`/`proteinB` columns.

---

## Step 8 — CAPRI quality classification

```bash
# Print CAPRI quality tier counts
ppinsight compare tutorials/all_scores.tsv --capri-quality

# Quality bar chart (stacked: high/medium/acceptable/incorrect)
ppinsight compare tutorials/all_scores.tsv --plot-type quality_bar \
    -o tutorials/sample_output/quality_bar.png
```

CAPRI thresholds follow Lensink et al. 2016
([DOI: 10.1002/prot.25007](https://doi.org/10.1002/prot.25007)) with
DockQ fallback from Basu & Wallner 2016
([DOI: 10.1371/journal.pone.0161879](https://doi.org/10.1371/journal.pone.0161879)).

---

## Step 9 — Filter by pair

```bash
ppinsight compare tutorials/all_scores.tsv --metric dockq \
    --pair 2UUY_rec:2UUY_lig -o tutorials/sample_output/violin_2uuy.png
```

---

## Step 10 — Scatter (model agreement)

Compare how two engines score the **same pairs** on the **same metric**:

```bash
ppinsight compare tutorials/all_scores.tsv --metric dockq \
    --plot-type scatter --models haddock lightdock \
    -o tutorials/sample_output/scatter_haddock_vs_lightdock.png
```

> Scatter requires both models to have scored overlapping pairs.
> If there's no overlap, the CLI will error with a clear message.

---

## Step 11 — Compare per-model score files

If you collected scores into separate files per engine:

```bash
ppinsight compare lightdock_scores.tsv haddock_scores.tsv \
    --metric dockq --names LightDock HADDOCK
```

---

## Step 12 — Rank engines per pair

`--rank` prints a per-pair ranking table showing which engine scored
best on each protein pair:

```bash
ppinsight compare tutorials/all_scores.tsv --metric dockq --rank
```

---

## Step 13 — Classify interactions

If your pairs file includes `label` values (`interaction` /
`non-interaction`), `--classify` produces a confusion-matrix summary:

```bash
ppinsight compare tutorials/all_scores.tsv --metric dockq --classify
```

> Override the decision threshold with `--threshold 0.5`.

---

## Step 14 — Check the run record (provenance)

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

## Full pipeline (from PDB structures)

If you have your own PDB files and want to run the complete pipeline:

```bash
# 1. Place PDBs in data/input/
cp receptor.pdb ligand.pdb data/input/

# 2. (Optional) Parse a protein interaction table into a pairs file
ppinsight parse data/input/pairs/my_table.tsv -o data/input/pairs/pairs.csv --stats

# 3. Dock (choose one or more engines)
ppinsight lightdock receptor ligand --steps 100
ppinsight haddock receptor ligand --run
ppinsight rosetta receptor ligand --n-runs 5

# 4. Collect scores (writes to data/output/scores/scores.tsv by default)
ppinsight collect data/output/lightdock_runs data/output/haddock_runs \
    data/output/rosetta_runs

# 5. Compare
ppinsight compare data/output/scores/scores.tsv --metric dockq
ppinsight compare data/output/scores/scores.tsv --capri-quality
ppinsight compare data/output/scores/scores.tsv --metric dockq --rank

# 6. (Optional) Evaluate against a native structure
ppinsight quality docked_model.pdb native.pdb
ppinsight quality docked_models_dir/ native.pdb --engine lightdock -o quality.tsv
```

> See `data/README.md` for details on the `data/` directory layout.

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `--metric X not found` | Run `--list-metrics` to see available metrics |
| Bar chart won't render | Needs ≥ 2 models or ≥ 2 pairs — try `--plot-type violin` |
| Scatter errors about alignment | Both `--models` must have scored the same pairs |
| `ModuleNotFoundError: DockQ` | Install with `pip install -e ".[quality]"` |
| `lightdock` not found | Included with ppinsight — run `pip install -e .` and verify: `which lightdock3_setup.py` |
