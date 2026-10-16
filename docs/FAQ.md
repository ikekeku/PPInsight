# PPInsight FAQ

Frequently asked questions about PPInsight — a toolkit for comparing
protein-protein docking results across multiple prediction engines.

---

## General

### What does PPInsight do?

PPInsight collects docking scores from supported docking engines,
normalises them, and visualises the results with a unified Python API and
command-line interface.  It can also assess prediction quality using the
CAPRI protocol (requires DockQ).

### Which Python version is required?

Python ≥ 3.11.

### How do I install PPInsight?

```bash
# Core install
pip install ppinsight

# With DockQ quality assessment
pip install "ppinsight[quality]"

# Everything
pip install "ppinsight[quality,dev]"
```

### Does `ppinsight fetch` accept gene names like `VEGFA`?

Yes. `ppinsight fetch` accepts:

- UniProt accessions (for example `P15692`)
- FASTA-style IDs (for example `sp|P15692|VEGFA_HUMAN`)
- common human gene/name inputs (for example `VEGFA` or `"VEGFA human"`)

Name-based inputs are resolved against reviewed human UniProt records
(`organism_id=9606`) before sequence/PDB download.

### What do `--fasta FILE` and `--csv FILE` write in `ppinsight fetch`?

- `--fasta FILE`: combined FASTA records returned by UniProt.
- `--csv FILE`: structured metadata columns:
   `ID`, `Name`, `Description`, `Sequence Length`, `Sequence`.

### How do I preview ambiguous fetch terms before downloading?

Use `--search`:

```bash
ppinsight fetch --search "Neuropilin-1 human"
```

This prints top reviewed-human UniProt matches (accession, entry name,
protein name, genes, organism, sequence length, evidence, annotation score)
and exits without downloading files.

### How do I avoid accidental overwrite of local fetched PDB files?

`ppinsight fetch` is safe-by-default: existing aliases in `--pdb-dir` are
kept and conflicting writes are skipped.

- Use `--force` to overwrite intentionally.
- Use `--remove ...` to delete mistaken local fetch outputs first.

### What does `--no-auto-filter` mean for HADDOCK or Rosetta?

Some accession-named PDB files are mixed co-complex depositions that contain
extra proteins beyond the accession in the filename. For HADDOCK and Rosetta,
PPInsight checks DBREF metadata and keeps only the chains mapped to the
requested accession by default.

Use `--no-auto-filter` only when you intentionally want to dock the full
deposited complex (all co-complex chains), not the accession-mapped subset.

Why this exists: accession-based runs are usually meant to evaluate the
specific protein partner named in your pair file. Auto-filtering keeps the
input aligned with that intent and avoids docking unrelated co-crystallized
chains by accident.

---

## CAPRI Quality Assessment

### What is the CAPRI protocol?

CAPRI (Critical Assessment of PRedicted Interactions) is a community-wide
blind assessment of protein-protein docking methods.  Each predicted complex
is classified into one of four quality tiers:

| Tier       | f(nat) | L-RMSD     | i-RMSD     |
|------------|--------|------------|------------|
| High       | ≥ 0.50 | ≤ 1.0 Å    | ≤ 1.0 Å    |
| Medium     | ≥ 0.30 | ≤ 5.0 Å    | ≤ 2.0 Å    |
| Acceptable | ≥ 0.10 | ≤ 10.0 Å   | ≤ 4.0 Å    |
| Incorrect  | —      | —          | —          |

For details see: Lensink MF *et al.*, *Proteins* 84(S1):323-348, 2016.
DOI: 10.1002/prot.25007

### How does PPInsight assign CAPRI tiers?

From each pose's **DockQ** only, using these thresholds (Basu & Wallner
2016, DOI: 10.1371/journal.pone.0161879):

| Tier       | DockQ              |
|------------|--------------------|
| High       | ≥ 0.80             |
| Medium     | 0.49 – 0.80        |
| Acceptable | 0.23 – 0.49        |
| Incorrect  | < 0.23             |

DockQ was designed to reproduce the f(nat) / RMSD classes above, so the two
agree for most poses; borderline poses can land in different tiers.  Using
DockQ alone keeps every plot on one scale.  DockQ measured against your
native (`quality_dockq`) is used when a pose has it; HADDOCK's own `dockq`
only when it does not.

**Non-negotiable:** `quality_bar_chart` and `capri_summary_table` require
DockQ rows.  Passing engine-only scores (e.g. HADDOCK energy, Rosetta
total_score) raises a `ValueError`.  This is intentional — CAPRI tiers are
only meaningful when evaluated against a reference structure.
`quality_bar_chart` also shows one protein pair at a time; use `--pair` when
the file has several.

### `quality_bar_chart` raises `ValueError: CAPRI tiers are assigned from DockQ`

You are passing energy scores (score, total_score, luciferin_score) rather
than DockQ.  Run DockQ against your native structure:

```bash
ppinsight quality all_scores.tsv native.pdb -o all_scores_quality.tsv
```

This reads the collected pose paths from `output_path` and appends
DockQ/CAPRI evaluation rows directly into a new unified scores file that
`compare` can consume.

### What does the `quality_bar_chart` show?

The percentage of each engine's poses that falls into each CAPRI tier, for
one protein pair.  Each bar is one engine, divided into four colour-coded
tier segments that add up to 100 %.  The legend gives each tier's DockQ
range, and the **n = X** label under each engine name is the number of
poses that were classified.

### DockQ is not installed — how do I get quality scores?

```bash
pip install "ppinsight[quality]"
# or:
pip install "dockq @ git+https://github.com/nrontsis/DockQ.git@update-to-numpy>2"
```

---

## Visualisation

### What plot types are available?

| `--plot-type`  | Description                                    | Requires              |
|----------------|------------------------------------------------|-----------------------|
| `violin`       | Full score distribution + data points          | any metric            |
| `strip`        | Every pose as a point, with a box summary      | any metric            |
| `roc`          | ROC curves with AUC                            | `label` column        |
| `scatter`      | Model-agreement scatter (use with `--models`)  | 2 models, same pairs  |
| `quality_bar`  | % of poses in each CAPRI tier (from DockQ)     | DockQ, one pair       |
| `cdf`          | Empirical CDF per engine                       | any metric            |
| `funnel`       | Rank audit: each engine's own score vs DockQ for the same poses (no `--metric`) | native-referenced DockQ rows from `ppinsight quality`, one pair |
| `topn`         | Rank audit: best DockQ within each engine's top N poses, for every N (no `--metric`) | same as `funnel` |

### What is the difference between `dockq` and `quality_dockq`?

`dockq` (and `fnat`, `irmsd`, `lrmsd`) are **HADDOCK's own** values from its
caprieval step.  HADDOCK measures them against the reference it was given.
PPInsight runs HADDOCK without one, so caprieval compares each model with
HADDOCK's own best-scoring model — that model always shows DockQ = 1.0.
They say how similar a pose is to HADDOCK's favourite pose, not how close it
is to the real complex.  Plots label them "(HADDOCK caprieval)" and the CLI
prints a note whenever you use them.

`quality_dockq` (and `quality_fnat`, `quality_irmsd`, `quality_lrmsd`) are
measured by `ppinsight quality` against the native structure you supply,
the same way for every engine.  Use these whenever you have a native:

```bash
ppinsight quality data/output/scores/scores.tsv data/input/native.pdb \
    -o data/output/scores/scores_with_dockq.tsv
ppinsight compare data/output/scores/scores_with_dockq.tsv --metric quality_dockq
```

When a file has both, CAPRI tiers (`quality_bar`, `--capri-quality`) use the
native-referenced value for each pose and fall back to HADDOCK's own value
only when no native-referenced value exists; the legend says which was used.

### What does the top-N plot show, and how do I read it?

For each engine, it ranks the poses by the engine's own score and plots the
best DockQ found among the top 1, 2, 3, … N poses.  It answers: "if I only
looked at this engine's N best-scored poses, how good is the best one?"
It uses the same input as `funnel`.

Read it left to right:

- **High at N = 1** — the engine's own top pick is already good; its
  scoring ranks a near-native pose first.  This is the result you want.
- **Low at small N, then a jump at large N** — the engine *did* produce a
  good pose, but its scoring ranked it low.  The N where the jump happens is
  that pose's rank, so you would have had to look that far down the list.
- **Flat below the acceptable line all the way** — the engine never
  produced an acceptable pose, so better ranking would not help.

Because every engine is compared at the same N, an engine that outputs
thousands of extra poor poses is not penalised the way pose-averaged numbers
(means, percentages) penalise it.  A curve only rises or stays flat: adding
poses can never make the best one worse.

### Why does the ROC plot refuse DockQ?

ROC asks whether a score separates pairs that bind from pairs that do not.
DockQ, Fnat and the RMSDs compare a pose with a native complex, and a pair
that does not bind has no native complex, so those values cannot exist for
half of the comparison.  ROC and `--classify` therefore need a score every
pair has — an engine's own score (`score`, `luciferin_score`, `i_sc`).  Each
of those belongs to one engine, so each ROC covers one engine.

### When is the ROC plot useful?

When you have a labelled screen: many pairs marked `interaction` /
`non-interaction` (for example from `ppinsight parse` on an interactome
table).  Docking scores are generally weak at separating binders from
non-binders, so an area under the curve (AUC) near 0.5 is a common, honest
result rather than a bug — measuring that is the point of the benchmark.

### When should I use `strip` instead of `violin`?

Use `strip` when you want to see every pose, which works well up to a few
hundred poses per engine.  Nothing is smoothed, so a handful of near-native
poses stays visible.  Use `violin` for large pose ensembles, where
thousands of points would merge into a solid block.

Violin outlines are smoothed with 0.75 × seaborn's default bandwidth: the
default blurred the gap between two clusters of poses (e.g. many incorrect
poses plus a near-native group), and half the default looked too bumpy.
Above 5,000 poses per engine, both plots draw a random sample of 5,000 dots
(the outline, box and n still use every pose) and always draw the rare good
poses on top; the legend says which.

### What is the funnel plot for, and why does it not take `--metric`?

The funnel plot answers one question: does an engine's own score put its
near-native poses at the top?  Each engine gets one panel.  Each point is
one pose, with its DockQ against the native structure on the x-axis and
the engine's own ranking score on the y-axis (better scores point down).
Scores and DockQ are matched by `pose_id`, so every point is a real pose.
Because the axes are fixed to that question, the plot takes no `--metric`.

To make the input, run `ppinsight quality` on the collected scores file:

```bash
ppinsight quality data/output/scores/scores.tsv data/input/native.pdb \
    -o data/output/scores/scores_with_dockq.tsv
ppinsight compare data/output/scores/scores_with_dockq.tsv --plot-type funnel
```

If DockQ was computed on a separate folder of copied or renamed pose files,
the DockQ rows get different `pose_id` values and the funnel plot stops
with an error instead of guessing which poses match.

### The legend overlaps the x-axis title in `quality_bar_chart`

This was fixed in the current release.  The legend is placed below the
x-axis title as a figure-level legend.  If you are using a custom figure
size, you can move it:

```python
fig = quality_bar_chart(df)
fig.legends[0].set_bbox_to_anchor((0.5, -0.05))
```

### How do I add a CAPRI threshold line to a CDF plot?

CAPRI threshold lines are drawn automatically when the metric is on the
DockQ scale (`metric="dockq"` or `metric="quality_dockq"`):

```python
from ppinsight.visualizer import cdf_plot
fig = cdf_plot(df, metric="dockq")
```

For other metrics, you can add them manually:

```python
ax = fig.axes[0]
ax.axvline(0.5, color="red", ls="--", label="Custom threshold")
ax.legend()
```

---

## Troubleshooting

| Error | Likely cause | Fix |
|-------|-------------|-----|
| `ValueError: CAPRI tiers are assigned from DockQ` | Passing engine scores to a CAPRI-only function | Add DockQ with `ppinsight quality` first |
| `ValueError: No data found for metric 'X'` | Wrong metric name or empty DataFrame | Check `score_type` values with `df["score_type"].unique()` |
| `ValueError: No shared protein pairs` | Two models have no pairs in common | Check proteinA/proteinB columns match |
| `ModuleNotFoundError: No module named 'dockq'` | DockQ not installed | `pip install ppinsight[quality]` |
