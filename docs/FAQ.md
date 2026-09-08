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

**Non-negotiable:** `quality_bar_chart` and `capri_summary_table` require
CAPRI metrics (fnat, irmsd, lrmsd, or dockq).  Passing engine-only scores
(e.g. HADDOCK energy, Rosetta total_score) will raise a `ValueError`.  This
is intentional — CAPRI tier classification is only meaningful when evaluated
against a native reference structure.

### `quality_bar_chart` raises `ValueError: quality_bar_chart requires CAPRI metrics`

You are passing energy scores (score, total_score, luciferin_score) rather
than quality metrics.  Run DockQ to obtain fnat/irmsd/lrmsd/dockq values:

```bash
ppinsight quality all_scores.tsv native.pdb -o all_scores_quality.tsv
```

This reads the collected pose paths from `output_path` and appends
DockQ/CAPRI evaluation rows directly into a new unified scores file that
`compare` can consume.

### What does the `quality_bar_chart` show?

Each bar represents one docking engine and is divided into four colour-coded
CAPRI tier segments, all scaled to 100 %.  The **N = X** label to the right
of each bar shows the total number of individual poses that were classified.

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
| `ridge`        | Ridgeline KDE per engine                       | any metric            |
| `roc`          | ROC curves with AUC                            | `label` column        |
| `scatter`      | Model-agreement scatter (use with `--models`)  | 2 models, same pairs  |
| `quality_bar`  | 100 % stacked pose-level CAPRI tier bar        | CAPRI metrics         |
| `cdf`          | Empirical CDF per engine                       | any metric            |
| `difference`   | Per-pair difference histogram (use `--models`) | 2 models, same pairs  |

### The legend overlaps the x-axis title in `quality_bar_chart`

This was fixed in the current release.  The legend is now placed below the
x-axis using `bbox_to_anchor=(0.5, -0.205)`.  If you are using a custom
figure size, you can adjust it:

```python
fig = quality_bar_chart(df)
ax = fig.axes[0]
ax.get_legend().set_bbox_to_anchor((0.5, -0.32))
```

### How do I add a CAPRI threshold line to a CDF plot?

CAPRI threshold lines are drawn automatically when `metric="dockq"`:

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
| `ValueError: quality_bar_chart requires CAPRI metrics` | Passing engine scores to a CAPRI-only function | Add fnat/irmsd/lrmsd/dockq columns first |
| `ValueError: No data found for metric 'X'` | Wrong metric name or empty DataFrame | Check `score_type` values with `df["score_type"].unique()` |
| `ValueError: No shared protein pairs` | Two models have no pairs in common | Check proteinA/proteinB columns match |
| `ModuleNotFoundError: No module named 'dockq'` | DockQ not installed | `pip install ppinsight[quality]` |
