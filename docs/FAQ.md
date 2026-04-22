# PPInsight FAQ

Frequently asked questions about PPInsight — a toolkit for comparing
protein-protein docking results across multiple prediction engines.

---

## General

### What does PPInsight do?

PPInsight collects docking scores from supported docking engines,
normalises them, and visualises the results with a unified Python API and
command-line interface.  It can also assess prediction quality using the
CAPRI protocol (requires DockQ) and predict binding affinity using PRODIGY.

### Which Python version is required?

Python ≥ 3.11.

### How do I install PPInsight?

```bash
# Core install
pip install ppinsight

# With DockQ quality assessment
pip install "ppinsight[quality]"

# With PRODIGY binding-affinity scoring
pip install "ppinsight[prodigy]"

# Everything
pip install "ppinsight[quality,prodigy,dev]"
```

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
ppinsight quality native.pdb pdbs/ -o quality.tsv
```

Then merge the quality scores into your unified scores file before calling
`quality_bar_chart`.

### What does the `quality_bar_chart` show?

Each bar represents one docking engine and is divided into four colour-coded
CAPRI tier segments, all scaled to 100 %.  The **N = X** label to the right
of each bar shows the total number of individual poses that were classified.

### DockQ is not installed — how do I get quality scores?

```bash
pip install "ppinsight[quality]"
# or:
pip install "dockq @ git+https://github.com/bjornwallner/DockQ.git"
```

---

## PRODIGY Binding Affinity

### What is PRODIGY?

PRODIGY (PROtein binDIng enerGY prediction) predicts the binding free energy
(ΔG, kcal/mol) and dissociation constant (Kd, M) of a protein-protein
complex from its 3-D structure.

Reference: Vangone A & Bonvin AMJJ. (2015). *eLife* 4:e07454.
DOI: 10.7554/eLife.07454

### How do I install PRODIGY?

```bash
pip install "ppinsight[prodigy]"
# or:
pip install prodigy-prot
```

### How do I score docked poses with PRODIGY?

```python
from ppinsight.prodigy import score_pdb, add_prodigy_to_scores
import pandas as pd

# Score a single PDB
result = score_pdb("complex.pdb", chains=["A", "B"])
print(result["prodigy_ddg"])   # ΔG in kcal/mol
print(result["prodigy_kd"])    # Kd in M

# Append to a unified scores DataFrame
scores = pd.read_csv("all_scores.tsv", sep="\t")
scores = add_prodigy_to_scores(scores, pdb_dir="pdbs/")
```

Or use the command-line tool:

```bash
ppinsight prodigy all_scores.tsv --pdb-dir pdbs/ --output scored.tsv
ppinsight prodigy all_scores.tsv --pdb-dir pdbs/ --top-n 5 --engine HADDOCK
```

### `score_pdb` returns `{"error": "no_contacts", …}` with NaN values

PRODIGY could not find inter-chain contacts in the PDB file.  Common causes:

1. **Single-chain structure** — PRODIGY requires at least two protein chains
   in the same file.  Make sure receptor and ligand chains are in one PDB.
2. **Chains too far apart** — the chains do not have inter-chain contacts
   within the 5.5 Å distance cutoff.  The structure may be a docking decoy
   where the two proteins are not in contact.
3. **Wrong chain IDs** — pass `chains=["A", "B"]` explicitly.

### `prodigy_ddg` vs `prodigy_kd` — which should I use for ranking?

For ranking docked poses, prefer `prodigy_ddg` (lower = stronger predicted
binding).  `prodigy_kd` is useful when you need to express the result in
concentration units for biological interpretation.

Both metrics have `higher_is_better=False` in `METRIC_METADATA`.

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
x-axis using `bbox_to_anchor=(0.5, -0.28)`.  If you are using a custom
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
| `ImportError: prodigy-prot is required` | Optional package missing | `pip install ppinsight[prodigy]` |
| `ValueError: No data found for metric 'X'` | Wrong metric name or empty DataFrame | Check `score_type` values with `df["score_type"].unique()` |
| `ValueError: No shared protein pairs` | Two models have no pairs in common | Check proteinA/proteinB columns match |
| `ModuleNotFoundError: No module named 'dockq'` | DockQ not installed | `pip install ppinsight[quality]` |
