# Scoring Evaluation & Cross-Model Comparison Framework

> **Status**: Design document — directly usable for implementation.
> **Depends on**: `METRIC_METADATA` (visualizer.py), `quality.py`, `normalize_scores()`

---

## 1. Scoring Metrics Per Engine

### 1.1 LightDock

| Aspect | Detail |
|---|---|
| **Primary scoring metric** | `luciferin_score` |
| **What PPInsight reads** | `luciferin_score` (from `rank_by_scoring.list` or swarm GSO output) |
| **Underlying function** | DFIRE statistical potential by default (`fastdfire` C implementation). 12+ alternative scoring functions available (DFIRE2, CPyDock, PISA, etc.) |
| **What it represents** | The "luciferin" value optimized by the Glowworm Swarm Optimization (GSO) algorithm. In biological terms it encodes a statistical pair-wise distance-dependent energy derived from known protein structures. |
| **Direction** | **Higher is better** — glowworms move toward brighter (higher-scoring) neighbors. |
| **Units** | Arbitrary (DFIRE potential units). Not REU, not kcal/mol. |
| **Ranking workflow** | `lgd_rank.py` aggregates per-swarm top poses → `rank_by_scoring.list` sorted descending by luciferin. |

**Source**: LightDock documentation — *Basics* and *Theory* pages (lightdock.org); Jiménez-García et al. (2018).

### 1.2 HADDOCK

| Aspect | Detail |
|---|---|
| **Primary scoring metric** | `score` (the HADDOCK score) |
| **What PPInsight reads** | `score`, plus CAPRI metrics (`dockq`, `irmsd`, `lrmsd`, `fnat`) from `capri_ss.tsv` |
| **Underlying function** | Linear combination of physics-based energy terms: |
| | $HS = w_{vdw} E_{vdw} + w_{elec} E_{elec} + w_{desolv} E_{desolv} + w_{bsa} BSA + w_{air} E_{air} + w_{sym} E_{sym}$ |
| **Weights vary by stage** | **rigidbody**: $0.01 E_{vdw} + 1.0 E_{elec} + 1.0 E_{desolv} - 0.01 BSA + 0.01 E_{air}$ |
| | **flexref**: $1.0 E_{vdw} + 1.0 E_{elec} + 1.0 E_{desolv} - 0.01 BSA + 0.1 E_{air}$ |
| | **emscoring / mdscoring**: $1.0 E_{vdw} + 0.2 E_{elec} + 1.0 E_{desolv} + 0.0 E_{air}$ |
| **Direction** | **Lower is better** — it is an energy-like score. |
| **Units** | Arbitrary weighted energy units (not directly kcal/mol). |
| **CAPRI metrics** | HADDOCK's `caprieval` module independently computes Fnat, I-RMSD, L-RMSD, DockQ, Global-RMSD against a reference. These are the *quality* metrics. |

**Source**: HADDOCK3 User Manual — *HADDOCK scoring function* (bonvinlab.org/haddock3-user-manual/haddocking.html); Dominguez, Boelens & Bonvin (2003).

**Critical caveat from HADDOCK documentation**: *"Do not take scores as proxies of binding affinity to compare different complexes — compare scores only within the same system/complex."* (Best Practice Guide → Clustering → Dos and Don'ts).

### 1.3 Rosetta (RosettaDock)

| Aspect | Detail |
|---|---|
| **Primary scoring metric** | `I_sc` (interface score) |
| **What PPInsight reads** | `i_sc` (or `dg_separated`), `total_score`, `irms`, `rms`, `cluster_size` |
| **Underlying function** | $I_{sc} = E_{total}^{complex} - E_{total}^{partner_A} - E_{total}^{partner_B}$ |
| | The total score itself is the Rosetta all-atom energy (ref2015 score function) comprising vdW, solvation, H-bond, electrostatic, rotamer, and pair-wise statistical terms. |
| **What it represents** | Interface score (`I_sc`) isolates the energetic contribution of the protein–protein interface. It is computed by scoring the complex, then rigidly separating the partners and scoring each in isolation, then subtracting. |
| **Direction** | **Lower is better** — more negative = stronger predicted binding. |
| **Units** | REU (Rosetta Energy Units). Not kcal/mol. |
| **Typical range** | RosettaDock documentation states: *"Typical values for I_sc of good decoys are in the range of −5 to −10."* |

**Source**: RosettaDock Protocol documentation (rosettacommons.org); Gray et al. (2003); Marze et al. (2018).

### 1.4 Cross-Engine Quality Metrics (DockQ / CAPRI)

All three engines can be evaluated post-hoc via the **DockQ** tool, which is engine-agnostic and requires only a predicted complex PDB + a native reference PDB.

| Metric | Direction | What it measures |
|---|---|---|
| `quality_dockq` | Higher is better (0–1) | Combined quality score |
| `quality_fnat` | Higher is better (0–1) | Fraction of native contacts |
| `quality_irmsd` | Lower is better (Å) | Interface RMSD |
| `quality_lrmsd` | Lower is better (Å) | Ligand RMSD |

PPInsight already implements this via `ppinsight.quality.evaluate_complex()`.

---

## 2. Quality Thresholds

### 2.1 DockQ/CAPRI — Documentation-Defined Thresholds ✅

These are the only universally accepted, publication-defined thresholds across docking engines:

| Category | DockQ threshold | Fnat | I-RMSD (Å) | L-RMSD (Å) |
|---|---|---|---|---|
| **High** | ≥ 0.80 | ≥ 0.50 | ≤ 1.0 | ≤ 1.0 |
| **Medium** | ≥ 0.49 | ≥ 0.30 | ≤ 2.0 | ≤ 5.0 |
| **Acceptable** | ≥ 0.23 | ≥ 0.10 | ≤ 4.0 | ≤ 10.0 |
| **Incorrect** | < 0.23 | < 0.10 | > 4.0 | > 10.0 |

**Source**: Basu & Wallner (2016), *PLoS ONE* 11(8), e0161879. Already implemented in PPInsight at `quality.py:CAPRI_THRESHOLDS`.

The full CAPRI protocol uses Fnat + L-RMSD + I-RMSD jointly (as in `visualizer.capri_quality()`). DockQ is a single scalar that approximates this multi-metric classification.

### 2.2 Engine-Native Scores — No Universal Thresholds ❌

**LightDock `luciferin_score`**:
> No documentation-defined quality thresholds exist. The LightDock documentation and publications do not define threshold values that classify a luciferin score as "good" or "bad" in absolute terms. Ranking is purely relative: top-N models sorted by descending score. The magnitude depends on the scoring function chosen (DFIRE, CPyDock, etc.) and the specific protein system.

**HADDOCK `score`**:
> No universal thresholds exist. The HADDOCK documentation explicitly warns: *"Do not take scores as proxies of binding affinity to compare different complexes."* The score is only meaningful for ranking models within a single docking run for a single system. The HADDOCK weights change per protocol stage, and optimal score ranges are system-dependent. HADDOCK's own quality assessment relies on CAPRI metrics (Fnat, I-RMSD, L-RMSD, DockQ) computed by the `caprieval` module — not on the HADDOCK score itself.

**Rosetta `I_sc`**:
> No formal thresholds are published. The RosettaDock documentation provides only a rough empirical guideline: *"Typical values for I_sc of good decoys are in the range of −5 to −10."* This is not a formal classification system and varies with system size, scoring function version (ref2015 vs. earlier), and prepacking quality. Rosetta scores are in REU (not kcal/mol) and are not calibrated for absolute binding affinity prediction.

### 2.3 Alternative Scoring Metrics That *Do* Offer Thresholds

Since engine-native scores lack thresholds, the only path to universal binning is through **DockQ/CAPRI metrics**, which are already implemented in PPInsight. For users who have a native reference structure, this is the recommended path.

For users who **lack** a native reference structure (the common real-world case), the scoring metrics can only be used for **within-run relative ranking** — there is no scientifically defensible way to assign absolute quality labels.

---

## 3. Cross-Model Normalization

### 3.1 The Problem

Engine-native scores are on incomparable scales:
- LightDock luciferin: positive values, higher = better, DFIRE units
- HADDOCK score: negative values, lower = better, weighted energy units
- Rosetta I_sc: negative values, lower = better, REU

Comparing raw values across engines is meaningless.

### 3.2 Recommended Approach: Proportional (0–100%) Quality Classification

**When a native reference is available**: Use DockQ/CAPRI classification, then compare at the **classification level**:

$$\text{Proportion}_{\text{category}} = \frac{N_{\text{models in category}}}{N_{\text{total models}}} \times 100\%$$

For each engine, compute the percentage of models in each CAPRI bin (high / medium / acceptable / incorrect). This allows direct cross-engine comparison as proportions.

**Example**: "LightDock produced 15% high-quality models vs. HADDOCK's 22% vs. Rosetta's 18%."

### 3.3 When No Native Reference Is Available

Use PPInsight's existing `normalize_scores(method='rank')` to convert engine-native scores to percentile ranks (0–1):

$$\text{rank}_i = \frac{\text{rank}(x_i)}{N}$$

This preserves ordering within each engine and places all scores on a common [0, 1] scale. Combined with the direction-aware flipping already in `normalize_scores()`, lower-is-better metrics are flipped so that 1.0 always means "best within this run."

**Important**: Rank normalization enables visual comparison of *distributions* but does **not** enable quality classification. A model at the 90th percentile in a poor run is still a poor model.

### 3.4 Summary of Normalization Options

| Scenario | Method | Output scale | Enables quality labels? |
|---|---|---|---|
| Native available | DockQ → CAPRI bin | 4 categories | ✅ Yes (high/med/acc/inc) |
| Native available | DockQ → proportion per bin | 0–100% | ✅ Yes |
| No native | `normalize_scores(method='rank')` | 0.0–1.0 | ❌ No — only relative |
| No native | `normalize_scores(method='minmax')` | 0.0–1.0 | ❌ No — only relative |

---

## 4. Visualization Design

### 4.1 Visualization A — Normalized Stacked Bar Chart (CAPRI Proportions)

**Purpose**: Compare the distribution of model quality across engines at a glance.

**Specification**:

| Element | Detail |
|---|---|
| **x-axis** | Engine name (one bar per engine) |
| **y-axis** | Proportion (0–100%) |
| **Segments** | 4 stacked segments per bar: high (green), medium (gold), acceptable (orange), incorrect (red) |
| **Annotation** | Total N (number of models) printed above each bar |
| **Data source** | DockQ scores → `classify_capri()` → count per category → proportion |
| **Requires** | Native reference structure |

**Implementation notes**:
- This is an extension of the existing `quality_bar_chart()` in `visualizer.py`.
- The existing function already plots CAPRI proportions per model. The new version groups by engine when multiple engines are present in a combined DataFrame.
- Color palette: use the existing `CAPRI_COLORS` or define a consistent 4-color mapping.

**Pseudocode sketch**:
```
for each engine:
    counts = df[df.model == engine].capri_quality.value_counts()
    proportions = counts / counts.sum()
    plot stacked bar
annotate total N above each bar
```

### 4.2 Visualization B — Cumulative Distribution Function (CDF) Plot

**Purpose**: Compare how scores (or DockQ values) distribute across engines. The CDF reveals not just averages but the full shape — are there a few excellent models, or a broad plateau of mediocre ones?

**Specification**:

| Element | Detail |
|---|---|
| **x-axis** | DockQ score (0–1) — or normalized engine-native score |
| **y-axis** | Cumulative proportion (0–1) |
| **Lines** | One line per engine, different color |
| **Reference lines** | Vertical dashed lines at DockQ = 0.23, 0.49, 0.80 (CAPRI thresholds) |
| **Legend** | Engine name + total N |
| **Requires** | DockQ scores (with native), or rank-normalized scores (without native) |

**Implementation notes**:
- Use `numpy.sort()` + `numpy.arange()` / `len` for empirical CDF.
- Or use `seaborn.ecdfplot()` for a one-liner.
- When plotting DockQ, the CAPRI threshold lines provide interpretive anchors.
- When plotting rank-normalized scores (no native), omit the threshold lines since the axis is not DockQ.

**Interpretation**: A CDF curve that rises steeply early (left side) means most models are low quality. A curve that stays low until the right side means the engine found many high-quality models. The engine whose CDF curve is furthest to the right "wins."

---

## 5. Key Interpretation Notes

### 5.1 Engine Scores ≠ Binding Affinity

None of the three engine-native scores (luciferin, HADDOCK score, Rosetta I_sc) are calibrated predictions of binding affinity (ΔG). They are **internal ranking functions** optimized for each engine's search algorithm:
- LightDock luciferin drives the GSO swarm dynamics.
- HADDOCK score combines physics terms with restraint satisfaction.
- Rosetta I_sc is a difference of all-atom energies.

Comparing absolute values across engines is not meaningful.

### 5.2 Within-Run Ranking Is Valid; Cross-Run Comparison Is Not

All three engines are designed so that **better-scoring models within a single run for a single target** are more likely to be near-native. This is their intended use. Comparing an I_sc from one target to an I_sc from a different target, or comparing an I_sc to a HADDOCK score, is not supported by the scoring functions.

### 5.3 DockQ/CAPRI Is the Common Currency

The only scientifically defensible cross-engine comparison uses DockQ/CAPRI, which:
- Is computed independently from any engine's scoring function.
- Requires a known native structure (the "ground truth").
- Has published, well-accepted thresholds (Basu & Wallner 2016).

PPInsight already implements this via `quality.py` and `visualizer.capri_quality()`.

### 5.4 The "No Native" Case

In real-world use, users typically do **not** have a native structure (that's why they're docking!). In this case:
- Engine-native scores are used for relative ranking only.
- Rank normalization can place scores on a common scale for visual comparison.
- **No quality labels should be assigned.** The stacked bar chart (§4.1) is not applicable; only the CDF plot (§4.2) with rank-normalized scores is valid.

### 5.5 Score Directions Are Non-Negotiable

`METRIC_METADATA` in `visualizer.py` is the single source of truth for score directions. The `normalize_scores()` function uses `higher_is_better` to flip lower-is-better metrics before normalization. Any new metric must be registered in `METRIC_METADATA` with its correct direction.

---

## Summary Table

| Engine | Primary Metric | Direction | Has Universal Thresholds? | Cross-Engine Comparable? |
|---|---|---|---|---|
| LightDock | `luciferin_score` | Higher = better | ❌ No | ❌ Not directly |
| HADDOCK | `score` | Lower = better | ❌ No | ❌ Not directly |
| Rosetta | `i_sc` / `dg_separated` | Lower = better | ❌ No (rough guide only) | ❌ Not directly |
| All (via DockQ) | `quality_dockq` | Higher = better | ✅ Yes (CAPRI) | ✅ Yes |

**Bottom line**: For cross-engine comparison, use DockQ + CAPRI classification. For within-engine ranking, use each engine's native score. For visual comparison of distributions, use rank normalization + CDF plot.
