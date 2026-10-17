# Scoring Evaluation & Cross-Model Comparison Framework

What each engine's score means, which thresholds exist, and how PPInsight
compares engines.

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

`ppinsight quality` computes these.

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

**Source**: Basu & Wallner (2016), *PLoS ONE* 11(8), e0161879.

The full CAPRI protocol uses Fnat, L-RMSD and I-RMSD jointly; DockQ
approximates it with one number. PPInsight assigns tiers from DockQ alone.

### 2.2 Engine-Native Scores — No Universal Thresholds ❌

**LightDock `luciferin_score`**:
> No documentation-defined quality thresholds exist. The LightDock documentation and publications do not define threshold values that classify a luciferin score as "good" or "bad" in absolute terms. Ranking is purely relative: top-N models sorted by descending score. The magnitude depends on the scoring function chosen (DFIRE, CPyDock, etc.) and the specific protein system.

**HADDOCK `score`**:
> No universal thresholds exist. The HADDOCK documentation explicitly warns: *"Do not take scores as proxies of binding affinity to compare different complexes."* The score is only meaningful for ranking models within a single docking run for a single system. The HADDOCK weights change per protocol stage, and optimal score ranges are system-dependent. HADDOCK's own quality assessment relies on CAPRI metrics (Fnat, I-RMSD, L-RMSD, DockQ) computed by the `caprieval` module — not on the HADDOCK score itself.

**Rosetta `I_sc`**:
> No formal thresholds are published. The RosettaDock documentation provides only a rough empirical guideline: *"Typical values for I_sc of good decoys are in the range of −5 to −10."* This is not a formal classification system and varies with system size, scoring function version (ref2015 vs. earlier), and prepacking quality. Rosetta scores are in REU (not kcal/mol) and are not calibrated for absolute binding affinity prediction.

### 2.3 With and Without a Native Structure

With a native structure, use DockQ/CAPRI to label pose quality. Without one,
engine scores support only **within-run ranking**; there is no defensible way
to assign quality labels.

### 2.4 Native-Reference Rank Audit

When a native reference is available, evaluate all score-traceable poses (or a
documented broad subset) with DockQ and summarize each engine separately:

| Field | Question answered |
|---|---|
| `poses_ranked` | How many poses did the engine score? |
| `poses_evaluated` | How many of those poses were available for DockQ? |
| `best_dockq` | Did the engine generate a near-native pose at all? |
| `engine_rank_of_best_dockq` | Where did that best-DockQ pose rank by the engine's own score? |
| `best_dockq_in_top_n` | How accurate was the engine's usual top-N selection? |
| `acceptable_pose_generated` | Did the engine generate at least one DockQ ≥ 0.23 pose? |

This is not a metric-versus-metric scatter plot. It audits one engine's ranking
against an independent native-reference measure for the same poses. PPInsight
writes this table with `ppinsight quality <scores.tsv> <native.pdb> --rank-summary -o <quality.tsv>`.

For HADDOCK, prefer PDBs from the most refined available stage in this order:
`emref`, then `flexref`, then `rigidbody`. HADDOCK commonly stores successful
model PDBs as `.pdb.gz`; PPInsight's quality evaluator can read these compressed
files without permanently unpacking the run directory.

---

## 3. Cross-Model Normalization

### 3.1 The Problem

Engine-native scores are on incomparable scales:
- LightDock luciferin: positive values, higher = better, DFIRE units
- HADDOCK score: negative values, lower = better, weighted energy units
- Rosetta I_sc: negative values, lower = better, REU

Comparing raw values across engines is meaningless.

### 3.2 With a Native Reference: CAPRI Tier Proportions

Compare the percentage of each engine's poses in each CAPRI tier (high /
medium / acceptable / incorrect), as the `quality_bar` plot does.

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

## 4. Key Interpretation Notes

### 4.1 Engine Scores ≠ Binding Affinity

None of the three engine-native scores (luciferin, HADDOCK score, Rosetta I_sc) are calibrated predictions of binding affinity (ΔG). They are **internal ranking functions** optimized for each engine's search algorithm:
- LightDock luciferin drives the GSO swarm dynamics.
- HADDOCK score combines physics terms with restraint satisfaction.
- Rosetta I_sc is a difference of all-atom energies.

Comparing absolute values across engines is not meaningful.

### 4.2 Within-Run Ranking Is Valid; Cross-Run Comparison Is Not

All three engines are designed so that **better-scoring models within a single run for a single target** are more likely to be near-native. This is their intended use. Comparing an I_sc from one target to an I_sc from a different target, or comparing an I_sc to a HADDOCK score, is not supported by the scoring functions.

### 4.3 DockQ/CAPRI Is the Common Currency

The only scientifically defensible cross-engine comparison uses DockQ/CAPRI, which:
- Is computed independently from any engine's scoring function.
- Requires a known native structure (the "ground truth").
- Has published, well-accepted thresholds (Basu & Wallner 2016).

### 4.4 The "No Native" Case

Without a native structure:
- Engine-native scores are used for relative ranking only.
- Rank normalization can place scores on a common scale for visual comparison (for example, a CDF).
- **No quality labels should be assigned**, so `quality_bar` does not apply.

### 4.5 Score Directions

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
