"""Build a unified scores file from ALL available data in the repo.

Combines:
  - LightDock walkthrough (2UUY_rec vs 2UUY_lig, 182 swarms)
  - HADDOCK caprieval     (e2aP_1F3G vs hpr_ensemble, 3 stages)
  - Rosetta global dock   (COL_D vs IMM_D, expected output)

Produces:
  examples/all_scores.tsv   — unified scores for all 3 engines + 3 pairs
  examples/demo_pairs.tsv   — pairs file with interaction labels (for ROC)
"""

import os
import sys
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rows: list[dict] = []

# ══════════════════════════════════════════════════════════════════
# 1. LightDock — from walkthrough_scores.tsv (already collected)
# ══════════════════════════════════════════════════════════════════
wt = os.path.join(ROOT, "walkthrough_scores.tsv")
if os.path.exists(wt):
    df_wt = pd.read_csv(wt, sep="\t")
    # Keep only LightDock rows (HADDOCK in walkthrough was from a
    # different caprieval — we'll re-read HADDOCK below for the
    # correct pair name).
    ld = df_wt[df_wt["model"] == "lightdock"].copy()
    for _, r in ld.iterrows():
        rows.append({
            "model": "lightdock",
            "score_type": r["score_type"],
            "score_value": r["score_value"],
            "proteinA": "2UUY_rec",
            "proteinB": "2UUY_lig",
        })
    print(f"  LightDock: {len(ld)} rows (2UUY_rec vs 2UUY_lig)")

# ══════════════════════════════════════════════════════════════════
# 2. HADDOCK — caprieval stages (e2aP_1F3G vs hpr_ensemble)
# ══════════════════════════════════════════════════════════════════
haddock_metrics = ["score", "dockq", "irmsd", "fnat", "lrmsd"]
for stage in ["9_caprieval", "5_caprieval", "2_caprieval"]:
    capri = os.path.join(ROOT, "examples", "haddock3", "run1-test",
                         stage, "capri_ss.tsv")
    if not os.path.exists(capri):
        continue
    df_h = pd.read_csv(capri, sep="\t")
    for metric in haddock_metrics:
        if metric not in df_h.columns:
            continue
        vals = pd.to_numeric(df_h[metric], errors="coerce").dropna()
        for v in vals:
            rows.append({
                "model": "haddock",
                "score_type": metric,
                "score_value": float(v),
                "proteinA": "e2aP_1F3G",
                "proteinB": "hpr_ensemble",
            })
    print(f"  HADDOCK ({stage}): {len(df_h)} structures")
    break  # Use best stage (9_caprieval = final)

# Also add HADDOCK data for the 2UUY pair (from walkthrough)
if os.path.exists(wt):
    df_wt = pd.read_csv(wt, sep="\t")
    hd = df_wt[df_wt["model"] == "haddock"].copy()
    for _, r in hd.iterrows():
        rows.append({
            "model": "haddock",
            "score_type": r["score_type"],
            "score_value": r["score_value"],
            "proteinA": "2UUY_rec",
            "proteinB": "2UUY_lig",
        })
    print(f"  HADDOCK (walkthrough 2UUY): {len(hd)} rows")

# ══════════════════════════════════════════════════════════════════
# 3. Rosetta — global docking score file (COL_D vs IMM_D)
# ══════════════════════════════════════════════════════════════════
rosetta_sc = os.path.join(
    ROOT, "examples", "rosetta", "Protein-Protein-Docking",
    "output_files", "expected_output", "score_global_dock.sc",
)
if os.path.exists(rosetta_sc):
    # Rosetta .sc files are whitespace-delimited, first line is
    # SEQUENCE:, second is header, rest are data.
    with open(rosetta_sc) as f:
        lines = [l for l in f if l.startswith("SCORE:") and "total_score" not in l]
        header_line = [l for l in open(rosetta_sc) if l.startswith("SCORE:") and "total_score" in l]
    if header_line:
        cols = header_line[0].split()[1:]  # skip "SCORE:"
        for line in lines:
            vals = line.split()[1:]
            rec = dict(zip(cols, vals))
            for metric in ["total_score", "I_sc", "Irms", "Fnat", "rms"]:
                if metric in rec:
                    try:
                        v = float(rec[metric])
                    except ValueError:
                        continue
                    rows.append({
                        "model": "rosetta",
                        "score_type": metric.lower(),
                        "score_value": v,
                        "proteinA": "COL_D",
                        "proteinB": "IMM_D",
                    })
    print(f"  Rosetta: {len(lines)} decoys (COL_D vs IMM_D)")

# ══════════════════════════════════════════════════════════════════
# Assemble and write
# ══════════════════════════════════════════════════════════════════
df = pd.DataFrame(rows)
print(f"\nTotal rows: {len(df)}")
print(f"Models: {sorted(df['model'].unique())}")
print(f"Pairs: {df[['proteinA','proteinB']].drop_duplicates().values.tolist()}")
print(f"Metrics: {sorted(df['score_type'].unique())}")

out_scores = os.path.join(ROOT, "examples", "all_scores.tsv")
df.to_csv(out_scores, sep="\t", index=False)
print(f"\nSaved {out_scores}")

# ── Pairs file with interaction labels ────────────────────────────
# 2UUY is a known interacting complex, e2aP/hpr is known interacting,
# COL_D/IMM_D is a known interacting complex (colicin D / immunity protein).
# We'll also add a fake non-interacting pair by pairing cross-complex
# proteins — this gives ROC something to classify against.
pairs_rows = [
    {"proteinA": "2UUY_rec", "proteinB": "2UUY_lig", "label": "interaction"},
    {"proteinA": "e2aP_1F3G", "proteinB": "hpr_ensemble", "label": "interaction"},
    {"proteinA": "COL_D", "proteinB": "IMM_D", "label": "interaction"},
]
pairs_df = pd.DataFrame(pairs_rows)
out_pairs = os.path.join(ROOT, "examples", "demo_pairs.tsv")
pairs_df.to_csv(out_pairs, sep="\t", index=False)
print(f"Saved {out_pairs}")
