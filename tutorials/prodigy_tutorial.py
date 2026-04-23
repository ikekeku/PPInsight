#!/usr/bin/env python3
"""PRODIGY binding-affinity tutorial for PPInsight.

This script demonstrates how to:
  1. Fetch a PDB structure by PDB ID using Bio.PDB.PDBList directly.
  2. Score the docked complex with PRODIGY to obtain predicted ΔG and Kd.
  3. Append PRODIGY scores to an existing unified scores DataFrame.
  4. Visualise the distribution of predicted ΔG values with a CDF plot.

The example uses the ErbB2/ErbB3 receptor complex (PDB: 1IVO, UniProt
P04626 + P21860), a well-characterised oncogenic protein-protein
interaction.

Note: For accession-based fetching (UniProt accession → PDB → chain
selection), use ``ppinsight.protein_fetch.get_uniprot_data``.

Prerequisites
-------------
    pip install "ppinsight[prodigy]"
    # biopython already in the core dependencies

Usage
-----
    cd PPInsight
    python tutorials/prodigy_tutorial.py
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

from Bio.PDB import PDBList

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

# ---------------------------------------------------------------------------
# Step 0: Check optional dependency
# ---------------------------------------------------------------------------

try:
    import prodigy_prot  # noqa: F401
except ImportError:
    print(
        "ERROR: prodigy-prot is not installed.\n"
        "Install it with:\n\n"
        "    pip install ppinsight[prodigy]\n",
        file=sys.stderr,
    )
    sys.exit(1)

# ---------------------------------------------------------------------------
# Step 1: Fetch PDB structure by PDB ID
# ---------------------------------------------------------------------------

PDB_ID = "1IVO"  # ErbB2/ErbB3 heterodimer
OUT_DIR = ROOT / "tutorials" / "prodigy_output"
OUT_DIR.mkdir(parents=True, exist_ok=True)

pdb_path = OUT_DIR / f"{PDB_ID}.pdb"

print(f"Step 1 — Fetching {PDB_ID} from the PDB …")
try:
    pdbl = PDBList()
    raw_path = pdbl.retrieve_pdb_file(PDB_ID, pdir=str(OUT_DIR), file_format="pdb")
    if raw_path and Path(raw_path).exists():
        shutil.copy(raw_path, pdb_path)
        print(f"  Saved to {pdb_path.relative_to(ROOT)}")
    else:
        raise FileNotFoundError(f"retrieve_pdb_file returned no file for {PDB_ID}")
except Exception as exc:
    print(f"  WARNING: Could not fetch {PDB_ID}: {exc}")
    print("  Proceeding with mock data for illustration purposes.")
    pdb_path = None

# ---------------------------------------------------------------------------
# Step 2: Score the complex with PRODIGY
# ---------------------------------------------------------------------------

from ppinsight.prodigy import score_pdb  # noqa: E402

print("\nStep 2 — Scoring with PRODIGY …")
if pdb_path and pdb_path.exists():
    # 1IVO has chains A (ErbB2) and B (ErbB3)
    result = score_pdb(pdb_path, chains=["A", "B"], temperature=25.0)

    if "error" in result:
        print(f"  PRODIGY error: {result['error']}")
        print("  (This can happen if chains are not in contact in the PDB.)")
    else:
        print(f"  PDB          : {PDB_ID}")
        print(f"  ΔG (kcal/mol): {result['prodigy_ddg']:.2f}")
        print(f"  Kd (M)       : {result['prodigy_kd']:.2e}")
        print(f"  NIS-aliphatic: {result['nis_a']:.3f}")
        print(f"  NIS-charged  : {result['nis_c']:.3f}")
        print(f"  Contacts     : {result['n_contacts']}")
else:
    print("  Skipped (PDB not available).")
    result = {"prodigy_ddg": -9.4, "prodigy_kd": 2.1e-7,
               "nis_a": 0.28, "nis_c": 0.19, "n_contacts": 48}
    print("  Using mock values for illustration:")
    print(f"    ΔG = {result['prodigy_ddg']:.2f} kcal/mol")
    print(f"    Kd = {result['prodigy_kd']:.2e} M")

# ---------------------------------------------------------------------------
# Step 3: Append PRODIGY scores to a scores DataFrame
# ---------------------------------------------------------------------------

import pandas as pd  # noqa: E402,I001
import numpy as np  # noqa: E402,I001
from ppinsight.prodigy import add_prodigy_to_scores  # noqa: E402,I001

print("\nStep 3 — Demonstrating add_prodigy_to_scores() …")

# Build a small fabricated scores DataFrame to illustrate the API.
# In real usage you would load this from all_scores.tsv.
rng = np.random.default_rng(0)
models = ["HADDOCK", "LightDock", "Rosetta"]
fabricated = []
for model in models:
    for i in range(10):
        fabricated.append({
            "model": model,
            "score_type": "score",
            "score_value": float(rng.normal(-80, 20)),
            "proteinA": "ErbB2",
            "proteinB": "ErbB3",
            "pdb": f"{model.lower()}_{i:02d}.pdb",
        })
demo_scores = pd.DataFrame(fabricated)

# Stub: use mocked score_pdb so we do not need real PDB files in the tutorial.
from unittest.mock import patch  # noqa: E402

ddg_values = iter(rng.normal(-9.0, 1.5, len(demo_scores["pdb"].unique())))

def _mock_score(pdb_path, **kwargs):
    return {"prodigy_ddg": float(next(ddg_values)), "prodigy_kd": 1e-7}

with patch("ppinsight.prodigy.score_pdb", side_effect=_mock_score):
    scored = add_prodigy_to_scores(demo_scores, pdb_dir=OUT_DIR)

summary = (
    scored.groupby("model")["prodigy_ddg"]
    .agg(["mean", "std", "count"])
    .rename(columns={"mean": "mean_ΔG", "std": "std_ΔG", "count": "n"})
    .sort_values("mean_ΔG")
)
print("\n  Per-engine PRODIGY ΔG summary:")
print(summary.to_string())

# ---------------------------------------------------------------------------
# Step 4: CDF plot of predicted ΔG
# ---------------------------------------------------------------------------

import matplotlib  # noqa: E402,I001
matplotlib.use("Agg")

from ppinsight.visualizer import cdf_plot  # noqa: E402

print("\nStep 4 — Generating CDF plot of predicted ΔG …")

# Reshape to long-format for the visualizer
long_df = scored[["model", "proteinA", "proteinB", "prodigy_ddg"]].copy()
long_df = long_df.rename(columns={"prodigy_ddg": "score_value"})
long_df["score_type"] = "prodigy_ddg"

plot_path = OUT_DIR / "prodigy_ddg_cdf.png"
cdf_plot(long_df, metric="prodigy_ddg", output=str(plot_path))
print(f"  Saved to {plot_path.relative_to(ROOT)}")

print("\nDone!  All outputs written to tutorials/prodigy_output/")
