"""
Rosetta decoy clustering — group docked structures by Cα-RMSD.

After a RosettaDock run produces thousands of decoys, many of them are
structurally near-identical.  Hierarchical clustering on pairwise Cα-RMSD
groups them into distinct binding modes.  The **best-scoring member** of
the **largest cluster** is the recommended prediction.

Algorithm
---------
1. Load scores, sort by total_score, take top *n* decoys.
2. Load each decoy PDB as a PyRosetta ``Pose``.
3. Compute an *n × n* pairwise Cα-RMSD matrix.
4. Convert to a condensed distance vector (``squareform``).
5. Hierarchical clustering (average linkage) with a distance cutoff
   (default 4 Å) using ``scipy.cluster.hierarchy``.
6. Renumber clusters so cluster 0 = largest (ties broken by best score).
7. Annotate the scores DataFrame with cluster membership.

Dependencies
------------
* **PyRosetta** — for ``CA_rmsd()``
* **scipy** — for hierarchical clustering
* **numpy** — array operations
* **pandas** — DataFrame manipulation
"""

from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Pairwise Cα-RMSD
# ---------------------------------------------------------------------------

def _pairwise_ca_rmsd(poses: Sequence, rmsd_func: Callable | None = None) -> np.ndarray:
    """Compute an *n × n* pairwise Cα-RMSD matrix over a list of Poses.

    Parameters
    ----------
    poses : sequence of pyrosetta.Pose (or any object understood by *rmsd_func*)
        The structures to compare.
    rmsd_func : callable, optional
        ``(pose_a, pose_b) -> float``.  Defaults to PyRosetta's ``CA_rmsd``.

    Returns
    -------
    np.ndarray
        Symmetric *n × n* matrix where entry *(i, j)* is the Cα-RMSD
        between ``poses[i]`` and ``poses[j]``.
    """
    if rmsd_func is None:
        from pyrosetta.rosetta.core.scoring import CA_rmsd  # type: ignore[import]
        rmsd_func = CA_rmsd

    n = len(poses)
    mat = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            rmsd = rmsd_func(poses[i], poses[j])
            mat[i, j] = rmsd
            mat[j, i] = rmsd
    return mat


# ---------------------------------------------------------------------------
# Default PDB loader
# ---------------------------------------------------------------------------

def _default_pose_loader(pdb_path: str | Path):
    """Load a PDB file as a PyRosetta Pose."""
    import pyrosetta  # type: ignore[import]
    return pyrosetta.pose_from_pdb(str(pdb_path))


# ---------------------------------------------------------------------------
# Main clustering entry point
# ---------------------------------------------------------------------------

def cluster_decoys(
    scores_df: pd.DataFrame,
    pdb_dir: str | Path,
    *,
    score_col: str = "total_score",
    top_n: int = 200,
    rmsd_cutoff: float = 4.0,
    pdb_pattern: str = "docked_{description}.pdb",
    pose_loader: Callable | None = None,
    rmsd_func: Callable | None = None,
) -> pd.DataFrame:
    """Cluster Rosetta decoys by Cα-RMSD and annotate the scores DataFrame.

    Parameters
    ----------
    scores_df : DataFrame
        Must have at least a *score_col* column and a ``description`` column
        (the standard Rosetta ``.sc`` identifier for each decoy).
    pdb_dir : path-like
        Directory containing the decoy PDB files.
    score_col : str
        Column to rank decoys by before taking the top *n*.  Lower is better.
    top_n : int
        How many top-scoring decoys to cluster (default 200).
    rmsd_cutoff : float
        RMSD threshold (Å) for ``fcluster(..., criterion='distance')``.
    pdb_pattern : str
        ``str.format``-style pattern with ``{description}`` placeholder,
        used to locate decoy PDB files.
    pose_loader : callable, optional
        ``(pdb_path) -> Pose``.  Defaults to ``pyrosetta.pose_from_pdb``.
        Useful for testing without PyRosetta.
    rmsd_func : callable, optional
        ``(pose_a, pose_b) -> float``.  Defaults to PyRosetta ``CA_rmsd``.
        Useful for testing without PyRosetta.

    Returns
    -------
    DataFrame
        Copy of the (top-*n* slice of) input with extra columns:

        * ``cluster`` — cluster id (0 = largest cluster)
        * ``cluster_size`` — number of members in that cluster
        * ``cluster_rank`` — rank within the cluster (0 = best)
        * ``is_top_of_cluster`` — True for the best-scoring decoy per cluster
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    if pose_loader is None:
        pose_loader = _default_pose_loader
    # rmsd_func=None is handled inside _pairwise_ca_rmsd

    pdb_dir = Path(pdb_dir)

    # ── 1. Sort & trim ───────────────────────────────────────────
    if score_col not in scores_df.columns:
        raise KeyError(
            f"Score column '{score_col}' not found in DataFrame. "
            f"Available: {list(scores_df.columns)}"
        )
    if "description" not in scores_df.columns:
        raise KeyError(
            "DataFrame must have a 'description' column identifying each decoy."
        )

    df = scores_df.sort_values(score_col, ascending=True).head(top_n).copy()
    df = df.reset_index(drop=True)

    if len(df) < 2:
        warnings.warn(
            f"Only {len(df)} decoy(s) available — skipping clustering.",
            stacklevel=2,
        )
        df["cluster"] = 0
        df["cluster_size"] = len(df)
        df["cluster_rank"] = range(len(df))
        df["is_top_of_cluster"] = True
        return df

    # ── 2. Load PDB poses ────────────────────────────────────────
    poses: list = []
    valid_idx: list[int] = []
    for idx, row in df.iterrows():
        desc = row["description"]
        # Try multiple naming conventions
        candidates = [
            pdb_dir / pdb_pattern.format(description=desc),
            pdb_dir / f"{desc}.pdb",
            pdb_dir / f"docked_{desc}.pdb",
            pdb_dir / f"decoy_{desc}.pdb",
        ]
        pdb_path = None
        for c in candidates:
            if c.is_file():
                pdb_path = c
                break
        if pdb_path is None:
            logger.warning("PDB file not found for decoy '%s' — skipping", desc)
            continue
        pose = pose_loader(pdb_path)
        poses.append(pose)
        valid_idx.append(idx)

    if len(poses) < 2:
        warnings.warn(
            f"Only {len(poses)} decoy PDB(s) found — skipping clustering.",
            stacklevel=2,
        )
        df["cluster"] = 0
        df["cluster_size"] = len(df)
        df["cluster_rank"] = range(len(df))
        df["is_top_of_cluster"] = True
        return df

    # Trim DataFrame to only the decoys we could load
    df = df.loc[valid_idx].reset_index(drop=True)

    # ── 3. Pairwise RMSD ────────────────────────────────────────
    rmsd_mat = _pairwise_ca_rmsd(poses, rmsd_func=rmsd_func)

    # ── 4. Hierarchical clustering ───────────────────────────────
    condensed = squareform(rmsd_mat, checks=False)
    Z = linkage(condensed, method="average")
    raw_labels = fcluster(Z, t=rmsd_cutoff, criterion="distance")
    # raw_labels are 1-based from scipy

    # ── 5. Renumber: 0 = largest cluster ─────────────────────────
    df["_raw_cluster"] = raw_labels
    cluster_sizes = df["_raw_cluster"].value_counts()  # descending
    # Break ties by best score in that cluster
    cluster_best = df.groupby("_raw_cluster")[score_col].min()
    sorted_clusters = (
        pd.DataFrame({"size": cluster_sizes, "best": cluster_best})
        .sort_values(["size", "best"], ascending=[False, True])
    )
    remap = {old: new for new, old in enumerate(sorted_clusters.index)}
    df["cluster"] = df["_raw_cluster"].map(remap)
    df.drop(columns="_raw_cluster", inplace=True)

    # ── 6. Per-cluster stats ─────────────────────────────────────
    df["cluster_size"] = df["cluster"].map(df["cluster"].value_counts())
    df["cluster_rank"] = df.groupby("cluster")[score_col].rank(
        method="first", ascending=True
    ).astype(int) - 1
    df["is_top_of_cluster"] = df["cluster_rank"] == 0

    return df


# ---------------------------------------------------------------------------
# Convenience wrappers
# ---------------------------------------------------------------------------

def best_of_largest_cluster(clustered_df: pd.DataFrame) -> pd.Series:
    """Return the single best-scoring member of cluster 0 (the largest).

    Parameters
    ----------
    clustered_df : DataFrame
        Output of :func:`cluster_decoys`.

    Returns
    -------
    Series
        The row corresponding to the best decoy in the largest cluster.
    """
    c0 = clustered_df[clustered_df["cluster"] == 0]
    if c0.empty:
        raise ValueError("No cluster 0 in the DataFrame")
    return c0.loc[c0["is_top_of_cluster"]].iloc[0]
