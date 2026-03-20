"""
Tests for Rosetta decoy clustering (ppinsight.rosetta.cluster).

These tests use **FakePose** objects and a custom RMSD function injected
via ``pose_loader`` / ``rmsd_func`` parameters, so neither PyRosetta nor
real PDB files are required for the unit-test suite.
"""

import math
import textwrap
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Import guards — skip the entire module if scipy is missing
# ---------------------------------------------------------------------------
scipy = pytest.importorskip("scipy")

from ppinsight.rosetta.cluster import (  # noqa: E402
    cluster_decoys,
    best_of_largest_cluster,
    _pairwise_ca_rmsd,
)


# ---------------------------------------------------------------------------
# Helpers: fake poses & RMSD
# ---------------------------------------------------------------------------

class FakePose:
    """A trivial stand-in for a PyRosetta Pose, identified by *coords*."""

    def __init__(self, x: float, y: float = 0.0, z: float = 0.0, name: str = ""):
        self.x, self.y, self.z = x, y, z
        self.name = name

    def __repr__(self):
        return f"FakePose({self.name}, x={self.x})"


def _fake_ca_rmsd(pose_a: FakePose, pose_b: FakePose) -> float:
    """Euclidean distance between two FakePoses (pretend Cα-RMSD)."""
    return math.sqrt(
        (pose_a.x - pose_b.x) ** 2
        + (pose_a.y - pose_b.y) ** 2
        + (pose_a.z - pose_b.z) ** 2
    )


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def cluster_a_coords():
    """Three poses near the origin → should form cluster A."""
    return [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]


@pytest.fixture()
def cluster_b_coords():
    """Three poses shifted +50 Å → should form cluster B."""
    return [(50.0, 0.0, 0.0), (51.0, 0.0, 0.0), (50.0, 1.0, 0.0)]


@pytest.fixture()
def scores_df(cluster_a_coords, cluster_b_coords):
    """DataFrame with 6 decoys: cluster A (better scores) + cluster B."""
    rows = []
    for i, (_x, _y, _z) in enumerate(cluster_a_coords):
        rows.append({
            "description": f"decoy_a{i}",
            "total_score": -100 + i * 5,  # -100, -95, -90
            "i_sc": -30 + i,
        })
    for i, (_x, _y, _z) in enumerate(cluster_b_coords):
        rows.append({
            "description": f"decoy_b{i}",
            "total_score": -50 + i * 5,  # -50, -45, -40
            "i_sc": -10 + i,
        })
    return pd.DataFrame(rows)


@pytest.fixture()
def fake_poses(cluster_a_coords, cluster_b_coords):
    """Map description → FakePose."""
    poses = {}
    for i, (x, y, z) in enumerate(cluster_a_coords):
        poses[f"decoy_a{i}"] = FakePose(x, y, z, name=f"decoy_a{i}")
    for i, (x, y, z) in enumerate(cluster_b_coords):
        poses[f"decoy_b{i}"] = FakePose(x, y, z, name=f"decoy_b{i}")
    return poses


@pytest.fixture()
def pdb_dir(tmp_path, fake_poses):
    """Create dummy PDB files so cluster_decoys can find them by name."""
    for desc in fake_poses:
        (tmp_path / f"{desc}.pdb").write_text("ATOM  dummy\nEND\n")
    return tmp_path


def _make_pose_loader(fake_poses):
    """Return a pose_loader that looks up FakePoses by PDB stem."""
    def _loader(pdb_path):
        stem = Path(pdb_path).stem
        return fake_poses[stem]
    return _loader


# ---------------------------------------------------------------------------
# Tests: cluster_decoys
# ---------------------------------------------------------------------------

class TestClusterDecoys:
    """Core clustering functionality."""

    def test_basic_clustering(self, scores_df, pdb_dir, fake_poses):
        """Two well-separated groups should yield 2 clusters."""
        result = cluster_decoys(
            scores_df, pdb_dir,
            rmsd_cutoff=4.0, top_n=200,
            pose_loader=_make_pose_loader(fake_poses),
            rmsd_func=_fake_ca_rmsd,
        )
        assert "cluster" in result.columns
        assert "cluster_size" in result.columns
        assert "cluster_rank" in result.columns
        assert "is_top_of_cluster" in result.columns

        n_clusters = result["cluster"].nunique()
        assert n_clusters == 2, f"Expected 2 clusters, got {n_clusters}"

    def test_largest_cluster_is_zero(self, scores_df, pdb_dir, fake_poses):
        """Cluster 0 should be the largest (tie-break by best score)."""
        result = cluster_decoys(
            scores_df, pdb_dir,
            rmsd_cutoff=4.0,
            pose_loader=_make_pose_loader(fake_poses),
            rmsd_func=_fake_ca_rmsd,
        )
        # Both clusters have 3 members → tie broken by best score.
        # Cluster A has best score -100 < -50 → should be cluster 0.
        c0 = result[result["cluster"] == 0]
        assert all(d.startswith("decoy_a") for d in c0["description"]), (
            "Cluster 0 should contain the better-scoring cluster A decoys"
        )

    def test_top_of_cluster_flag(self, scores_df, pdb_dir, fake_poses):
        """Each cluster should have exactly one is_top_of_cluster=True row."""
        result = cluster_decoys(
            scores_df, pdb_dir,
            rmsd_cutoff=4.0,
            pose_loader=_make_pose_loader(fake_poses),
            rmsd_func=_fake_ca_rmsd,
        )
        tops = result[result["is_top_of_cluster"]]
        n_clusters = result["cluster"].nunique()
        assert len(tops) == n_clusters

    def test_too_few_decoys(self, pdb_dir, fake_poses):
        """A single decoy should get cluster 0 without errors."""
        single = pd.DataFrame([{
            "description": "decoy_a0",
            "total_score": -100,
            "i_sc": -30,
        }])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = cluster_decoys(
                single, pdb_dir,
                rmsd_cutoff=4.0,
                pose_loader=_make_pose_loader(fake_poses),
                rmsd_func=_fake_ca_rmsd,
            )
        assert len(result) == 1
        assert result.iloc[0]["cluster"] == 0

    def test_missing_score_column(self, scores_df, pdb_dir, fake_poses):
        """Requesting a non-existent score column should raise KeyError."""
        with pytest.raises(KeyError, match="nonexistent"):
            cluster_decoys(
                scores_df, pdb_dir,
                score_col="nonexistent",
                pose_loader=_make_pose_loader(fake_poses),
                rmsd_func=_fake_ca_rmsd,
            )

    def test_missing_description_column(self, pdb_dir, fake_poses):
        """Missing 'description' column should raise KeyError."""
        bad_df = pd.DataFrame({"total_score": [-100, -90]})
        with pytest.raises(KeyError, match="description"):
            cluster_decoys(
                bad_df, pdb_dir,
                pose_loader=_make_pose_loader(fake_poses),
                rmsd_func=_fake_ca_rmsd,
            )

    def test_top_n_limits_decoys(self, scores_df, pdb_dir, fake_poses):
        """top_n=3 should only cluster the 3 best-scoring decoys."""
        result = cluster_decoys(
            scores_df, pdb_dir,
            top_n=3, rmsd_cutoff=4.0,
            pose_loader=_make_pose_loader(fake_poses),
            rmsd_func=_fake_ca_rmsd,
        )
        assert len(result) == 3
        # The 3 best scores are from cluster A (-100, -95, -90)
        assert all(d.startswith("decoy_a") for d in result["description"])

    def test_cluster_rank_ordering(self, scores_df, pdb_dir, fake_poses):
        """Within each cluster, cluster_rank should be 0-indexed by score."""
        result = cluster_decoys(
            scores_df, pdb_dir,
            rmsd_cutoff=4.0,
            pose_loader=_make_pose_loader(fake_poses),
            rmsd_func=_fake_ca_rmsd,
        )
        for cid in result["cluster"].unique():
            subset = result[result["cluster"] == cid].sort_values("cluster_rank")
            ranks = list(subset["cluster_rank"])
            assert ranks == list(range(len(ranks)))


# ---------------------------------------------------------------------------
# Tests: best_of_largest_cluster
# ---------------------------------------------------------------------------

class TestBestOfLargestCluster:
    """Convenience function to extract the single best prediction."""

    def test_returns_best_in_cluster_zero(self, scores_df, pdb_dir, fake_poses):
        """Should return the best-scoring member of the largest cluster."""
        clustered = cluster_decoys(
            scores_df, pdb_dir,
            rmsd_cutoff=4.0,
            pose_loader=_make_pose_loader(fake_poses),
            rmsd_func=_fake_ca_rmsd,
        )
        best = best_of_largest_cluster(clustered)
        assert best["cluster"] == 0
        assert best["is_top_of_cluster"] == True  # noqa: E712

    def test_empty_cluster_raises(self):
        """Should raise if there is no cluster 0."""
        df = pd.DataFrame({
            "cluster": [1, 2],
            "is_top_of_cluster": [True, True],
            "total_score": [-10, -20],
        })
        with pytest.raises(ValueError, match="No cluster 0"):
            best_of_largest_cluster(df)


# ---------------------------------------------------------------------------
# Tests: _parse_rosetta_clustered (collect_scores integration)
# ---------------------------------------------------------------------------

class TestParseRosettaClustered:
    """Verify that collect_scores prefers clustered_scores.csv."""

    def test_parse_clustered_csv(self, tmp_path):
        """_parse_rosetta should pick up clustered_scores.csv."""
        from ppinsight.collect_scores import _parse_rosetta

        csv_path = tmp_path / "clustered_scores.csv"
        csv_path.write_text(textwrap.dedent("""\
            description,total_score,i_sc,cluster,cluster_size,cluster_rank,is_top_of_cluster
            decoy_1,-100,-30,0,5,0,True
            decoy_2,-95,-28,0,5,1,False
            decoy_3,-50,-10,1,3,0,True
        """))

        result = _parse_rosetta(str(tmp_path), pair=("recA", "ligB"))

        assert len(result) > 0
        # Only cluster reps (decoy_1 and decoy_3)
        assert "source" in result.columns
        assert all(result["source"] == "cluster")

    def test_clustered_takes_priority_over_sc(self, tmp_path):
        """When both clustered_scores.csv and .sc exist, clustered wins."""
        from ppinsight.collect_scores import _parse_rosetta

        csv_path = tmp_path / "clustered_scores.csv"
        csv_path.write_text(textwrap.dedent("""\
            description,total_score,i_sc,cluster,cluster_size,cluster_rank,is_top_of_cluster
            decoy_1,-100,-30,0,5,0,True
        """))

        sc_path = tmp_path / "scores.sc"
        sc_path.write_text(textwrap.dedent("""\
            SEQUENCE:
            SCORE: total_score I_sc description
            SCORE: -80.0 -20.0 decoy_sc1
        """))

        result = _parse_rosetta(str(tmp_path))
        assert "source" in result.columns
        assert all(result["source"] == "cluster")


# ---------------------------------------------------------------------------
# Tests: registry detection
# ---------------------------------------------------------------------------

class TestRegistryDetectsClusteredScores:
    """The Rosetta detector should recognize clustered_scores.csv."""

    def test_detect_with_clustered_csv(self, tmp_path):
        """Directory with clustered_scores.csv should be detected as Rosetta."""
        (tmp_path / "clustered_scores.csv").write_text("dummy\n")

        from ppinsight.collect_scores import detect_engine
        assert detect_engine(str(tmp_path)) == "rosetta"


# ---------------------------------------------------------------------------
# Tests: METRIC_METADATA entries
# ---------------------------------------------------------------------------

class TestMetricMetadata:
    """New clustering-related metrics should appear in METRIC_METADATA."""

    def test_dg_separated_in_metadata(self):
        from ppinsight.visualizer import METRIC_METADATA
        assert "dg_separated" in METRIC_METADATA
        assert METRIC_METADATA["dg_separated"]["higher_is_better"] is False

    def test_cluster_size_in_metadata(self):
        from ppinsight.visualizer import METRIC_METADATA
        assert "cluster_size" in METRIC_METADATA
        assert METRIC_METADATA["cluster_size"]["higher_is_better"] is True


# ---------------------------------------------------------------------------
# Tests: pairwise RMSD helper
# ---------------------------------------------------------------------------

class TestPairwiseCaRmsd:
    """Unit test for the pairwise RMSD matrix computation."""

    def test_symmetric_matrix(self):
        """Matrix should be symmetric with zeros on the diagonal."""
        poses = [FakePose(0), FakePose(5), FakePose(10)]
        mat = _pairwise_ca_rmsd(poses, rmsd_func=_fake_ca_rmsd)

        assert mat.shape == (3, 3)
        np.testing.assert_array_equal(np.diag(mat), [0, 0, 0])
        np.testing.assert_array_almost_equal(mat, mat.T)
        # Known distances: 0→5=5, 0→10=10, 5→10=5
        assert abs(mat[0, 1] - 5.0) < 1e-6
        assert abs(mat[0, 2] - 10.0) < 1e-6
        assert abs(mat[1, 2] - 5.0) < 1e-6

    def test_single_pose(self):
        """A single pose should give a 1×1 zero matrix."""
        poses = [FakePose(0)]
        mat = _pairwise_ca_rmsd(poses, rmsd_func=_fake_ca_rmsd)
        assert mat.shape == (1, 1)
        assert mat[0, 0] == 0.0
