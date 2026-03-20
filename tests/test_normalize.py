"""Tests for score normalization, aggregation, clustering, and new plot types."""

import pandas as pd
import numpy as np
import os
import pytest

from ppinsight.visualizer import (
    normalize_scores,
    get_metric_direction,
    METRIC_METADATA,
    violin_plot,
    score_heatmap,
    roc_curve_plot,
    rank_comparison_scatter,
)
from ppinsight.collect_scores import (
    _parse_lightdock_clusters,
    _parse_lightdock,
    aggregate_scores,
)

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for CI


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def multi_model_scores():
    """Scores from two models on two different scales."""
    return pd.DataFrame({
        "model": ["lightdock"] * 4 + ["haddock"] * 4,
        "score_type": ["luciferin_score"] * 4 + ["score"] * 4,
        "score_value": [10, 20, 30, 40, -100, -80, -60, -40],
        "proteinA": ["A"] * 8,
        "proteinB": ["B"] * 8,
    })


@pytest.fixture
def single_score_type():
    """Two models sharing the same score_type but different scales."""
    return pd.DataFrame({
        "model": ["m1", "m1", "m1", "m2", "m2", "m2"],
        "score_type": ["score"] * 6,
        "score_value": [0.0, 50.0, 100.0, -10.0, 0.0, 10.0],
        "proteinA": ["X"] * 6,
        "proteinB": ["Y"] * 6,
    })


@pytest.fixture
def labeled_scores():
    """Multi-model scores with interaction labels for plot tests."""
    rows = []
    for m in ("lightdock", "haddock"):
        st = "luciferin_score" if m == "lightdock" else "score"
        for pA, pB, label in [("A", "B", "interaction"), ("C", "D", "non-interaction")]:
            for val in range(5):
                sv = val * 10 if label == "interaction" else val * 2
                if m == "haddock":
                    sv = -sv  # HADDOCK: lower = better
                rows.append({
                    "model": m, "score_type": st, "score_value": float(sv),
                    "proteinA": pA, "proteinB": pB, "label": label,
                })
    return pd.DataFrame(rows)


@pytest.fixture
def lightdock_cluster_dir(tmp_path):
    """Create a fake LightDock output with cluster.repr files."""
    for i in range(3):
        swarm = tmp_path / f"swarm_{i}"
        swarm.mkdir()
        # cluster.repr format: cluster_id:population:best_scoring:glowworm_id:pdb_file
        (swarm / "cluster.repr").write_text(
            f"0:5:{20.0 - i * 3:.5f}:{i}:lightdock_{i}.pdb\n"
            f"1:2:{10.0 - i * 2:.5f}:{i + 10}:lightdock_{i + 10}.pdb\n"
        )
        # Also create a gso file for fallback testing
        (swarm / "gso_10.out").write_text(
            f"# header\n"
            f"(1.0, 2.0) 0 0 5.0 1 0.5 {15.0 + i:.3f}\n"
            f"(3.0, 4.0) 0 0 4.0 2 0.4 {12.0 + i:.3f}\n"
        )
    return tmp_path


@pytest.fixture
def lightdock_no_cluster_dir(tmp_path):
    """LightDock output WITHOUT cluster.repr (only gso_*.out)."""
    for i in range(2):
        swarm = tmp_path / f"swarm_{i}"
        swarm.mkdir()
        (swarm / "gso_10.out").write_text(
            f"# header\n"
            f"(1.0, 2.0) 0 0 5.0 1 0.5 {15.0 + i:.3f}\n"
        )
    return tmp_path


# =========================================================================
# Direction-aware normalization
# =========================================================================

class TestMetricMetadata:
    def test_known_metrics_have_direction(self):
        assert get_metric_direction("luciferin_score") is True
        assert get_metric_direction("score") is False  # HADDOCK
        assert get_metric_direction("dockq") is True
        assert get_metric_direction("irmsd") is False
        assert get_metric_direction("interface_score") is False  # Rosetta

    def test_unknown_metric_defaults_to_higher(self):
        assert get_metric_direction("unknown_metric_xyz") is True

    def test_new_rosetta_metrics_direction(self):
        """Audit #4/#10: New Rosetta metrics must be lower-is-better."""
        assert get_metric_direction("i_sc") is False
        assert get_metric_direction("total_score") is False
        assert get_metric_direction("irms") is False
        assert get_metric_direction("rms") is False

    def test_all_metadata_entries_have_required_keys(self):
        """Every entry in METRIC_METADATA must have higher_is_better and description."""
        for name, meta in METRIC_METADATA.items():
            assert "higher_is_better" in meta, f"{name} missing higher_is_better"
            assert "description" in meta, f"{name} missing description"
            assert isinstance(meta["higher_is_better"], bool)


class TestDirectionAwareMinMax:
    def test_flips_lower_is_better(self, multi_model_scores):
        """After direction-aware minmax, higher always means better binding."""
        result = normalize_scores(
            multi_model_scores, method="minmax", per_model=True,
            direction_aware=True,
        )
        # HADDOCK 'score' is lower-is-better → should be flipped
        # Original: -100, -80, -60, -40 (best = -100)
        # After flip: 100, 80, 60, 40 → minmax: 1.0, 0.67, 0.33, 0.0
        haddock = result[result["model"] == "haddock"]
        # The row that was originally -100 (best) should now be 1.0
        assert haddock["score_value"].max() == pytest.approx(1.0)
        assert haddock["score_value"].min() == pytest.approx(0.0)

    def test_no_flip_when_disabled(self, multi_model_scores):
        """direction_aware=False should NOT flip anything."""
        result = normalize_scores(
            multi_model_scores, method="minmax", per_model=True,
            direction_aware=False,
        )
        # Both groups should span [0, 1] without flipping
        for (model, st), grp in result.groupby(["model", "score_type"]):
            assert grp["score_value"].min() == pytest.approx(0.0)
            assert grp["score_value"].max() == pytest.approx(1.0)

    def test_lightdock_not_flipped(self, multi_model_scores):
        """LightDock luciferin_score is higher-is-better → no flip needed."""
        result = normalize_scores(
            multi_model_scores, method="minmax", per_model=True,
            direction_aware=True,
        )
        ld = result[result["model"] == "lightdock"]
        # Original order: 10, 20, 30, 40 → minmax: 0, 0.33, 0.67, 1.0
        # 40 was the best original → should still be 1.0
        assert ld["score_value"].max() == pytest.approx(1.0)


class TestMinMaxBasic:
    def test_per_model(self, multi_model_scores):
        result = normalize_scores(
            multi_model_scores, method="minmax", per_model=True,
            direction_aware=False,
        )
        for (model, st), grp in result.groupby(["model", "score_type"]):
            assert grp["score_value"].min() == pytest.approx(0.0)
            assert grp["score_value"].max() == pytest.approx(1.0)

    def test_global(self, single_score_type):
        result = normalize_scores(
            single_score_type, method="minmax", per_model=False,
            direction_aware=False,
        )
        vals = result["score_value"]
        assert vals.min() == pytest.approx(0.0)
        assert vals.max() == pytest.approx(1.0)

    def test_constant_values(self):
        df = pd.DataFrame({
            "model": ["m1"] * 3,
            "score_type": ["s"] * 3,
            "score_value": [5.0, 5.0, 5.0],
        })
        result = normalize_scores(df, method="minmax", direction_aware=False)
        assert (result["score_value"] == 0.5).all()


class TestZScore:
    def test_per_model(self, multi_model_scores):
        result = normalize_scores(
            multi_model_scores, method="zscore", per_model=True,
            direction_aware=False,
        )
        for (model, st), grp in result.groupby(["model", "score_type"]):
            assert grp["score_value"].mean() == pytest.approx(0.0, abs=1e-10)
            assert grp["score_value"].std(ddof=0) == pytest.approx(1.0, abs=1e-10)

    def test_constant_values(self):
        df = pd.DataFrame({
            "model": ["m1"] * 3,
            "score_type": ["s"] * 3,
            "score_value": [7.0, 7.0, 7.0],
        })
        result = normalize_scores(df, method="zscore", direction_aware=False)
        assert (result["score_value"] == 0.0).all()


class TestRank:
    def test_per_model(self, multi_model_scores):
        result = normalize_scores(
            multi_model_scores, method="rank", per_model=True,
            direction_aware=False,
        )
        for (model, st), grp in result.groupby(["model", "score_type"]):
            assert grp["score_value"].min() > 0
            assert grp["score_value"].max() <= 1.0

    def test_ordering_preserved(self, single_score_type):
        result = normalize_scores(
            single_score_type, method="rank", per_model=True,
            direction_aware=False,
        )
        for (model, st), grp in result.groupby(["model", "score_type"]):
            assert grp["score_value"].is_monotonic_increasing


class TestEdgeCases:
    def test_bad_method_raises(self, multi_model_scores):
        with pytest.raises(ValueError, match="method must be"):
            normalize_scores(multi_model_scores, method="bad")

    def test_non_numeric_coerced(self):
        df = pd.DataFrame({
            "model": ["m1", "m1"],
            "score_type": ["s", "s"],
            "score_value": ["hello", "10"],
        })
        result = normalize_scores(df, method="minmax", direction_aware=False)
        assert pd.isna(result["score_value"].iloc[0])

    def test_does_not_modify_input(self, multi_model_scores):
        original = multi_model_scores["score_value"].tolist()
        _ = normalize_scores(multi_model_scores, method="minmax")
        assert multi_model_scores["score_value"].tolist() == original

    def test_no_score_type_column(self):
        df = pd.DataFrame({
            "model": ["m1", "m1", "m1"],
            "score_value": [10.0, 20.0, 30.0],
        })
        result = normalize_scores(df, method="minmax", direction_aware=False)
        assert result["score_value"].min() == pytest.approx(0.0)
        assert result["score_value"].max() == pytest.approx(1.0)


# =========================================================================
# LightDock clustering parser
# =========================================================================

class TestLightdockClusters:
    def test_reads_cluster_repr(self, lightdock_cluster_dir):
        df = _parse_lightdock_clusters(str(lightdock_cluster_dir))
        # 3 swarms × 2 clusters each = 6 rows
        assert len(df) == 6
        assert "cluster_id" in df.columns
        assert "cluster_pop" in df.columns
        assert "swarm" in df.columns

    def test_best_cluster_score(self, lightdock_cluster_dir):
        df = _parse_lightdock_clusters(str(lightdock_cluster_dir))
        # swarm_0, cluster 0 should have score 20.0
        best = df[(df["swarm"] == "swarm_0") & (df["cluster_id"] == 0)]
        assert best["score_value"].iloc[0] == pytest.approx(20.0)

    def test_falls_back_without_clusters(self, lightdock_no_cluster_dir):
        """If no cluster.repr exists, falls back to raw gso parsing."""
        df = _parse_lightdock_clusters(str(lightdock_no_cluster_dir))
        assert len(df) > 0
        assert "cluster_id" not in df.columns  # raw parser doesn't add this


# =========================================================================
# Score aggregation
# =========================================================================

class TestAggregateScores:
    @pytest.fixture
    def many_poses(self):
        """10 poses for one pair from one model."""
        return pd.DataFrame({
            "model": ["lightdock"] * 10,
            "score_type": ["luciferin_score"] * 10,
            "score_value": list(range(1, 11)),  # 1..10
            "proteinA": ["A"] * 10,
            "proteinB": ["B"] * 10,
        })

    def test_best(self, many_poses):
        result = aggregate_scores(many_poses, strategy="best")
        assert len(result) == 1
        assert result["score_value"].iloc[0] == 10.0

    def test_topN_mean(self, many_poses):
        result = aggregate_scores(many_poses, strategy="topN_mean", n=3)
        assert len(result) == 1
        # Top 3: 10, 9, 8 → mean = 9.0
        assert result["score_value"].iloc[0] == pytest.approx(9.0)

    def test_median(self, many_poses):
        result = aggregate_scores(many_poses, strategy="median")
        assert result["score_value"].iloc[0] == pytest.approx(5.5)

    def test_mean(self, many_poses):
        result = aggregate_scores(many_poses, strategy="mean")
        assert result["score_value"].iloc[0] == pytest.approx(5.5)

    def test_bad_strategy_raises(self, many_poses):
        with pytest.raises(ValueError, match="strategy must be"):
            aggregate_scores(many_poses, strategy="unknown")

    def test_preserves_labels(self):
        df = pd.DataFrame({
            "model": ["m1"] * 3,
            "score_type": ["s"] * 3,
            "score_value": [1.0, 2.0, 3.0],
            "proteinA": ["A"] * 3,
            "proteinB": ["B"] * 3,
            "label": ["interaction"] * 3,
        })
        result = aggregate_scores(df, strategy="best")
        assert "label" in result.columns
        assert result["label"].iloc[0] == "interaction"

    def test_best_lower_is_better(self):
        """For lower-is-better metrics like HADDOCK score, 'best' = lowest."""
        df = pd.DataFrame({
            "model": ["haddock"] * 5,
            "score_type": ["score"] * 5,  # HADDOCK score: lower is better
            "score_value": [-120.0, -80.0, -150.0, -50.0, -100.0],
            "proteinA": ["A"] * 5,
            "proteinB": ["B"] * 5,
        })
        result = aggregate_scores(df, strategy="best")
        assert len(result) == 1
        # Best HADDOCK score is the most negative (lowest) value
        assert result["score_value"].iloc[0] == pytest.approx(-150.0)

    def test_topN_lower_is_better(self):
        """For lower-is-better metrics, topN_mean picks the N lowest scores."""
        df = pd.DataFrame({
            "model": ["rosetta"] * 5,
            "score_type": ["i_sc"] * 5,  # Rosetta I_sc: lower is better
            "score_value": [-30.0, -10.0, -25.0, -5.0, -20.0],
            "proteinA": ["A"] * 5,
            "proteinB": ["B"] * 5,
        })
        result = aggregate_scores(df, strategy="topN_mean", n=3)
        # Top 3 (best = most negative): -30, -25, -20 → mean = -25.0
        assert result["score_value"].iloc[0] == pytest.approx(-25.0)

    def test_best_higher_is_better(self):
        """For higher-is-better metrics like luciferin_score, 'best' = highest."""
        df = pd.DataFrame({
            "model": ["lightdock"] * 5,
            "score_type": ["luciferin_score"] * 5,
            "score_value": [5.0, 15.0, 10.0, 20.0, 8.0],
            "proteinA": ["A"] * 5,
            "proteinB": ["B"] * 5,
        })
        result = aggregate_scores(df, strategy="best")
        assert result["score_value"].iloc[0] == pytest.approx(20.0)

    def test_mixed_directions_per_group(self):
        """Aggregation should handle mixed directions in a single DataFrame."""
        df = pd.DataFrame({
            "model": ["lightdock"] * 3 + ["haddock"] * 3,
            "score_type": ["luciferin_score"] * 3 + ["score"] * 3,
            "score_value": [10.0, 20.0, 15.0,  -50.0, -100.0, -150.0],
            "proteinA": ["A"] * 6,
            "proteinB": ["B"] * 6,
        })
        result = aggregate_scores(df, strategy="best")
        assert len(result) == 2
        ld = result[result["score_type"] == "luciferin_score"]
        hd = result[result["score_type"] == "score"]
        # LightDock: higher is better → 20.0
        assert ld["score_value"].iloc[0] == pytest.approx(20.0)
        # HADDOCK: lower is better → -150.0
        assert hd["score_value"].iloc[0] == pytest.approx(-150.0)


# =========================================================================
# New plot types (smoke tests — verify they produce figures without error)
# =========================================================================

class TestViolinPlot:
    def test_basic(self, labeled_scores):
        fig = violin_plot(labeled_scores, metric="luciferin_score")
        assert fig is not None

    def test_split_by_label(self, labeled_scores):
        fig = violin_plot(labeled_scores, metric="luciferin_score", split_by_label=True)
        assert fig is not None

    def test_save_to_file(self, labeled_scores, tmp_path):
        out = str(tmp_path / "violin.png")
        violin_plot(labeled_scores, metric="luciferin_score", output=out)
        assert os.path.isfile(out)


class TestScoreHeatmap:
    def test_basic(self, labeled_scores):
        fig = score_heatmap(labeled_scores, metric="luciferin_score")
        assert fig is not None

    def test_save_to_file(self, labeled_scores, tmp_path):
        out = str(tmp_path / "heatmap.png")
        score_heatmap(labeled_scores, metric="luciferin_score", output=out)
        assert os.path.isfile(out)


class TestROCCurve:
    def test_basic(self, labeled_scores):
        fig = roc_curve_plot(labeled_scores, metric="luciferin_score")
        assert fig is not None

    def test_save_to_file(self, labeled_scores, tmp_path):
        out = str(tmp_path / "roc.png")
        roc_curve_plot(labeled_scores, metric="luciferin_score", output=out)
        assert os.path.isfile(out)


class TestRankScatter:
    def test_basic(self):
        df = pd.DataFrame({
            "model": ["m1"] * 3 + ["m2"] * 3,
            "score_type": ["s"] * 6,
            "score_value": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            "proteinA": ["A", "B", "C"] * 2,
            "proteinB": ["X", "Y", "Z"] * 2,
        })
        fig = rank_comparison_scatter(df, metric="s", model_x="m1", model_y="m2")
        assert fig is not None

    def test_no_common_pairs_raises(self):
        df = pd.DataFrame({
            "model": ["m1", "m2"],
            "score_type": ["s", "s"],
            "score_value": [1.0, 2.0],
            "proteinA": ["A", "C"],
            "proteinB": ["B", "D"],
        })
        with pytest.raises(ValueError, match="No common pairs"):
            rank_comparison_scatter(df, metric="s", model_x="m1", model_y="m2")
