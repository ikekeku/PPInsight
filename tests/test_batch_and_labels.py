"""Tests for parse_pairs and batch_dock modules."""


import pandas as pd
import pytest

from ppinsight import parse_pairs

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_table(tmp_path):
    """Create a minimal RTK-interactome-style TSV."""
    tsv = tmp_path / "table.tsv"
    # Mimics the real format: leading tab, header rows, data rows
    tsv.write_text(
        "\tProteinA\t\tProtein B\t\t\t\t\t\t\t\t\n"
        "\tFamily\tMember\tInteractions\t\t\t\t\tNon-Interactions\t\t\t\n"
        "\tEphrin\tEPHA1\tERBB2 [6]\tFGFR2 [7]\t\t\t\tEGFR [7]\t\t\t\n"
        "\t\tEPHA2\tEGFR [5][8]\tERBB2 [5]\t\t\t\tMERTK [7]\t\t\t\n"
        "\tErbB\tEGFR\tFGFR1 [5]\t\t\t\t\t\t\t\t\n"
    )
    return str(tsv)


# ---------------------------------------------------------------------------
# parse_pairs
# ---------------------------------------------------------------------------

class TestParsePairs:
    def test_basic(self, sample_table):
        df = parse_pairs.parse_interaction_table(sample_table)
        assert "proteinA" in df.columns
        assert "proteinB" in df.columns
        assert "label" in df.columns
        assert "family" in df.columns
        assert len(df) > 0

    def test_interactions_and_non_interactions(self, sample_table):
        df = parse_pairs.parse_interaction_table(sample_table)
        labels = set(df["label"])
        assert "interaction" in labels
        assert "non-interaction" in labels

    def test_strips_references(self, sample_table):
        df = parse_pairs.parse_interaction_table(sample_table)
        # No brackets should remain in proteinB
        for val in df["proteinB"]:
            assert "[" not in val
            assert "]" not in val

    def test_references_extracted(self, sample_table):
        df = parse_pairs.parse_interaction_table(sample_table)
        # EGFR [5][8] should have references "5,8"
        egfr_row = df[(df["proteinA"] == "EPHA2") & (df["proteinB"] == "EGFR")]
        assert len(egfr_row) == 1
        assert "5" in egfr_row.iloc[0]["references"]
        assert "8" in egfr_row.iloc[0]["references"]

    def test_family_propagation(self, sample_table):
        df = parse_pairs.parse_interaction_table(sample_table)
        # EPHA2 is under Ephrin family (continuation row with no family)
        epha2_rows = df[df["proteinA"] == "EPHA2"]
        assert (epha2_rows["family"] == "Ephrin").all()

    def test_real_table(self):
        """Test against the real RTK interactome table if present."""
        import os
        real_path = os.path.join(
            os.path.dirname(__file__), "..",
            "RTK Interactome Tomasz Draft Table 1. Discovered cross-family "
            "interactions for PPInsight - Sheet1.tsv"
        )
        if not os.path.isfile(real_path):
            pytest.skip("Real RTK table not found")
        df = parse_pairs.parse_interaction_table(real_path)
        assert len(df) > 100  # should have 180+ pairs
        assert (df["label"] == "interaction").sum() > 50
        assert (df["label"] == "non-interaction").sum() > 20

    def test_empty_table(self, tmp_path):
        tsv = tmp_path / "empty.tsv"
        # Table with headers but no data rows at all
        tsv.write_text(
            "\tProteinA\t\tProtein B\t\t\t\t\t\t\t\t\n"
            "\tFamily\tMember\tInteractions\t\t\t\t\tNon-Interactions\t\t\t\n"
        )
        with pytest.raises(ValueError, match="No protein pairs"):
            parse_pairs.parse_interaction_table(str(tsv))


class TestParsePairsCLI:
    def test_basic(self, sample_table, tmp_path):
        out = str(tmp_path / "pairs.csv")
        parse_pairs.main([sample_table, "-o", out])
        df = pd.read_csv(out)
        assert len(df) > 0
        assert "label" in df.columns

    def test_tsv_output(self, sample_table, tmp_path):
        out = str(tmp_path / "pairs.tsv")
        parse_pairs.main([sample_table, "-o", out])
        df = pd.read_csv(out, sep="\t")
        assert len(df) > 0

    def test_stats_flag(self, sample_table, tmp_path, capsys):
        out = str(tmp_path / "pairs.csv")
        parse_pairs.main([sample_table, "-o", out, "--stats"])
        captured = capsys.readouterr()
        assert "Interactions:" in captured.out
        assert "Non-interactions:" in captured.out


# ---------------------------------------------------------------------------
# batch_dock (unit tests — no actual docking)
# ---------------------------------------------------------------------------

class TestBatchDock:
    def test_dry_run(self, tmp_path):
        from ppinsight import batch_dock

        # Create a minimal pairs DataFrame
        pairs_df = pd.DataFrame({
            "proteinA": ["EPHA1", "EPHA2"],
            "proteinB": ["ERBB2", "EGFR"],
            "label": ["interaction", "interaction"],
        })

        results = batch_dock.batch_dock(
            pairs_df,
            engines=["lightdock"],
            pdb_dir=str(tmp_path),  # no PDBs here
            dry_run=True,
        )
        assert len(results) == 2
        assert (results["status"] == "dry_run").all()

    def test_missing_pdbs(self, tmp_path):
        from ppinsight import batch_dock

        pairs_df = pd.DataFrame({
            "proteinA": ["FAKE_A"],
            "proteinB": ["FAKE_B"],
            "label": ["interaction"],
        })

        results = batch_dock.batch_dock(
            pairs_df,
            engines=["lightdock"],
            pdb_dir=str(tmp_path),
        )
        assert len(results) == 1
        assert results.iloc[0]["status"] == "pdb_missing_A"

    def test_limit(self, tmp_path):
        from ppinsight import batch_dock

        pairs_df = pd.DataFrame({
            "proteinA": ["A", "B", "C"],
            "proteinB": ["X", "Y", "Z"],
            "label": ["interaction"] * 3,
        })

        results = batch_dock.batch_dock(
            pairs_df,
            engines=["lightdock"],
            pdb_dir=str(tmp_path),
            limit=1,
            dry_run=True,
        )
        assert len(results) == 1


# ---------------------------------------------------------------------------
# Visualizer label-aware features
# ---------------------------------------------------------------------------

class TestCompareScoresByLabel:
    @pytest.fixture
    def labeled_scores(self):
        return pd.DataFrame({
            "model": ["lightdock"] * 4 + ["haddock"] * 4,
            "score_type": ["dockq"] * 8,
            "score_value": [0.8, 0.7, 0.3, 0.2, 0.9, 0.85, 0.4, 0.25],
            "proteinA": ["A", "A", "B", "B"] * 2,
            "proteinB": ["X", "X", "Y", "Y"] * 2,
            "label": ["interaction", "interaction",
                       "non-interaction", "non-interaction"] * 2,
        })

    def test_basic(self, labeled_scores):
        import matplotlib

        from ppinsight.visualizer import compare_scores_by_label
        matplotlib.use("Agg")
        fig = compare_scores_by_label(labeled_scores, "dockq", output=None)
        # Should produce a figure without error
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_save_to_file(self, labeled_scores, tmp_path):
        import matplotlib

        from ppinsight.visualizer import compare_scores_by_label
        matplotlib.use("Agg")
        out = str(tmp_path / "labeled.png")
        fig = compare_scores_by_label(labeled_scores, "dockq", output=out)
        assert os.path.isfile(out)
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_no_label_column(self):
        import matplotlib

        from ppinsight.visualizer import compare_scores_by_label
        matplotlib.use("Agg")
        df = pd.DataFrame({
            "model": ["lightdock"],
            "score_type": ["dockq"],
            "score_value": [0.8],
        })
        with pytest.raises(ValueError, match="label"):
            compare_scores_by_label(df, "dockq")


class TestClassificationSummary:
    @pytest.fixture
    def labeled_scores(self):
        return pd.DataFrame({
            "model": ["lightdock"] * 4 + ["haddock"] * 4,
            "score_type": ["dockq"] * 8,
            "score_value": [0.8, 0.7, 0.3, 0.2, 0.9, 0.85, 0.4, 0.25],
            "proteinA": ["A", "A", "B", "B"] * 2,
            "proteinB": ["X", "X", "Y", "Y"] * 2,
            "label": ["interaction", "interaction",
                       "non-interaction", "non-interaction"] * 2,
        })

    def test_basic(self, labeled_scores):
        from ppinsight.visualizer import classification_summary
        result = classification_summary(labeled_scores, "dockq", threshold=0.5)
        assert "model" in result.columns
        assert "TP" in result.columns
        assert "accuracy" in result.columns
        assert len(result) == 2  # two models

    def test_auto_threshold(self, labeled_scores):
        from ppinsight.visualizer import classification_summary
        result = classification_summary(labeled_scores, "dockq")
        assert result["threshold"].iloc[0] is not None

    def test_lower_is_better(self, labeled_scores):
        from ppinsight.visualizer import classification_summary
        result = classification_summary(
            labeled_scores, "dockq", threshold=0.5, higher_is_better=False
        )
        assert len(result) == 2

    def test_auto_direction_dockq(self, labeled_scores):
        """DockQ is higher-is-better; auto-detection should pick that up."""
        from ppinsight.visualizer import classification_summary
        # Don't pass higher_is_better — let it auto-detect.
        # With threshold=0.5 and dockq (higher=better), interaction pairs
        # (scores 0.8, 0.7 / 0.9, 0.85 > 0.5) should be TP,
        # non-interaction pairs (0.3, 0.2 / 0.4, 0.25 < 0.5) should be TN.
        result = classification_summary(labeled_scores, "dockq", threshold=0.5)
        for _, row in result.iterrows():
            assert row["TP"] == 1   # one interaction pair, correctly classified
            assert row["TN"] == 1   # one non-interaction pair, correctly classified
            assert row["FP"] == 0
            assert row["FN"] == 0

    def test_auto_direction_haddock_score_type(self):
        """HADDOCK 'score' is lower-is-better; auto-detect must flip direction."""
        from ppinsight.visualizer import classification_summary
        df = pd.DataFrame({
            "model": ["haddock"] * 4,
            "score_type": ["score"] * 4,
            # Lower score = better binding.  Interacting pair should have
            # the lower (more negative) score.
            "score_value": [-120.0, -110.0, -20.0, -10.0],
            "proteinA": ["A", "A", "B", "B"],
            "proteinB": ["X", "X", "Y", "Y"],
            "label": ["interaction", "interaction",
                       "non-interaction", "non-interaction"],
        })
        # threshold at -65 (median).  With lower-is-better auto-detected:
        # interaction pair mean = -115 ≤ -65 → predicted interacting (TP)
        # non-interaction pair mean = -15 > -65 → predicted non-interacting (TN)
        result = classification_summary(df, "score", threshold=-65.0)
        assert result["TP"].iloc[0] == 1
        assert result["TN"].iloc[0] == 1
        assert result["FP"].iloc[0] == 0
        assert result["FN"].iloc[0] == 0


import os  # noqa: E402
