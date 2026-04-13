"""
Tests for the visualizer module.

Covers:
- to_plot (per-model file loading)
- load_scores (unified scores file)
- compare_scores (bar chart from frames)
- compare_scores_unified (bar chart from unified file)
- available_metrics / available_pairs
- CLI main()
- Edge cases and error handling
"""

import csv
import os

import pandas as pd
import pytest

from ppinsight import visualizer as vis

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def per_model_csvs(tmp_path):
    """Three tiny CSV files that mimic per-model docking score output."""
    paths = []
    for i in range(1, 4):
        p = tmp_path / f"model_{i}.csv"
        with open(p, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["run", "score"])
            for j in range(5):
                writer.writerow([j, -(10 + i) + j * 0.1])
        paths.append(str(p))
    return paths


@pytest.fixture
def per_model_tsvs(tmp_path):
    """Two TSV files mimicking HADDOCK capri_ss output (subset of columns)."""
    paths = []
    for i, name in enumerate(["haddock", "lightdock"], start=1):
        p = tmp_path / f"{name}_scores.tsv"
        with open(p, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["model", "score", "irmsd", "fnat", "dockq"])
            for j in range(4):
                writer.writerow([
                    f"pose_{j}",
                    -(100 + i * 10) + j,
                    1.0 + j * 0.5,
                    0.8 - j * 0.1,
                    0.7 - j * 0.05,
                ])
        paths.append(str(p))
    return paths


@pytest.fixture
def unified_scores_csv(tmp_path):
    """A unified scores.csv with rows for two models and two metrics."""
    p = tmp_path / "scores.csv"
    with open(p, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "proteinA", "proteinB", "model",
            "score_type", "score_value", "output_path",
        ])
        # HADDOCK dockq
        for v in [0.54, 0.82, 0.11]:
            writer.writerow(["2UUY_rec", "2UUY_lig", "HADDOCK", "dockq", v, "/path/h"])
        # Rosetta dockq
        for v in [0.45, 0.38]:
            writer.writerow(["2UUY_rec", "2UUY_lig", "Rosetta", "dockq", v, "/path/r"])
        # HADDOCK irmsd
        for v in [1.98, 1.03]:
            writer.writerow(["2UUY_rec", "2UUY_lig", "HADDOCK", "irmsd", v, "/path/h"])
        # Rosetta irmsd
        for v in [3.50]:
            writer.writerow(["2UUY_rec", "2UUY_lig", "Rosetta", "irmsd", v, "/path/r"])
        # Different protein pair
        writer.writerow(["1ABC_rec", "1ABC_lig", "HADDOCK", "dockq", 0.60, "/path/h2"])
    return str(p)


@pytest.fixture
def unified_scores_tsv(tmp_path):
    """Same data but as a .tsv."""
    p = tmp_path / "scores.tsv"
    with open(p, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["proteinA", "proteinB", "model", "score_type", "score_value"])
        writer.writerow(["2UUY_rec", "2UUY_lig", "HADDOCK", "dockq", 0.54])
        writer.writerow(["2UUY_rec", "2UUY_lig", "Rosetta", "dockq", 0.45])
    return str(p)


# ---------------------------------------------------------------------------
# to_plot tests
# ---------------------------------------------------------------------------

class TestToPlot:
    def test_not_a_list(self):
        """Passing a non-list raises TypeError."""
        with pytest.raises(TypeError, match="Format your input as a list"):
            vis.to_plot("single_string.csv")

    def test_loads_csv(self, per_model_csvs):
        frames = vis.to_plot(per_model_csvs)
        assert len(frames) == 3
        for df in frames:
            assert isinstance(df, pd.DataFrame)
            assert "score" in df.columns

    def test_loads_tsv(self, per_model_tsvs):
        frames = vis.to_plot(per_model_tsvs)
        assert len(frames) == 2
        assert "dockq" in frames[0].columns


# ---------------------------------------------------------------------------
# load_scores tests
# ---------------------------------------------------------------------------

class TestLoadScores:
    def test_csv(self, unified_scores_csv):
        df = vis.load_scores(unified_scores_csv)
        assert "model" in df.columns
        assert "score_value" in df.columns
        assert len(df) == 9  # 3+2+2+1+1

    def test_tsv(self, unified_scores_tsv):
        df = vis.load_scores(unified_scores_tsv)
        assert len(df) == 2

    def test_missing_columns(self, tmp_path):
        """File without required columns raises ValueError."""
        p = tmp_path / "bad.csv"
        p.write_text("col_a,col_b\n1,2\n")
        with pytest.raises(ValueError, match="missing required column"):
            vis.load_scores(str(p))


# ---------------------------------------------------------------------------
# compare_scores tests (per-model frames)
# ---------------------------------------------------------------------------

class TestCompareScores:
    def test_basic(self, per_model_csvs, monkeypatch):
        import matplotlib.pyplot as _plt
        monkeypatch.setattr(_plt, "show", lambda: None)

        frames = vis.to_plot(per_model_csvs)
        names = ["Model A", "Model B", "Model C"]
        fig = vis.compare_scores(frames, names, "score", "Test Title")
        assert fig is not None
        _plt.close(fig)

    def test_missing_score_type(self, per_model_csvs):
        frames = vis.to_plot(per_model_csvs)
        with pytest.raises(LookupError, match="does not exist"):
            vis.compare_scores(frames, ["a", "b", "c"], "nonexistent", "t")

    def test_mismatched_lengths(self, per_model_csvs):
        frames = vis.to_plot(per_model_csvs)
        with pytest.raises(ValueError, match="must match"):
            vis.compare_scores(frames, ["only_one"], "score", "t")

    def test_save_to_file(self, per_model_csvs, tmp_path, monkeypatch):
        import matplotlib.pyplot as _plt
        monkeypatch.setattr(_plt, "show", lambda: None)

        frames = vis.to_plot(per_model_csvs)
        out = str(tmp_path / "plot.png")
        fig = vis.compare_scores(
            frames, ["A", "B", "C"], "score", output=out,
        )
        assert os.path.isfile(out)
        _plt.close(fig)

    def test_tsv_metric(self, per_model_tsvs, monkeypatch):
        """Can plot a HADDOCK-style metric from TSV files."""
        import matplotlib.pyplot as _plt
        monkeypatch.setattr(_plt, "show", lambda: None)

        frames = vis.to_plot(per_model_tsvs)
        fig = vis.compare_scores(frames, ["HADDOCK", "LightDock"], "dockq")
        assert fig is not None
        _plt.close(fig)


# ---------------------------------------------------------------------------
# compare_scores_unified tests
# ---------------------------------------------------------------------------

class TestCompareScoresUnified:
    def test_basic(self, unified_scores_csv, monkeypatch):
        import matplotlib.pyplot as _plt
        monkeypatch.setattr(_plt, "show", lambda: None)

        df = vis.load_scores(unified_scores_csv)
        fig = vis.compare_scores_unified(df, "dockq")
        assert fig is not None
        _plt.close(fig)

    def test_filter_pair(self, unified_scores_csv, monkeypatch):
        import matplotlib.pyplot as _plt
        monkeypatch.setattr(_plt, "show", lambda: None)

        df = vis.load_scores(unified_scores_csv)
        fig = vis.compare_scores_unified(
            df, "dockq", pair=("2UUY_rec", "2UUY_lig"),
        )
        assert fig is not None
        _plt.close(fig)

    def test_no_matching_rows(self, unified_scores_csv):
        df = vis.load_scores(unified_scores_csv)
        with pytest.raises(ValueError, match="No rows match"):
            vis.compare_scores_unified(df, "nonexistent_metric")

    def test_save_to_file(self, unified_scores_csv, tmp_path, monkeypatch):
        import matplotlib.pyplot as _plt
        monkeypatch.setattr(_plt, "show", lambda: None)

        df = vis.load_scores(unified_scores_csv)
        out = str(tmp_path / "unified.png")
        fig = vis.compare_scores_unified(df, "dockq", output=out)
        assert os.path.isfile(out)
        _plt.close(fig)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

class TestHelpers:
    def test_available_metrics(self, unified_scores_csv):
        df = vis.load_scores(unified_scores_csv)
        metrics = vis.available_metrics(df)
        assert "dockq" in metrics
        assert "irmsd" in metrics

    def test_available_pairs(self, unified_scores_csv):
        df = vis.load_scores(unified_scores_csv)
        pairs = vis.available_pairs(df)
        assert len(pairs) == 2
        assert ("2uuy_rec", "2uuy_lig") in pairs or ("2UUY_rec", "2UUY_lig") in pairs


# ---------------------------------------------------------------------------
# CLI main()
# ---------------------------------------------------------------------------

class TestCLI:
    def test_help(self, capsys):
        with pytest.raises(SystemExit) as exc:
            vis.main(["--help"])
        assert exc.value.code == 0

    def test_unified_file(self, unified_scores_csv, tmp_path, monkeypatch):
        import matplotlib.pyplot as _plt
        monkeypatch.setattr(_plt, "show", lambda: None)

        out = str(tmp_path / "cli.png")
        vis.main([unified_scores_csv, "--metric", "dockq", "--output", out])
        assert os.path.isfile(out)

    def test_per_model_files(self, per_model_csvs, tmp_path, monkeypatch):
        import matplotlib.pyplot as _plt
        monkeypatch.setattr(_plt, "show", lambda: None)

        out = str(tmp_path / "multi.png")
        vis.main(per_model_csvs + ["--metric", "score", "--output", out])
        assert os.path.isfile(out)

    def test_list_metrics(self, unified_scores_csv, capsys):
        vis.main([unified_scores_csv, "--metric", "dockq", "--list-metrics"])
        captured = capsys.readouterr()
        assert "dockq" in captured.out

    def test_list_pairs(self, unified_scores_csv, capsys):
        vis.main([unified_scores_csv, "--metric", "dockq", "--list-pairs"])
        captured = capsys.readouterr()
        assert "2uuy_rec" in captured.out.lower()
