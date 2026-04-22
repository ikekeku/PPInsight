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
import sys
from pathlib import Path

import pandas as pd
import pytest

from ppinsight import visualizer as vis

_EXAMPLES_DIR = Path(__file__).resolve().parents[1] / "examples"

# ---------------------------------------------------------------------------
# Import build_fabricated_scores cleanly from the examples script.
# The script's side effects (OUT.mkdir, apply_theme, matplotlib backend) are
# all guarded inside main(), so a plain import is safe.
# ---------------------------------------------------------------------------
if str(_EXAMPLES_DIR.parent) not in sys.path:
    sys.path.insert(0, str(_EXAMPLES_DIR.parent))

import importlib.util as _ilu

_spec = _ilu.spec_from_file_location(
    "generate_example_plots",
    _EXAMPLES_DIR / "generate_example_plots.py",
)
_example_module = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_example_module)
build_fabricated_scores = _example_module.build_fabricated_scores
del _spec, _example_module, _ilu

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


# ---------------------------------------------------------------------------
# CAPRI pose-level classification
# ---------------------------------------------------------------------------

class TestCapriPoseClassification:
    def test_summary_uses_pose_level_rows(self):
        df = pd.DataFrame({
            "model": ["HADDOCK"] * 8,
            "score_type": ["dockq", "fnat", "irmsd", "lrmsd"] * 2,
            "score_value": [0.82, 0.72, 1.0, 2.4, 0.18, 0.06, 7.5, 12.0],
            "proteinA": ["A"] * 8,
            "proteinB": ["B"] * 8,
            "run_id": ["run_1"] * 8,
            "pose_id": ["pose_1"] * 4 + ["pose_2"] * 4,
        })

        summary = vis.capri_summary_table(df)

        assert int(summary.loc[0, "high"]) == 1
        assert int(summary.loc[0, "incorrect"]) == 1
        assert int(summary.loc[0, "total"]) == 2

    def test_classify_capri_propagates_pose_level_labels(self):
        df = pd.DataFrame({
            "model": ["HADDOCK"] * 8,
            "score_type": ["dockq", "fnat", "irmsd", "lrmsd"] * 2,
            "score_value": [0.82, 0.72, 1.0, 2.4, 0.18, 0.06, 7.5, 12.0],
            "proteinA": ["A"] * 8,
            "proteinB": ["B"] * 8,
            "run_id": ["run_1"] * 8,
            "pose_id": ["pose_1"] * 4 + ["pose_2"] * 4,
        })

        classified = vis.classify_capri(df)

        pose_1 = classified[classified["pose_id"] == "pose_1"]
        pose_2 = classified[classified["pose_id"] == "pose_2"]
        assert set(pose_1["capri_quality"]) == {"high"}
        assert set(pose_2["capri_quality"]) == {"incorrect"}


class TestMetricPresentation:
    def test_metric_display_name_uses_dockq_case(self):
        assert vis.get_metric_display_name("dockq") == "DockQ"

    def test_ridge_plot_respects_dockq_bounds(self, monkeypatch):
        import matplotlib.pyplot as _plt

        monkeypatch.setattr(_plt, "show", lambda: None)

        df = pd.DataFrame({
            "model": ["HADDOCK"] * 4 + ["Rosetta"] * 4,
            "score_type": ["dockq"] * 8,
            "score_value": [0.22, 0.45, 0.78, 0.96, 0.18, 0.33, 0.61, 0.88],
        })

        fig = vis.ridge_plot(df, metric="dockq")
        xmin, xmax = fig.axes[-1].get_xlim()
        assert xmin >= 0.0
        assert xmax <= 1.0
        assert fig.axes[-1].get_xlabel() == "DockQ"
        expected_thresholds = {
            round(threshold, 2) for threshold, _, _ in vis.DOCKQ_CAPRI_THRESHOLDS
        }
        assert fig.legends
        legend_labels = {text.get_text() for text in fig.legends[0].get_texts()}
        assert "Median" in legend_labels
        assert {
            f"CAPRI {label} ({threshold:.2f})"
            for threshold, label, _ in vis.DOCKQ_CAPRI_THRESHOLDS
        }.issubset(legend_labels)
        for ax in fig.axes:
            constant_x_lines = [
                line for line in ax.lines
                if len({round(float(x), 6) for x in line.get_xdata()}) == 1
            ]
            dashed_thresholds = {
                round(float(line.get_xdata()[0]), 2)
                for line in constant_x_lines
                if line.get_linestyle() == "--"
            }
            median_markers = [line for line in ax.lines if line.get_marker() == "o"]
            assert dashed_thresholds == expected_thresholds
            assert len(median_markers) == 1
        _plt.close(fig)

    def test_violin_plot_adds_faint_horizontal_gridlines(self, monkeypatch):
        import matplotlib.pyplot as _plt

        monkeypatch.setattr(_plt, "show", lambda: None)

        df = pd.DataFrame({
            "model": ["HADDOCK"] * 4 + ["Rosetta"] * 4,
            "score_type": ["dockq"] * 8,
            "score_value": [0.22, 0.45, 0.78, 0.96, 0.18, 0.33, 0.61, 0.88],
        })

        fig = vis.violin_plot(df, metric="dockq")
        fig.canvas.draw()
        gridlines = [
            line for line in fig.axes[0].yaxis.get_gridlines()
            if line.get_visible()
        ]
        assert gridlines
        assert fig.axes[0].get_axisbelow()
        _plt.close(fig)


class TestExampleGalleryData:
    def test_build_fabricated_scores_is_deterministic(self):
        df1 = build_fabricated_scores()
        df2 = build_fabricated_scores()
        pd.testing.assert_frame_equal(df1, df2)

    def test_build_fabricated_scores_keeps_dockq_in_bounds(self):
        df = build_fabricated_scores()
        dockq = df.loc[df["score_type"] == "dockq", "score_value"]
        assert ((dockq >= 0.0) & (dockq <= 1.0)).all()
