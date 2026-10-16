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

# ---------------------------------------------------------------------------
# Import build_fabricated_scores directly from the examples script.
# The script's side effects (OUT.mkdir, apply_theme, matplotlib backend) are
# all guarded inside main(), so a plain import is safe.
# ---------------------------------------------------------------------------
_EXAMPLES_DIR = Path(__file__).resolve().parents[1] / "examples"
sys.path.insert(0, str(_EXAMPLES_DIR))
try:
    from generate_example_plots import build_fabricated_scores  # noqa: E402
finally:
    sys.path.pop(0)

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

    def test_metric_overview_table(self, unified_scores_csv):
        df = vis.load_scores(unified_scores_csv)
        overview = vis.metric_overview_table(df)
        assert not overview.empty
        assert {"model", "pair", "rows", "metric_count", "metrics"} <= set(
            overview.columns
        )
        assert overview["metric_count"].max() >= 1


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

    def test_unified_file_pair_filter(self, unified_scores_csv, tmp_path, monkeypatch):
        import matplotlib.pyplot as _plt
        monkeypatch.setattr(_plt, "show", lambda: None)

        out = str(tmp_path / "pair.png")
        vis.main([
            unified_scores_csv,
            "--metric", "dockq",
            "--pair", "2UUY_rec:2UUY_lig",
            "--output", out,
        ])
        assert os.path.isfile(out)

    def test_unified_file_auto_saves_without_output(
        self,
        unified_scores_csv,
        tmp_path,
        monkeypatch,
    ):
        import matplotlib.pyplot as _plt

        monkeypatch.chdir(tmp_path)
        monkeypatch.setattr(
            _plt,
            "show",
            lambda: (_ for _ in ()).throw(AssertionError("show should not be called")),
        )

        vis.main([unified_scores_csv, "--metric", "dockq"])

        plots_dir = tmp_path / "data" / "output" / "plots"
        assert plots_dir.exists()
        assert list(plots_dir.glob("*.png"))

    def test_list_metrics(self, unified_scores_csv, capsys):
        vis.main([unified_scores_csv, "--metric", "dockq", "--list-metrics"])
        captured = capsys.readouterr()
        assert "dockq" in captured.out

    def test_list_pairs(self, unified_scores_csv, capsys):
        vis.main([unified_scores_csv, "--metric", "dockq", "--list-pairs"])
        captured = capsys.readouterr()
        assert "2uuy_rec" in captured.out.lower()

    def test_table_without_metric_prints_overview(self, unified_scores_csv, capsys):
        vis.main([unified_scores_csv, "--table"])
        captured = capsys.readouterr()
        assert "Metric-Agnostic Overview" in captured.out


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
        assert vis.get_metric_display_name("quality_dockq") == "DockQ"

    def test_haddock_dockq_is_named_as_haddocks_own(self):
        assert vis.get_metric_display_name("dockq") == "DockQ (HADDOCK caprieval)"
        note = vis.haddock_caprieval_note(["dockq", "score"])
        assert note is not None and "quality_dockq" in note
        assert vis.haddock_caprieval_note(["quality_dockq"]) is None

    @pytest.mark.parametrize("metric", ["dockq", "quality_dockq"])
    def test_strip_plot_respects_dockq_bounds(self, monkeypatch, metric):
        import matplotlib.pyplot as _plt

        monkeypatch.setattr(_plt, "show", lambda: None)

        df = pd.DataFrame({
            "model": ["haddock"] * 4 + ["rosetta"] * 3,
            "score_type": [metric] * 7,
            "score_value": [0.22, 0.45, 0.78, 0.96, 0.18, 0.33, 0.61],
        })

        fig = vis.strip_plot(df, metric=metric)
        ax = fig.axes[0]
        ymin, ymax = ax.get_ylim()
        assert (ymin, ymax) == (0.0, 1.0)
        assert ax.get_ylabel() == vis.get_metric_display_name(metric)
        tick_labels = [t.get_text() for t in ax.get_xticklabels()]
        assert tick_labels == ["HADDOCK3\n(n = 4)", "RosettaDock\n(n = 3)"]
        dashed_thresholds = {
            round(float(line.get_ydata()[0]), 2)
            for line in ax.lines
            if line.get_linestyle() == "--"
            and len({round(float(y), 6) for y in line.get_ydata()}) == 1
        }
        assert dashed_thresholds == {
            round(threshold, 2) for threshold, _, _ in vis.DOCKQ_CAPRI_THRESHOLDS
        }
        legend_labels = {text.get_text() for text in fig.legends[0].get_texts()}
        assert {
            f"CAPRI {label} ({threshold:.2f})"
            for threshold, label, _ in vis.DOCKQ_CAPRI_THRESHOLDS
        }.issubset(legend_labels)
        _plt.close(fig)

    def test_strip_plot_draws_every_pose(self, monkeypatch):
        import matplotlib.pyplot as _plt

        monkeypatch.setattr(_plt, "show", lambda: None)

        df = pd.DataFrame({
            "model": ["haddock"] * 5,
            "score_type": ["score"] * 5,
            "score_value": [-50.0, -40.0, -30.0, -20.0, -10.0],
        })

        fig = vis.strip_plot(df, metric="score")
        points = sum(
            len(collection.get_offsets()) for collection in fig.axes[0].collections
        )
        assert points == 5
        _plt.close(fig)

    def test_violin_plot_shows_counts_and_visible_points(self, monkeypatch):
        import matplotlib.pyplot as _plt

        monkeypatch.setattr(_plt, "show", lambda: None)

        df = pd.DataFrame({
            "model": ["haddock"] * 4 + ["rosetta"] * 4,
            "score_type": ["quality_dockq"] * 8,
            "score_value": [0.22, 0.45, 0.78, 0.96, 0.18, 0.33, 0.61, 0.88],
        })

        fig = vis.violin_plot(df, metric="quality_dockq")
        ax = fig.axes[0]
        tick_labels = [t.get_text() for t in ax.get_xticklabels()]
        assert tick_labels == ["HADDOCK3\n(n = 4)", "RosettaDock\n(n = 4)"]
        strip_points = [
            c for c in ax.collections
            if len(c.get_offsets()) == 4 and c.get_alpha() is not None
        ]
        assert strip_points
        assert all(c.get_alpha() >= 0.3 for c in strip_points)
        assert all(c.get_zorder() > 1 for c in strip_points)
        assert fig.legends
        from matplotlib.collections import PolyCollection
        bodies = [c for c in ax.collections if isinstance(c, PolyCollection)]
        assert bodies
        assert all(
            abs(face[3] - vis.VIOLIN_FILL_ALPHA) < 1e-6
            for body in bodies for face in body.get_facecolors()
        )
        _plt.close(fig)

    def test_capri_line_colors_match_quality_bar_tiers(self):
        for _, tier, color in vis.DOCKQ_CAPRI_THRESHOLDS:
            assert color == vis.CAPRI_TIER_COLORS[tier]

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


def _funnel_scores(pose_ids_for_dockq=None):
    """Unified scores for one target: engine scores + quality_dockq rows."""
    pose_ids = ["p1", "p2", "p3", "p4"]
    dockq_ids = pose_ids_for_dockq or pose_ids
    rows = []
    for pose_id, score in zip(pose_ids, [-90.0, -70.0, -50.0, -30.0], strict=True):
        rows.append({
            "model": "haddock", "score_type": "score", "score_value": score,
            "proteinA": "A", "proteinB": "B", "pose_id": pose_id,
        })
    for pose_id, luc in zip(pose_ids, [5.0, 25.0, 15.0, 10.0], strict=True):
        rows.append({
            "model": "lightdock", "score_type": "luciferin_score",
            "score_value": luc, "proteinA": "A", "proteinB": "B",
            "pose_id": pose_id,
        })
    for model in ("haddock", "lightdock"):
        for pose_id, dockq in zip(dockq_ids, [0.85, 0.40, 0.10, 0.05], strict=True):
            rows.append({
                "model": model, "score_type": "quality_dockq",
                "score_value": dockq, "proteinA": "A", "proteinB": "B",
                "pose_id": pose_id,
            })
    return pd.DataFrame(rows)


class TestFunnelPlot:
    def test_one_panel_per_engine_with_dockq_on_x(self, monkeypatch):
        import matplotlib.pyplot as _plt

        monkeypatch.setattr(_plt, "show", lambda: None)

        fig = vis.funnel_plot(_funnel_scores())
        panels = [ax for ax in fig.axes if ax.get_title()]
        assert len(panels) == 2
        assert panels[0].get_title().startswith("HADDOCK3")
        assert "best DockQ 0.85 at engine rank 1" in panels[0].get_title()
        assert panels[0].get_xlim() == (0.0, 1.0)
        # LightDock: higher luciferin is better, so the axis is inverted
        # to keep "better" pointing down in every panel.
        assert panels[1].yaxis_inverted()
        assert not panels[0].yaxis_inverted()
        assert panels[0].get_ylabel() == "Score (best at bottom)"
        _plt.close(fig)

    def test_requires_native_dockq_rows(self):
        df = _funnel_scores()
        df = df[df["score_type"] != "quality_dockq"]
        with pytest.raises(ValueError, match="quality_dockq"):
            vis.funnel_plot(df)

    def test_rejects_mismatched_pose_ids(self):
        df = _funnel_scores(pose_ids_for_dockq=["x1", "x2", "x3", "x4"])
        with pytest.raises(ValueError, match="pose_id"):
            vis.funnel_plot(df)

    def test_rejects_multiple_pairs(self):
        df = _funnel_scores()
        other = df.copy()
        other["proteinA"] = "C"
        with pytest.raises(ValueError, match="one target"):
            vis.funnel_plot(pd.concat([df, other], ignore_index=True))

    def test_cli_rejects_metric(self, tmp_path, capsys):
        path = tmp_path / "merged.tsv"
        _funnel_scores().to_csv(path, sep="\t", index=False)
        with pytest.raises(SystemExit) as exc_info:
            vis.main([str(path), "--plot-type", "funnel", "--metric", "score"])
        assert exc_info.value.code == 2
        assert "does not take --metric" in capsys.readouterr().err

    def test_cli_prints_rank_audit_and_saves(self, tmp_path, capsys):
        path = tmp_path / "merged.tsv"
        out = tmp_path / "funnel.png"
        _funnel_scores().to_csv(path, sep="\t", index=False)
        vis.main([str(path), "--plot-type", "funnel", "-o", str(out)])
        stdout = capsys.readouterr().out
        assert "Rank audit" in stdout
        assert "HADDOCK3" in stdout
        assert out.exists()


class TestTopNPlot:
    def test_table_tracks_best_dockq_by_engine_rank(self):
        table = vis.topn_table(_funnel_scores())
        haddock = table[table["model"] == "haddock"].set_index("N")["best_dockq"]
        lightdock = table[table["model"] == "lightdock"].set_index("N")["best_dockq"]
        # HADDOCK ranks its 0.85 pose first (lowest score).
        assert haddock[1] == pytest.approx(0.85)
        # LightDock ranks its 0.40 pose first; 0.85 is only reached at N = 4.
        assert lightdock[1] == pytest.approx(0.40)
        assert lightdock[4] == pytest.approx(0.85)

    def test_plot_has_one_curve_per_engine_on_dockq_scale(self, monkeypatch):
        import matplotlib.pyplot as _plt

        monkeypatch.setattr(_plt, "show", lambda: None)
        fig = vis.topn_plot(_funnel_scores())
        ax = fig.axes[0]
        assert ax.get_ylim() == (0.0, 1.0)
        assert ax.get_xscale() == "log"
        labels = {t.get_text() for t in fig.legends[0].get_texts()}
        assert any(label.startswith("HADDOCK3") for label in labels)
        assert any(label.startswith("LightDock") for label in labels)
        _plt.close(fig)

    def test_cli_rejects_metric(self, tmp_path, capsys):
        path = tmp_path / "merged.tsv"
        _funnel_scores().to_csv(path, sep="\t", index=False)
        with pytest.raises(SystemExit) as exc_info:
            vis.main([str(path), "--plot-type", "topn", "--metric", "score"])
        assert exc_info.value.code == 2
        assert "topn does not take --metric" in capsys.readouterr().err


class TestReferenceBasedGuards:
    def test_roc_refuses_dockq(self):
        df = pd.DataFrame({
            "model": ["haddock"] * 4,
            "score_type": ["quality_dockq"] * 4,
            "score_value": [0.8, 0.6, 0.2, 0.1],
            "proteinA": ["A", "B", "C", "D"], "proteinB": ["W", "X", "Y", "Z"],
            "label": ["interaction"] * 2 + ["non-interaction"] * 2,
        })
        with pytest.raises(ValueError, match="reference native structure"):
            vis.roc_curve_plot(df, metric="quality_dockq")

    def test_capri_tiers_prefer_native_dockq_over_haddocks_own(self):
        df = pd.DataFrame({
            "model": ["haddock", "haddock", "lightdock"],
            "score_type": ["dockq", "quality_dockq", "quality_dockq"],
            "score_value": [1.0, 0.1, 0.9],
            "pose_id": ["h1", "h1", "l1"],
        })
        table = vis._capri_pose_table(df).set_index("model")["capri_quality"]
        # HADDOCK's own 1.0 is ignored when a native-referenced value exists,
        # and LightDock is not hidden by HADDOCK's column.
        assert table["haddock"] == "incorrect"
        assert table["lightdock"] == "high"


class TestLargeRuns:
    def test_strip_draws_capped_subset_but_keeps_best_poses(self, monkeypatch):
        import matplotlib.pyplot as _plt
        import numpy as np

        monkeypatch.setattr(_plt, "show", lambda: None)
        n = vis.MAX_DRAWN_POSES + 1000
        values = np.linspace(0.0, 0.2, n)
        values[-1] = 0.95  # one rare near-native pose
        df = pd.DataFrame({
            "model": ["lightdock"] * n,
            "score_type": ["quality_dockq"] * n,
            "score_value": values,
        })
        fig = vis.strip_plot(df, metric="quality_dockq")
        ax = fig.axes[0]
        drawn = np.concatenate([c.get_offsets()[:, 1] for c in ax.collections])
        # A uniform sample of MAX_DRAWN_POSES, plus the rare acceptable pose
        # if the sample missed it.
        assert vis.MAX_DRAWN_POSES <= len(drawn) <= vis.MAX_DRAWN_POSES + 1
        assert drawn.max() == pytest.approx(0.95)
        assert "5,000 sampled" in ax.get_xticklabels()[0].get_text()
        legend_labels = {t.get_text() for t in fig.legends[0].get_texts()}
        assert any("always drawn" in label for label in legend_labels)
        _plt.close(fig)


class TestDisplayNames:
    def test_scatter_accepts_display_names(self, monkeypatch):
        import matplotlib.pyplot as _plt

        monkeypatch.setattr(_plt, "show", lambda: None)
        df = pd.DataFrame({
            "model": ["haddock", "rosetta"] * 3,
            "score_type": ["quality_dockq"] * 6,
            "score_value": [0.2, 0.3, 0.5, 0.4, 0.8, 0.9],
            "proteinA": ["A", "A", "C", "C", "E", "E"],
            "proteinB": ["B", "B", "D", "D", "F", "F"],
        })
        fig = vis.model_agreement_scatter(
            df, metric="quality_dockq", model_x="HADDOCK3", model_y="RosettaDock",
        )
        ax = fig.axes[0]
        assert ax.get_xlabel() == "HADDOCK3 — mean DockQ per pair"
        assert ax.get_xlim() == ax.get_ylim() == (0.0, 1.0)
        _plt.close(fig)

    @staticmethod
    def _scatter_df(n_pairs):
        rows = []
        for i in range(n_pairs):
            for model, offset in (("haddock", 0.10), ("rosetta", 0.12)):
                value = offset + 0.02 * i
                rows.append({
                    "model": model, "score_type": "quality_dockq",
                    "score_value": value, "proteinA": f"R{i}", "proteinB": f"L{i}",
                })
        return pd.DataFrame(rows)

    def test_scatter_warns_below_minimum_pairs(self, monkeypatch):
        import matplotlib.pyplot as _plt

        monkeypatch.setattr(_plt, "show", lambda: None)
        df = self._scatter_df(vis.MIN_SCATTER_PAIRS - 1)
        with pytest.warns(UserWarning, match="fewer than 10 pairs|Only 9"):
            fig = vis.model_agreement_scatter(
                df, metric="quality_dockq", model_x="haddock", model_y="rosetta",
            )
        assert any("unstable" in t.get_text() for t in fig.axes[0].texts)
        _plt.close(fig)

    def test_scatter_does_not_warn_at_minimum_pairs(self, monkeypatch):
        import warnings

        import matplotlib.pyplot as _plt

        monkeypatch.setattr(_plt, "show", lambda: None)
        df = self._scatter_df(vis.MIN_SCATTER_PAIRS)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            fig = vis.model_agreement_scatter(
                df, metric="quality_dockq", model_x="haddock", model_y="rosetta",
            )
        assert not any("unstable" in t.get_text() for t in fig.axes[0].texts)
        _plt.close(fig)

    def test_scatter_cli_prints_plain_warning(self, tmp_path, capsys):
        path = tmp_path / "scores.tsv"
        self._scatter_df(3).to_csv(path, sep="\t", index=False)
        vis.main([
            str(path), "--metric", "quality_dockq", "--plot-type", "scatter",
            "--models", "haddock", "rosetta", "-o", str(tmp_path / "s.png"),
        ])
        err = capsys.readouterr().err
        assert "WARNING: Only 3 protein pair(s)" in err
        assert (tmp_path / "s.png").exists()

    def test_quality_bar_puts_lowercase_n_under_engine(self, monkeypatch):
        import matplotlib.pyplot as _plt

        monkeypatch.setattr(_plt, "show", lambda: None)
        df = pd.DataFrame({
            "model": ["haddock"] * 3,
            "score_type": ["quality_dockq"] * 3,
            "score_value": [0.1, 0.5, 0.9],
            "pose_id": ["a", "b", "c"],
        })
        fig = vis.quality_bar_chart(df)
        ax = fig.axes[0]
        assert [t.get_text() for t in ax.get_xticklabels()] == ["HADDOCK3\n(n = 3)"]
        assert not any("N=" in t.get_text() for t in ax.texts)
        _plt.close(fig)

    def test_quality_bar_uses_dockq_tiers_and_names_the_pair(self, monkeypatch):
        import matplotlib.pyplot as _plt

        monkeypatch.setattr(_plt, "show", lambda: None)
        df = pd.DataFrame({
            "model": ["haddock"] * 3 + ["rosetta"] * 3,
            "score_type": ["quality_dockq"] * 6,
            "score_value": [0.1, 0.5, 0.9, 0.3, 0.3, 0.05],
            "proteinA": ["A"] * 6, "proteinB": ["B"] * 6,
            "pose_id": ["a", "b", "c", "d", "e", "f"],
        })
        fig = vis.quality_bar_chart(df)
        ax = fig.axes[0]
        assert ax.get_title() == (
            "Percentage of poses in each CAPRI tier, by engine (A vs B)"
        )
        labels = [t.get_text() for t in fig.legends[0].get_texts()]
        assert labels == [
            "Incorrect (DockQ < 0.23)",
            "Acceptable (0.23 ≤ DockQ < 0.49)",
            "Medium (0.49 ≤ DockQ < 0.80)",
            "High (DockQ ≥ 0.80)",
        ]
        assert fig.legends[0].get_title().get_text() == (
            "CAPRI tier (DockQ measured against the native structure)"
        )
        _plt.close(fig)

    def test_quality_bar_refuses_multiple_pairs(self):
        df = pd.DataFrame({
            "model": ["haddock"] * 2,
            "score_type": ["quality_dockq"] * 2,
            "score_value": [0.1, 0.9],
            "proteinA": ["A", "C"], "proteinB": ["B", "D"],
            "pose_id": ["a", "b"],
        })
        with pytest.raises(ValueError, match="one target"):
            vis.quality_bar_chart(df)

    def test_bar_chart_says_mean(self, monkeypatch):
        import matplotlib.pyplot as _plt

        monkeypatch.setattr(_plt, "show", lambda: None)
        frames = [pd.DataFrame({"quality_dockq": [0.2, 0.4]})]
        fig = vis.compare_scores(frames, ["haddock"], "quality_dockq")
        ax = fig.axes[0]
        assert ax.get_ylim() == (0.0, 1.0)
        assert ax.get_ylabel() == "Mean DockQ (± 1 SD)"
        assert ax.get_title() == "Mean DockQ by engine"
        _plt.close(fig)


class TestExampleGalleryData:
    def test_build_fabricated_scores_is_deterministic(self):
        df1 = build_fabricated_scores()
        df2 = build_fabricated_scores()
        pd.testing.assert_frame_equal(df1, df2)

    def test_build_fabricated_scores_keeps_dockq_in_bounds(self):
        df = build_fabricated_scores()
        dockq = df.loc[df["score_type"] == "quality_dockq", "score_value"]
        assert not dockq.empty
        assert ((dockq >= 0.0) & (dockq <= 1.0)).all()

    def test_build_fabricated_scores_has_no_dockq_for_non_binders(self):
        df = build_fabricated_scores()
        non_binders = df[df["label"] == "non-interaction"]
        assert not non_binders.empty
        assert not non_binders["score_type"].map(vis.is_reference_based).any()
