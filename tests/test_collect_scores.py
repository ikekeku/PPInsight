"""Tests for the collect_scores module."""

import os
import textwrap

import pandas as pd
import pytest

from ppinsight import collect_scores

# ---------------------------------------------------------------------------
# Fixtures — tiny fake output directories
# ---------------------------------------------------------------------------

@pytest.fixture
def haddock_dir(tmp_path):
    """Create a minimal HADDOCK caprieval output directory."""
    root = tmp_path / "haddock_run"
    capri_dir = root / "9_caprieval"
    capri_dir.mkdir(parents=True)
    tsv = capri_dir / "capri_ss.tsv"
    tsv.write_text(textwrap.dedent("""\
        model\tscore\tdockq\tirmsd\tfnat\tlrmsd
        ../model_1.pdb\t-112.5\t0.72\t1.8\t0.65\t3.2
        ../model_2.pdb\t-98.3\t0.68\t2.1\t0.59\t4.1
    """))
    return str(root)


@pytest.fixture
def lightdock_dir(tmp_path):
    """Create a minimal LightDock simulation directory."""
    root = tmp_path / "lightdock_sim"
    swarm = root / "swarm_0"
    swarm.mkdir(parents=True)
    gso = swarm / "gso_100.out"
    gso.write_text(textwrap.dedent("""\
        #Coordinates  RecID  LigID  Luciferin  Neighbor  Vision  Scoring
        (0.0, 0.0, 0.0)  0  0  27.83  6  0.32  18.66
        (1.0, 1.0, 1.0)  0  0  25.10  4  0.16  15.42
    """))
    return str(root)


@pytest.fixture
def rosetta_dir(tmp_path):
    """Create a minimal Rosetta output directory."""
    root = tmp_path / "rosetta_run"
    root.mkdir()
    csv_path = root / "docking_scores.csv"
    csv_path.write_text(
        "run,total_score,i_sc\n"
        "1,-109.8,-9.8\n"
        "2,-111.2,-11.2\n"
        "3,-107.5,-7.5\n"
    )
    return str(root)


@pytest.fixture
def rosetta_sc_dir(tmp_path):
    """Create a Rosetta output directory with a native .sc score file."""
    root = tmp_path / "rosetta_sc_run"
    root.mkdir()
    sc_path = root / "score.sc"
    sc_path.write_text(textwrap.dedent("""\
        SEQUENCE:
        SCORE: total_score I_sc Irms rms Fnat description
        SCORE:  -185.32  -12.45  1.80  3.50  0.72  decoy_1
        SCORE:  -172.10  -8.90  2.40  5.10  0.55  decoy_2
        SCORE:  -190.55  -15.22  1.20  2.80  0.85  decoy_3
    """))
    return str(root)


@pytest.fixture
def haddock_cluster_dir(tmp_path):
    """Create a HADDOCK run directory with clustfcc + post-cluster caprieval."""
    root = tmp_path / "haddock_cluster_run"
    # Pre-clustering caprieval
    capri_pre = root / "2_caprieval"
    capri_pre.mkdir(parents=True)
    (capri_pre / "capri_ss.tsv").write_text(textwrap.dedent("""\
        model\tscore\tdockq\tirmsd\tfnat\tlrmsd
        ../model_1.pdb\t-100.0\t0.60\t2.5\t0.50\t4.0
    """))
    # Clustering step
    clustfcc = root / "7_clustfcc"
    clustfcc.mkdir(parents=True)
    # Post-clustering caprieval (cluster-level results)
    capri_post = root / "9_caprieval"
    capri_post.mkdir(parents=True)
    (capri_post / "capri_ss.tsv").write_text(textwrap.dedent("""\
        model\tscore\tdockq\tirmsd\tfnat\tlrmsd
        ../cluster_1_best.pdb\t-130.5\t0.80\t1.2\t0.78\t2.1
        ../cluster_2_best.pdb\t-95.0\t0.55\t3.0\t0.42\t5.5
    """))
    return str(root)


# ---------------------------------------------------------------------------
# detect_engine
# ---------------------------------------------------------------------------

class TestDetectEngine:
    def test_haddock(self, haddock_dir):
        assert collect_scores.detect_engine(haddock_dir) == "haddock"

    def test_haddock_cluster(self, haddock_cluster_dir):
        assert collect_scores.detect_engine(haddock_cluster_dir) == "haddock"

    def test_lightdock(self, lightdock_dir):
        assert collect_scores.detect_engine(lightdock_dir) == "lightdock"

    def test_rosetta(self, rosetta_dir):
        assert collect_scores.detect_engine(rosetta_dir) == "rosetta"

    def test_rosetta_sc(self, rosetta_sc_dir):
        """Detect Rosetta by .sc score file."""
        assert collect_scores.detect_engine(rosetta_sc_dir) == "rosetta"

    def test_unknown(self, tmp_path):
        with pytest.raises(ValueError, match="Cannot detect"):
            collect_scores.detect_engine(str(tmp_path))


# ---------------------------------------------------------------------------
# Individual parsers
# ---------------------------------------------------------------------------

class TestParseHaddock:
    def test_basic(self, haddock_dir):
        df = collect_scores._parse_haddock(haddock_dir)
        assert set(df["score_type"].unique()) == {
            "score", "dockq", "irmsd", "fnat", "lrmsd",
        }
        # 2 models × 5 metrics = 10 rows
        assert len(df) == 10
        assert (df["model"] == "haddock").all()
        assert {"pose_id", "output_path", "source_file", "pose_rank"} <= set(df.columns)
        assert df["pose_id"].nunique() == 2

    def test_pair(self, haddock_dir):
        df = collect_scores._parse_haddock(haddock_dir, pair=("recA", "ligB"))
        assert (df["proteinA"] == "recA").all()
        assert (df["proteinB"] == "ligB").all()


class TestParseLightdock:
    def test_basic(self, lightdock_dir):
        df = collect_scores._parse_lightdock(lightdock_dir)
        assert len(df) == 2
        assert (df["score_type"] == "luciferin_score").all()
        assert (df["model"] == "lightdock").all()
        assert {"pose_id", "source_file", "pose_rank", "swarm"} <= set(df.columns)
        assert df["pose_id"].tolist() == ["swarm_0:lightdock_0", "swarm_0:lightdock_1"]

    def test_custom_label(self, lightdock_dir):
        df = collect_scores._parse_lightdock(lightdock_dir, label="ld_run1")
        assert (df["model"] == "ld_run1").all()


class TestParseRosetta:
    def test_basic(self, rosetta_dir):
        df = collect_scores._parse_rosetta(rosetta_dir)
        assert len(df) == 6
        assert set(df["score_type"].unique()) == {"i_sc", "total_score"}
        assert df["pose_id"].nunique() == 3

    def test_legacy_csv_maps_to_total_score(self, tmp_path):
        root = tmp_path / "rosetta_legacy_csv"
        root.mkdir()
        (root / "docking_scores.csv").write_text("run,score\n1,-9.8\n2,-11.2\n")

        df = collect_scores._parse_rosetta(str(root))

        assert len(df) == 2
        assert set(df["score_type"].unique()) == {"total_score"}

    def test_missing_csv(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="docking_scores.csv"):
            collect_scores._parse_rosetta(str(tmp_path))

    def test_sc_file(self, rosetta_sc_dir):
        """Native .sc files produce multiple metrics."""
        df = collect_scores._parse_rosetta(rosetta_sc_dir)
        # 3 decoys × 5 metrics = 15 rows
        assert len(df) == 15
        metrics = set(df["score_type"].unique())
        assert "i_sc" in metrics
        assert "total_score" in metrics
        assert "irms" in metrics
        assert "rms" in metrics
        assert "fnat" in metrics
        assert {"pose_id", "source_file", "pose_rank"} <= set(df.columns)
        assert df["pose_id"].nunique() == 3

    def test_sc_file_values(self, rosetta_sc_dir):
        """Check that I_sc values are correctly parsed from .sc file."""
        df = collect_scores._parse_rosetta(rosetta_sc_dir)
        i_sc_vals = sorted(df[df["score_type"] == "i_sc"]["score_value"].tolist())
        assert i_sc_vals == pytest.approx([-15.22, -12.45, -8.90])

    def test_sc_preferred_over_csv(self, tmp_path):
        """When both .sc and .csv exist, .sc takes precedence."""
        root = tmp_path / "mixed_rosetta"
        root.mkdir()
        (root / "docking_scores.csv").write_text("run,score\n1,-5.0\n")
        (root / "score.sc").write_text(textwrap.dedent("""\
            SEQUENCE:
            SCORE: total_score I_sc description
            SCORE:  -100.0  -10.0  decoy_1
        """))
        df = collect_scores._parse_rosetta(str(root))
        assert "i_sc" in df["score_type"].values


# ---------------------------------------------------------------------------
# HADDOCK cluster parsing (#7)
# ---------------------------------------------------------------------------

class TestParseHaddockClusters:
    def test_cluster_auto_detection(self, haddock_cluster_dir):
        """Default (no_clusters=False) picks up cluster-level results."""
        df = collect_scores._parse_haddock(haddock_cluster_dir)
        # Should get cluster-level results (2 models × 5 metrics = 10)
        assert len(df) == 10
        # Check it tagged as cluster source
        assert (df["source"] == "cluster").all()
        # Verify values from post-clustering caprieval
        dockq_vals = df[df["score_type"] == "dockq"]["score_value"].tolist()
        assert 0.80 in dockq_vals
        assert 0.55 in dockq_vals

    def test_no_clusters_flag(self, haddock_cluster_dir):
        """With no_clusters=True, falls back to per-model results."""
        df = collect_scores._parse_haddock(
            haddock_cluster_dir, no_clusters=True
        )
        # Should get per-model results from latest caprieval (still 9_caprieval)
        assert len(df) > 0
        # No 'source' column when using per-model path
        assert (
            "source" not in df.columns
            or not (df.get("source", "") == "cluster").all()
        )

    def test_fallback_without_clustfcc(self, haddock_dir):
        """Without clustfcc dir, falls back to per-model parsing."""
        df = collect_scores._parse_haddock(haddock_dir)
        assert len(df) == 10  # same as original test
        assert "source" not in df.columns


# ---------------------------------------------------------------------------
# collect()
# ---------------------------------------------------------------------------

class TestCollect:
    def test_mixed(self, haddock_dir, lightdock_dir):
        df, prov = collect_scores.collect([haddock_dir, lightdock_dir])
        assert set(df["model"].unique()) == {"haddock", "lightdock"}
        assert len(prov) >= 1  # at least one provenance entry
        assert "run_id" in df.columns
        assert df["run_id"].notna().all()

    def test_labels_override(self, haddock_dir, lightdock_dir):
        df, _prov = collect_scores.collect(
            [haddock_dir, lightdock_dir], labels=["h1", "ld1"]
        )
        assert set(df["model"].unique()) == {"h1", "ld1"}

    def test_label_count_mismatch(self, haddock_dir):
        with pytest.raises(ValueError, match="labels"):
            collect_scores.collect([haddock_dir], labels=["a", "b"])


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

class TestCLI:
    def test_default_output_auto_named(self, haddock_dir, tmp_path, monkeypatch):
        """Without -o, collect should auto-select a descriptive scores path."""
        monkeypatch.chdir(tmp_path)

        collect_scores.main([haddock_dir])

        out_dir = tmp_path / "data" / "output" / "scores"
        written = sorted(out_dir.glob("*.tsv"))
        assert written
        assert any(p.name.startswith("scores_") for p in written)

    def test_default_output_auto_avoids_overwrite(
        self,
        haddock_dir,
        tmp_path,
        monkeypatch,
    ):
        """Repeated no--output runs should create suffixed files, not overwrite."""
        monkeypatch.chdir(tmp_path)

        collect_scores.main([haddock_dir])
        collect_scores.main([haddock_dir])

        out_dir = tmp_path / "data" / "output" / "scores"
        written = sorted(out_dir.glob("*.tsv"))
        assert len(written) >= 2

    def test_basic_tsv(self, haddock_dir, tmp_path):
        out = str(tmp_path / "out.tsv")
        collect_scores.main([haddock_dir, "-o", out])
        assert os.path.isfile(out)
        df = pd.read_csv(out, sep="\t")
        assert "model" in df.columns
        assert "score_value" in df.columns
        assert "run_id" in df.columns

    def test_provenance_sidecar(self, haddock_dir, tmp_path):
        """CLI should write a .provenance.json sidecar next to the scores."""
        import json
        out = str(tmp_path / "out.tsv")
        collect_scores.main([haddock_dir, "-o", out])
        sidecar = out + ".provenance.json"
        assert os.path.isfile(sidecar)
        with open(sidecar) as fh:
            prov = json.load(fh)
        assert len(prov) >= 1
        # Every entry should have the standard keys
        for _rid, meta in prov.items():
            assert "engine" in meta
            assert "source_dir" in meta
            assert "collected_at" in meta

    def test_csv_output(self, rosetta_dir, tmp_path):
        out = str(tmp_path / "out.csv")
        collect_scores.main([rosetta_dir, "-o", out])
        df = pd.read_csv(out)
        assert len(df) == 6
        assert set(df["score_type"].unique()) == {"i_sc", "total_score"}

    def test_pair_flag(self, lightdock_dir, tmp_path):
        out = str(tmp_path / "out.tsv")
        collect_scores.main([lightdock_dir, "-o", out, "--pair", "recA:ligB"])
        df = pd.read_csv(out, sep="\t")
        assert (df["proteinA"] == "recA").all()
        assert (df["proteinB"] == "ligB").all()

    def test_summary_flag(self, haddock_dir, tmp_path, capsys):
        """Summary is now always printed (no flag needed)."""
        out = str(tmp_path / "out.tsv")
        collect_scores.main([haddock_dir, "-o", out])
        captured = capsys.readouterr()
        assert "Summary" in captured.out

    def test_bad_pair_format(self, haddock_dir, tmp_path):
        with pytest.raises(SystemExit):
            collect_scores.main([haddock_dir, "--pair", "nocolon"])

    def test_pairs_requires_nonempty_pair_context(
        self,
        lightdock_dir,
        tmp_path,
        capsys,
    ):
        """--pairs should fail fast when collected rows have blank protein names."""
        pairs = tmp_path / "pairs.tsv"
        pairs.write_text("proteinA\tproteinB\tlabel\nrecA\tligB\tinteraction\n")
        out = str(tmp_path / "out.tsv")

        with pytest.raises(SystemExit) as exc:
            collect_scores.main([lightdock_dir, "-o", out, "--pairs", str(pairs)])
        assert exc.value.code == 2

        captured = capsys.readouterr()
        assert "--pairs requires populated proteinA/proteinB" in captured.err
        assert "pass --pair" in captured.err

    def test_pairs_annotation_with_pair_flag(self, lightdock_dir, tmp_path):
        """When --pair is provided, --pairs annotation should label rows."""
        pairs = tmp_path / "pairs.tsv"
        pairs.write_text("proteinA\tproteinB\tlabel\nrecA\tligB\tinteraction\n")
        out = str(tmp_path / "out.tsv")

        collect_scores.main([
            lightdock_dir,
            "-o", out,
            "--pair", "recA:ligB",
            "--pairs", str(pairs),
        ])
        df = pd.read_csv(out, sep="\t")
        assert (df["label"] == "interaction").all()


# ---------------------------------------------------------------------------
# Unit tests for _has_nonempty_pair_context
# ---------------------------------------------------------------------------

class TestHasNonemptyPairContext:
    """Tests for the _has_nonempty_pair_context helper."""

    def _df(self, rows):
        return pd.DataFrame(rows, columns=["proteinA", "proteinB"])

    def test_returns_true_when_one_row_has_both(self):
        df = self._df([("A", "B"), ("", "")])
        assert collect_scores._has_nonempty_pair_context(df) is True

    def test_returns_false_when_no_row_has_both(self):
        """proteinA and proteinB populated on *different* rows — must return False."""
        df = self._df([("A", ""), ("", "B")])
        assert collect_scores._has_nonempty_pair_context(df) is False

    def test_returns_false_when_all_blank(self):
        df = self._df([("", ""), ("", "")])
        assert collect_scores._has_nonempty_pair_context(df) is False

    def test_returns_false_when_columns_missing(self):
        df = pd.DataFrame({"model": ["haddock"]})
        assert collect_scores._has_nonempty_pair_context(df) is False

    def test_handles_nan_values(self):
        df = pd.DataFrame({"proteinA": [None, "A"], "proteinB": ["B", None]})
        assert collect_scores._has_nonempty_pair_context(df) is False

    def test_handles_whitespace_only(self):
        df = self._df([("  ", "B"), ("A", "  ")])
        assert collect_scores._has_nonempty_pair_context(df) is False
