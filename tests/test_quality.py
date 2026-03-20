"""Tests for ppinsight.quality – DockQ-based quality assessment."""

import math
import os
import tempfile

import pandas as pd
import pytest

from ppinsight.quality import (
    CAPRI_ORDER,
    CAPRI_THRESHOLDS,
    add_quality_to_scores,
    capri_success_rate,
    capri_summary,
    classify_capri,
    evaluate_complex,
    evaluate_directory,
    _find_model_files,
)


# ---------------------------------------------------------------------------
# Helpers: fake DockQ loader / runner for fully deterministic tests
# ---------------------------------------------------------------------------

class _FakeModel:
    """Mimics a BioPython Model returned by DockQ's load_PDB."""

    def __init__(self, path: str):
        self.path = path

    def get_chains(self):
        return []


def _fake_load_pdb(path, **kw):
    return _FakeModel(path)


def _make_fake_runner(dockq: float = 0.55,
                      fnat: float = 0.6,
                      irmsd: float = 2.5,
                      lrmsd: float = 5.0):
    """Return a fake ``run_on_all_native_interfaces`` that returns
    deterministic metrics without touching real PDB files."""

    def _runner(model, native, **kw):
        result_dict = {
            "AB": {
                "DockQ": dockq,
                "F1": fnat * 0.9,
                "iRMSD": irmsd,
                "LRMSD": lrmsd,
                "fnat": fnat,
                "nat_correct": 50,
                "nat_total": 100,
                "fnonnat": 0.1,
                "nonnat_count": 5,
                "model_total": 55,
                "clashes": 2,
                "len1": 200,
                "len2": 80,
                "class1": "receptor",
                "class2": "ligand",
                "is_het": False,
                "chain1": "A",
                "chain2": "B",
                "chain_map": kw.get("chain_map", {"A": "A", "B": "B"}),
            }
        }
        import numpy as np
        return result_dict, np.float32(dockq)

    return _runner


# ---------------------------------------------------------------------------
# classify_capri
# ---------------------------------------------------------------------------

class TestClassifyCapri:
    def test_high(self):
        assert classify_capri(0.85) == "high"

    def test_boundary_high(self):
        assert classify_capri(0.80) == "high"

    def test_medium(self):
        assert classify_capri(0.55) == "medium"

    def test_boundary_medium(self):
        assert classify_capri(0.49) == "medium"

    def test_acceptable(self):
        assert classify_capri(0.30) == "acceptable"

    def test_boundary_acceptable(self):
        assert classify_capri(0.23) == "acceptable"

    def test_incorrect(self):
        assert classify_capri(0.10) == "incorrect"

    def test_zero(self):
        assert classify_capri(0.0) == "incorrect"

    def test_perfect(self):
        assert classify_capri(1.0) == "high"

    def test_just_below_acceptable(self):
        assert classify_capri(0.229) == "incorrect"

    def test_capri_order(self):
        """CAPRI_ORDER must go from best to worst."""
        assert CAPRI_ORDER == ["high", "medium", "acceptable", "incorrect"]

    def test_thresholds_descending(self):
        """Thresholds must be in descending order for the classify logic."""
        vals = list(CAPRI_THRESHOLDS.values())
        assert vals == sorted(vals, reverse=True)


# ---------------------------------------------------------------------------
# evaluate_complex
# ---------------------------------------------------------------------------

class TestEvaluateComplex:
    def test_basic(self, tmp_path):
        """evaluate_complex with injectable fakes returns expected keys."""
        model = str(tmp_path / "model.pdb")
        native = str(tmp_path / "native.pdb")
        # Create dummy files (fakes don't read them)
        for p in (model, native):
            open(p, "w").close()

        result = evaluate_complex(
            model, native,
            _load_fn=_fake_load_pdb,
            _run_fn=_make_fake_runner(dockq=0.65, fnat=0.7, irmsd=1.8, lrmsd=3.5),
        )

        assert result["DockQ"] == pytest.approx(0.65, abs=1e-4)
        assert result["fnat"] == pytest.approx(0.7, abs=1e-4)
        assert result["iRMSD"] == pytest.approx(1.8, abs=1e-4)
        assert result["LRMSD"] == pytest.approx(3.5, abs=1e-4)
        assert result["capri_class"] == "medium"
        assert result["model_path"] == model
        assert result["native_path"] == native
        assert "interfaces" in result
        assert result["n_interfaces"] == 1

    def test_high_quality(self, tmp_path):
        model = str(tmp_path / "m.pdb")
        native = str(tmp_path / "n.pdb")
        for p in (model, native):
            open(p, "w").close()

        result = evaluate_complex(
            model, native,
            _load_fn=_fake_load_pdb,
            _run_fn=_make_fake_runner(dockq=0.95),
        )
        assert result["capri_class"] == "high"

    def test_incorrect_quality(self, tmp_path):
        model = str(tmp_path / "m.pdb")
        native = str(tmp_path / "n.pdb")
        for p in (model, native):
            open(p, "w").close()

        result = evaluate_complex(
            model, native,
            _load_fn=_fake_load_pdb,
            _run_fn=_make_fake_runner(dockq=0.05),
        )
        assert result["capri_class"] == "incorrect"

    def test_chain_map_forwarded(self, tmp_path):
        """chain_map kwarg is passed through to the runner."""
        model = str(tmp_path / "m.pdb")
        native = str(tmp_path / "n.pdb")
        for p in (model, native):
            open(p, "w").close()

        received = {}

        def capturing_runner(model_struct, native_struct, **kw):
            received.update(kw)
            return _make_fake_runner()(model_struct, native_struct, **kw)

        evaluate_complex(
            model, native,
            chain_map={"X": "A", "Y": "B"},
            _load_fn=_fake_load_pdb,
            _run_fn=capturing_runner,
        )
        assert received["chain_map"] == {"X": "A", "Y": "B"}


# ---------------------------------------------------------------------------
# evaluate_directory
# ---------------------------------------------------------------------------

class TestEvaluateDirectory:
    def _make_model_dir(self, tmp_path, n=3, subdir="models"):
        """Create a directory with n dummy PDB files."""
        d = tmp_path / subdir
        d.mkdir()
        for i in range(n):
            (d / f"model_{i}.pdb").write_text(f"ATOM dummy {i}")
        native = tmp_path / "native.pdb"
        native.write_text("ATOM native")
        return str(d), str(native)

    def test_basic(self, tmp_path):
        model_dir, native_path = self._make_model_dir(tmp_path, n=3)
        df = evaluate_directory(
            model_dir, native_path,
            _load_fn=_fake_load_pdb,
            _run_fn=_make_fake_runner(dockq=0.55),
        )
        assert len(df) == 3
        assert "DockQ" in df.columns
        assert "capri_class" in df.columns
        assert all(df["capri_class"] == "medium")

    def test_empty_directory(self, tmp_path):
        d = tmp_path / "empty"
        d.mkdir()
        native = tmp_path / "native.pdb"
        native.write_text("ATOM")
        with pytest.raises(FileNotFoundError, match="No model PDB files"):
            evaluate_directory(str(d), str(native),
                               _load_fn=_fake_load_pdb,
                               _run_fn=_make_fake_runner())

    def test_error_handling(self, tmp_path):
        """Models that fail evaluation get capri_class='error'."""
        d = tmp_path / "models"
        d.mkdir()
        (d / "good.pdb").write_text("ATOM good")
        (d / "bad.pdb").write_text("ATOM bad")
        native = tmp_path / "native.pdb"
        native.write_text("ATOM native")

        call_count = [0]

        def flaky_runner(model, native, **kw):
            call_count[0] += 1
            if call_count[0] == 1:
                raise RuntimeError("parse failure")
            return _make_fake_runner(dockq=0.60)(model, native, **kw)

        df = evaluate_directory(
            str(d), str(native),
            _load_fn=_fake_load_pdb,
            _run_fn=flaky_runner,
        )
        assert len(df) == 2
        assert "error" in df["capri_class"].values
        assert "medium" in df["capri_class"].values

    def test_glob_pattern(self, tmp_path):
        """Custom glob_pattern restricts which files are found."""
        d = tmp_path / "run"
        d.mkdir()
        (d / "docked_1.pdb").write_text("ATOM")
        (d / "docked_2.pdb").write_text("ATOM")
        (d / "setup.pdb").write_text("ATOM")  # should be excluded
        native = tmp_path / "native.pdb"
        native.write_text("ATOM")

        df = evaluate_directory(
            str(d), str(native),
            glob_pattern="docked_*.pdb",
            _load_fn=_fake_load_pdb,
            _run_fn=_make_fake_runner(),
        )
        assert len(df) == 2


# ---------------------------------------------------------------------------
# _find_model_files
# ---------------------------------------------------------------------------

class TestFindModelFiles:
    def test_lightdock_layout(self, tmp_path):
        for i in range(3):
            swarm = tmp_path / f"swarm_{i}"
            swarm.mkdir()
            (swarm / f"lightdock_{i}.pdb").write_text("ATOM")
        # Also need a setup.json so detect_engine can find it, but we
        # bypass detection by passing engine explicitly
        files = _find_model_files(str(tmp_path), engine="lightdock")
        assert len(files) == 3

    def test_rosetta_layout(self, tmp_path):
        out = tmp_path / "output_files"
        out.mkdir()
        (out / "decoy_001.pdb").write_text("ATOM")
        (out / "decoy_002.pdb").write_text("ATOM")
        files = _find_model_files(str(tmp_path), engine="rosetta")
        assert len(files) == 2

    def test_generic_fallback(self, tmp_path):
        (tmp_path / "a.pdb").write_text("ATOM")
        sub = tmp_path / "sub"
        sub.mkdir()
        (sub / "b.pdb").write_text("ATOM")
        files = _find_model_files(str(tmp_path), engine=None)
        assert len(files) == 2

    def test_custom_glob(self, tmp_path):
        (tmp_path / "model_1.pdb").write_text("ATOM")
        (tmp_path / "model_2.pdb").write_text("ATOM")
        (tmp_path / "native.pdb").write_text("ATOM")
        files = _find_model_files(str(tmp_path), glob_pattern="model_*.pdb")
        assert len(files) == 2


# ---------------------------------------------------------------------------
# capri_summary / capri_success_rate
# ---------------------------------------------------------------------------

class TestCapriSummary:
    def _make_df(self, classes):
        return pd.DataFrame({"capri_class": classes, "DockQ": range(len(classes))})

    def test_basic(self):
        df = self._make_df(["high", "medium", "medium", "incorrect", "acceptable"])
        counts = capri_summary(df)
        assert counts["high"] == 1
        assert counts["medium"] == 2
        assert counts["acceptable"] == 1
        assert counts["incorrect"] == 1

    def test_all_categories_present(self):
        df = self._make_df(["high"])
        counts = capri_summary(df)
        # All four standard categories should be present
        for cat in CAPRI_ORDER:
            assert cat in counts

    def test_success_rate(self):
        df = self._make_df(["high", "medium", "acceptable", "incorrect", "incorrect"])
        rate = capri_success_rate(df)
        assert rate == pytest.approx(3 / 5)

    def test_success_rate_perfect(self):
        df = self._make_df(["high", "high", "medium"])
        assert capri_success_rate(df) == pytest.approx(1.0)

    def test_success_rate_zero(self):
        df = self._make_df(["incorrect", "incorrect"])
        assert capri_success_rate(df) == pytest.approx(0.0)

    def test_success_rate_empty(self):
        df = pd.DataFrame({"capri_class": []})
        assert math.isnan(capri_success_rate(df))

    def test_with_errors(self):
        df = self._make_df(["high", "error", "incorrect"])
        counts = capri_summary(df)
        assert counts["high"] == 1
        assert counts.get("error", 0) == 1


# ---------------------------------------------------------------------------
# add_quality_to_scores
# ---------------------------------------------------------------------------

class TestAddQualityToScores:
    def test_basic(self):
        existing = pd.DataFrame({
            "model": ["ldock"], "score_type": ["luciferin_score"],
            "score_value": [42.0], "proteinA": ["A"], "proteinB": ["B"],
        })
        quality = pd.DataFrame({
            "model_path": ["m.pdb"],
            "DockQ": [0.65],
            "fnat": [0.7],
            "iRMSD": [1.8],
            "LRMSD": [3.5],
            "capri_class": ["medium"],
        })
        result = add_quality_to_scores(existing, quality, model_label="ldock")

        # Original row + 4 quality rows
        assert len(result) == 5
        quality_rows = result[result["score_type"].str.startswith("quality_")]
        assert len(quality_rows) == 4
        assert set(quality_rows["score_type"]) == {
            "quality_DockQ", "quality_fnat", "quality_iRMSD", "quality_LRMSD",
        }

    def test_empty_quality(self):
        existing = pd.DataFrame({
            "model": ["x"], "score_type": ["s"],
            "score_value": [1.0], "proteinA": [""], "proteinB": [""],
        })
        quality = pd.DataFrame(columns=["model_path", "DockQ", "fnat", "iRMSD", "LRMSD"])
        result = add_quality_to_scores(existing, quality)
        assert len(result) == 1  # no quality rows added


# ---------------------------------------------------------------------------
# METRIC_METADATA integration
# ---------------------------------------------------------------------------

class TestMetricMetadataIntegration:
    def test_quality_metrics_registered(self):
        """DockQ quality metrics should be in METRIC_METADATA."""
        from ppinsight.visualizer import METRIC_METADATA
        # These are the metrics we add in the visualizer update
        expected = {"quality_dockq", "quality_fnat", "quality_irmsd", "quality_lrmsd"}
        for m in expected:
            assert m in METRIC_METADATA, f"{m} not in METRIC_METADATA"

    def test_quality_dockq_higher_is_better(self):
        from ppinsight.visualizer import get_metric_direction
        assert get_metric_direction("quality_dockq") is True

    def test_quality_irmsd_lower_is_better(self):
        from ppinsight.visualizer import get_metric_direction
        assert get_metric_direction("quality_irmsd") is False


# ---------------------------------------------------------------------------
# CLI smoke test
# ---------------------------------------------------------------------------

class TestCLI:
    def test_single_model(self, tmp_path, monkeypatch):
        """CLI in single-model mode prints DockQ results."""
        model = tmp_path / "model.pdb"
        native = tmp_path / "native.pdb"
        model.write_text("ATOM dummy")
        native.write_text("ATOM native")

        # Monkeypatch the DockQ functions to avoid needing real PDBs
        import ppinsight.quality as qmod

        monkeypatch.setattr(qmod, "load_PDB", _fake_load_pdb)
        orig_eval = qmod.evaluate_complex

        def patched_eval(mp, np_, **kw):
            return orig_eval(mp, np_, _load_fn=_fake_load_pdb,
                             _run_fn=_make_fake_runner(dockq=0.72), **kw)

        monkeypatch.setattr(qmod, "evaluate_complex", patched_eval)

        from ppinsight.quality import main
        # Should not raise
        main([str(model), str(native)])

    def test_directory_mode(self, tmp_path, monkeypatch):
        """CLI in directory mode with --summary."""
        d = tmp_path / "run"
        d.mkdir()
        for i in range(3):
            (d / f"model_{i}.pdb").write_text("ATOM")
        native = tmp_path / "native.pdb"
        native.write_text("ATOM")
        output = tmp_path / "quality.csv"

        import ppinsight.quality as qmod

        monkeypatch.setattr(qmod, "load_PDB", _fake_load_pdb)
        orig_eval_dir = qmod.evaluate_directory

        def patched_eval_dir(rd, np_, **kw):
            return orig_eval_dir(rd, np_, _load_fn=_fake_load_pdb,
                                 _run_fn=_make_fake_runner(dockq=0.55), **kw)

        monkeypatch.setattr(qmod, "evaluate_directory", patched_eval_dir)

        from ppinsight.quality import main
        main([str(d), str(native), "-o", str(output), "--summary"])
        assert output.exists()
        df = pd.read_csv(str(output))
        assert len(df) == 3


# ---------------------------------------------------------------------------
# DockQ-unavailable graceful degradation
# ---------------------------------------------------------------------------

class TestDockQUnavailable:
    """Verify that the module degrades gracefully when DockQ is missing."""

    def test_classify_capri_works_without_dockq(self, monkeypatch):
        """CAPRI helpers don't need DockQ at all."""
        import ppinsight.quality as qmod
        monkeypatch.setattr(qmod, "_DOCKQ_AVAILABLE", False)
        # classify_capri is pure Python – should still work
        assert qmod.classify_capri(0.9) == "high"

    def test_evaluate_complex_raises_without_dockq(self, tmp_path, monkeypatch):
        """evaluate_complex raises ImportError when DockQ is missing
        and no injectable fakes are provided."""
        import ppinsight.quality as qmod
        monkeypatch.setattr(qmod, "_DOCKQ_AVAILABLE", False)
        monkeypatch.setattr(
            qmod, "_DOCKQ_IMPORT_ERROR",
            "DockQ is not installed. Install with: pip install 'ppinsight[quality]'"
        )

        model = tmp_path / "m.pdb"
        native = tmp_path / "n.pdb"
        model.write_text("ATOM")
        native.write_text("ATOM")

        with pytest.raises(ImportError, match="ppinsight\\[quality\\]"):
            qmod.evaluate_complex(str(model), str(native))

    def test_evaluate_complex_works_with_fakes_even_without_dockq(
        self, tmp_path, monkeypatch
    ):
        """Injectable _load_fn/_run_fn bypass the DockQ requirement."""
        import ppinsight.quality as qmod
        monkeypatch.setattr(qmod, "_DOCKQ_AVAILABLE", False)

        model = tmp_path / "m.pdb"
        native = tmp_path / "n.pdb"
        model.write_text("ATOM")
        native.write_text("ATOM")

        # Both fakes provided → should NOT check for DockQ
        result = qmod.evaluate_complex(
            str(model), str(native),
            _load_fn=_fake_load_pdb,
            _run_fn=_make_fake_runner(dockq=0.50),
        )
        assert result["DockQ"] == pytest.approx(0.50)

    def test_cli_exits_with_message_without_dockq(self, monkeypatch, capsys):
        """CLI main() exits with code 1 and a helpful message."""
        import ppinsight.quality as qmod
        monkeypatch.setattr(qmod, "_DOCKQ_AVAILABLE", False)
        monkeypatch.setattr(
            qmod, "_DOCKQ_IMPORT_ERROR",
            "DockQ is not installed. Install with: pip install 'ppinsight[quality]'"
        )

        with pytest.raises(SystemExit) as exc_info:
            qmod.main(["dummy.pdb", "native.pdb"])

        assert exc_info.value.code == 1
        captured = capsys.readouterr()
        assert "ppinsight[quality]" in captured.err
