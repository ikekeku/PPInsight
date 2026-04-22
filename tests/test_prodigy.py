"""Tests for ppinsight.prodigy."""

from __future__ import annotations

import math
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

try:
    from prodigy_prot.modules.prodigy import Prodigy as _ProdigyClass  # noqa: F401
    PRODIGY_AVAILABLE = True
except ImportError:
    try:
        from prodigy_prot.predict_IC import Prodigy as _ProdigyClass  # noqa: F401

        PRODIGY_AVAILABLE = True
    except ImportError:
        PRODIGY_AVAILABLE = False

requires_prodigy = pytest.mark.skipif(
    not PRODIGY_AVAILABLE,
    reason="prodigy-prot not installed (pip install ppinsight[prodigy])",
)

from ppinsight.prodigy import (  # noqa: E402
    _require_prodigy,
    add_prodigy_to_scores,
    score_directory,
    score_pdb,
)


class TestRequireProdigy:
    @requires_prodigy
    def test_returns_class_when_available(self):
        cls = _require_prodigy()
        assert cls is not None
        assert callable(cls)

    def test_raises_import_error_when_missing(self):
        with patch.dict("sys.modules", {"prodigy_prot": None,
                                        "prodigy_prot.modules": None,
                                        "prodigy_prot.modules.prodigy": None,
                                        "prodigy_prot.predict_IC": None}):
            with pytest.raises(ImportError, match="prodigy-prot"):
                _require_prodigy()


class TestScorePdb:
    def test_missing_file_returns_nan(self, tmp_path):
        result = score_pdb(tmp_path / "nonexistent.pdb")
        assert math.isnan(result["prodigy_ddg"])
        assert math.isnan(result["prodigy_kd"])
        assert "file_not_found" in result.get("error", "")

    @requires_prodigy
    def test_real_pdb_file(self):
        pdb = (
            Path(__file__).resolve().parent.parent
            / "examples" / "ppinsight_data" / "input_files" / "2UUY_rec.pdb"
        )
        if not pdb.exists():
            pytest.skip(f"Example PDB not found: {pdb}")
        result = score_pdb(pdb)
        assert "prodigy_ddg" in result
        assert "prodigy_kd" in result

    def test_mocked_successful_scoring(self, tmp_path):
        pdb_file = tmp_path / "complex.pdb"
        pdb_file.write_text("ATOM  …\n")

        mock_runner = MagicMock()
        mock_runner.ba_val = -9.5
        mock_runner.kd_val = 1.2e-7
        mock_runner.nis_a = 0.32
        mock_runner.nis_c = 0.18
        mock_runner.bins = {"CC": 5, "CP": 3, "AC": 2, "AA": 1, "PP": 0, "AP": 1}

        MockProdigy = MagicMock(return_value=mock_runner)
        mock_chain = MagicMock()
        mock_chain.id = "A"
        mock_model = MagicMock()
        mock_model.get_chains.return_value = [mock_chain]
        mock_struct = MagicMock()
        mock_struct.__iter__ = MagicMock(return_value=iter([mock_model]))

        with patch("ppinsight.prodigy._require_prodigy", return_value=MockProdigy), \
             patch("ppinsight.prodigy._parse_pdb", return_value=mock_struct):
            result = score_pdb(pdb_file)

        assert result["prodigy_ddg"] == pytest.approx(-9.5)
        assert result["prodigy_kd"] == pytest.approx(1.2e-7)
        assert result["nis_a"] == pytest.approx(0.32)
        assert result["n_contacts"] == 12

    def test_no_contacts_returns_nan(self, tmp_path):
        pdb_file = tmp_path / "single_chain.pdb"
        pdb_file.write_text("ATOM  …\n")

        mock_runner = MagicMock()
        mock_runner.predict.side_effect = ValueError(
            "No contacts found for selection"
        )
        MockProdigy = MagicMock(return_value=mock_runner)
        mock_chain = MagicMock()
        mock_chain.id = "A"
        mock_model = MagicMock()
        mock_model.get_chains.return_value = [mock_chain]
        mock_struct = MagicMock()
        mock_struct.__iter__ = MagicMock(return_value=iter([mock_model]))

        with patch("ppinsight.prodigy._require_prodigy", return_value=MockProdigy), \
             patch("ppinsight.prodigy._parse_pdb", return_value=mock_struct):
            result = score_pdb(pdb_file)

        assert math.isnan(result["prodigy_ddg"])
        assert math.isnan(result["prodigy_kd"])
        assert result["error"] == "no_contacts"


class TestScoreDirectory:
    def test_empty_directory_returns_empty_df(self, tmp_path):
        with pytest.warns(UserWarning, match="No files"):
            result = score_directory(tmp_path)
        assert result.empty

    def test_returns_dataframe_with_pdb_column(self, tmp_path):
        (tmp_path / "a.pdb").write_text("ATOM  …\n")
        (tmp_path / "b.pdb").write_text("ATOM  …\n")

        with patch("ppinsight.prodigy.score_pdb") as mock_score:
            mock_score.return_value = {
                "prodigy_ddg": -9.0, "prodigy_kd": 1e-7,
                "nis_a": 0.3, "nis_c": 0.2, "n_contacts": 10,
            }
            result = score_directory(tmp_path)

        assert "pdb" in result.columns
        assert len(result) == 2
        assert set(result["pdb"]) == {"a.pdb", "b.pdb"}


class TestAddProdigyToScores:
    def _make_scores(self) -> pd.DataFrame:
        return pd.DataFrame({
            "model": ["HADDOCK"] * 4 + ["LightDock"] * 4,
            "score_type": ["score"] * 8,
            "score_value": [-100, -90, -80, -70, 10, 12, 14, 16],
            "proteinA": ["rec"] * 8,
            "proteinB": ["lig"] * 8,
            "pdb": ["h1.pdb", "h2.pdb", "h3.pdb", "h4.pdb",
                    "l1.pdb", "l2.pdb", "l3.pdb", "l4.pdb"],
        })

    def test_adds_prodigy_columns(self, tmp_path):
        scores = self._make_scores()
        for pdb in scores["pdb"].unique():
            (tmp_path / pdb).write_text("ATOM  …\n")

        with patch("ppinsight.prodigy.score_pdb") as mock_score:
            mock_score.return_value = {"prodigy_ddg": -8.5, "prodigy_kd": 5e-7}
            result = add_prodigy_to_scores(scores, pdb_dir=tmp_path)

        assert "prodigy_ddg" in result.columns
        assert "prodigy_kd" in result.columns
        assert not result["prodigy_ddg"].isna().any()
        assert (result["prodigy_ddg"] - (-8.5)).abs().max() < 1e-6

    def test_engine_filter(self, tmp_path):
        scores = self._make_scores()
        for pdb in scores["pdb"].unique():
            (tmp_path / pdb).write_text("ATOM  …\n")

        with patch("ppinsight.prodigy.score_pdb") as mock_score:
            mock_score.return_value = {"prodigy_ddg": -9.0, "prodigy_kd": 1e-7}
            result = add_prodigy_to_scores(
                scores, pdb_dir=tmp_path, engine="HADDOCK"
            )

        haddock = result[result["model"] == "HADDOCK"]
        lightdock = result[result["model"] == "LightDock"]
        assert not haddock["prodigy_ddg"].isna().any()
        assert lightdock["prodigy_ddg"].isna().all()

    def test_no_pdb_column_warns(self, tmp_path):
        scores = self._make_scores().drop(columns=["pdb"])
        with pytest.warns(UserWarning, match="'pdb' column"):
            result = add_prodigy_to_scores(scores, pdb_dir=tmp_path)
        assert result["prodigy_ddg"].isna().all()

    def test_top_n_requires_metric(self, tmp_path):
        scores = self._make_scores()
        with pytest.raises(ValueError, match="metric must be specified"):
            add_prodigy_to_scores(scores, pdb_dir=tmp_path, top_n=2)

    def test_top_n_lower_is_better_uses_nsmallest(self, tmp_path):
        """For lower-is-better metrics (e.g. HADDOCK 'score'), nsmallest selects best."""
        scores = pd.DataFrame({
            "model": ["HADDOCK"] * 4,
            "score_type": ["score"] * 4,   # "score" → higher_is_better=False in METRIC_METADATA
            "score_value": [-100.0, -90.0, -80.0, -70.0],
            "proteinA": ["rec"] * 4,
            "proteinB": ["lig"] * 4,
            "pdb": ["h1.pdb", "h2.pdb", "h3.pdb", "h4.pdb"],
        })
        for pdb in scores["pdb"]:
            (tmp_path / pdb).write_text("ATOM  …\n")

        with patch("ppinsight.prodigy.score_pdb") as mock_score:
            mock_score.return_value = {"prodigy_ddg": -8.5, "prodigy_kd": 5e-7}
            result = add_prodigy_to_scores(
                scores, pdb_dir=tmp_path, top_n=2, metric="score"
            )

        # Only the 2 poses with the lowest (best) scores should be scored.
        scored = result[~result["prodigy_ddg"].isna()]
        assert len(scored) == 2
        assert set(scored["score_value"]) == {-100.0, -90.0}

    def test_top_n_higher_is_better_uses_nlargest(self, tmp_path):
        """For higher-is-better metrics (e.g. luciferin_score), nlargest selects best."""
        scores = pd.DataFrame({
            "model": ["LightDock"] * 4,
            "score_type": ["luciferin_score"] * 4,
            "score_value": [10.0, 12.0, 14.0, 16.0],
            "proteinA": ["rec"] * 4,
            "proteinB": ["lig"] * 4,
            "pdb": ["l1.pdb", "l2.pdb", "l3.pdb", "l4.pdb"],
        })
        for pdb in scores["pdb"]:
            (tmp_path / pdb).write_text("ATOM  …\n")

        with patch("ppinsight.prodigy.score_pdb") as mock_score:
            mock_score.return_value = {"prodigy_ddg": -7.0, "prodigy_kd": 1e-6}
            result = add_prodigy_to_scores(
                scores, pdb_dir=tmp_path, top_n=2, metric="luciferin_score"
            )

        # Only the 2 poses with the highest (best) scores should be scored.
        scored = result[~result["prodigy_ddg"].isna()]
        assert len(scored) == 2
        assert set(scored["score_value"]) == {14.0, 16.0}

    def test_top_n_filters_by_metric_in_long_format(self, tmp_path):
        """Ranking uses only the specified score_type rows, not all long-format rows."""
        scores = pd.DataFrame({
            "model": ["HADDOCK"] * 8,
            # 4 poses each with two score_type rows
            "score_type": ["score", "clash_score"] * 4,
            "score_value": [-100.0, 5.0, -90.0, 3.0, -80.0, 8.0, -70.0, 2.0],
            "proteinA": ["rec"] * 8,
            "proteinB": ["lig"] * 8,
            "pdb": ["h1.pdb", "h1.pdb", "h2.pdb", "h2.pdb",
                    "h3.pdb", "h3.pdb", "h4.pdb", "h4.pdb"],
        })
        for pdb in ["h1.pdb", "h2.pdb", "h3.pdb", "h4.pdb"]:
            (tmp_path / pdb).write_text("ATOM  …\n")

        with patch("ppinsight.prodigy.score_pdb") as mock_score:
            mock_score.return_value = {"prodigy_ddg": -8.0, "prodigy_kd": 1e-7}
            result = add_prodigy_to_scores(
                scores, pdb_dir=tmp_path, top_n=2, metric="score"
            )

        # Best 2 by HADDOCK 'score' (lower-is-better) → h1 (-100), h2 (-90)
        scored = result[~result["prodigy_ddg"].isna()]
        scored_pdbs = set(scored["pdb"])
        assert scored_pdbs == {"h1.pdb", "h2.pdb"}
        # All rows (both score_types) for those PDFs get prodigy results
        assert len(scored) == 4  # 2 score_types × 2 PDFs

