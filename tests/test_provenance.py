"""Tests for the provenance module (run-level metadata tracking)."""

import datetime
import json
import os

from ppinsight.provenance import (
    extract_run_metadata,
    make_run_id,
    read_sidecar,
    sidecar_path,
    write_sidecar,
)

# ---------------------------------------------------------------------------
# make_run_id
# ---------------------------------------------------------------------------

class TestMakeRunId:
    def test_with_pair(self):
        ts = datetime.datetime(2026, 4, 2, 15, 23, 0)
        rid = make_run_id("lightdock", ("2UUY_rec", "2UUY_lig"), timestamp=ts)
        assert rid == "lightdock_2UUY_rec_2UUY_lig_20260402T152300"

    def test_without_pair(self):
        ts = datetime.datetime(2026, 1, 1, 0, 0, 0)
        rid = make_run_id("haddock", timestamp=ts)
        assert rid == "haddock_20260101T000000"

    def test_auto_timestamp(self):
        """Without an explicit timestamp, the run_id still generates."""
        rid = make_run_id("rosetta", ("A", "B"))
        assert rid.startswith("rosetta_A_B_")
        # Should contain a date-like suffix
        assert len(rid.split("_")) >= 4


# ---------------------------------------------------------------------------
# extract_run_metadata
# ---------------------------------------------------------------------------

class TestExtractMetadata:
    def test_lightdock_with_setup_json(self, tmp_path):
        setup = {"glowworms": 200, "swarms": 50, "receptor_pdb": "rec.pdb"}
        (tmp_path / "setup.json").write_text(json.dumps(setup))
        meta = extract_run_metadata("lightdock", str(tmp_path))
        assert meta["engine"] == "lightdock"
        assert meta["source_dir"] == str(tmp_path.resolve())
        assert "collected_at" in meta
        assert "hostname" in meta
        assert meta["engine_meta"]["glowworms"] == 200

    def test_lightdock_no_setup(self, tmp_path):
        meta = extract_run_metadata("lightdock", str(tmp_path))
        assert "note" in meta["engine_meta"]

    def test_haddock_with_log(self, tmp_path):
        log_text = (
            "[2025-11-17 12:54:06,096 cli INFO] \n"
            "Starting HADDOCK3 v2025.9.1 on 2025-11-17 12:54:00\n"
            "Python 3.9.25 (main, Nov  3 2025, 22:33:05)\n"
            "[2025-11-17 12:54:49,851 libworkflow INFO] "
            "Reading instructions step 0_topoaa\n"
            "[2025-11-17 12:54:49,851 libworkflow INFO] "
            "Reading instructions step 1_rigidbody\n"
        )
        (tmp_path / "log").write_text(log_text)
        meta = extract_run_metadata("haddock", str(tmp_path))
        assert meta["engine_meta"]["haddock_version"] == "v2025.9.1"
        assert meta["engine_meta"]["run_started"] == "2025-11-17 12:54:00"
        assert "0_topoaa" in meta["engine_meta"]["steps"]
        assert "1_rigidbody" in meta["engine_meta"]["steps"]

    def test_haddock_no_log(self, tmp_path):
        meta = extract_run_metadata("haddock", str(tmp_path))
        assert "note" in meta["engine_meta"]

    def test_rosetta_with_flags(self, tmp_path):
        (tmp_path / "flag_global_docking").write_text(
            "-nstruct 100\n-partners A_B\n"
        )
        meta = extract_run_metadata("rosetta", str(tmp_path))
        assert "flag_global_docking" in meta["engine_meta"]
        assert "-nstruct 100" in meta["engine_meta"]["flag_global_docking"]

    def test_unknown_engine(self, tmp_path):
        meta = extract_run_metadata("unknown_engine", str(tmp_path))
        assert meta["engine"] == "unknown_engine"
        assert meta["engine_meta"] == {}


# ---------------------------------------------------------------------------
# Sidecar I/O
# ---------------------------------------------------------------------------

class TestSidecar:
    def test_path(self):
        assert sidecar_path("scores.tsv") == "scores.tsv.provenance.json"

    def test_write_and_read(self, tmp_path):
        scores_file = str(tmp_path / "scores.tsv")
        prov = {
            "run_001": {"engine": "lightdock", "source_dir": "/foo"},
            "run_002": {"engine": "haddock", "source_dir": "/bar"},
        }
        out = write_sidecar(scores_file, prov)
        assert os.path.isfile(out)

        loaded = read_sidecar(scores_file)
        assert loaded["run_001"]["engine"] == "lightdock"
        assert loaded["run_002"]["engine"] == "haddock"

    def test_merge_preserves_old(self, tmp_path):
        """Writing new entries should not overwrite existing ones."""
        scores_file = str(tmp_path / "scores.tsv")
        write_sidecar(scores_file, {"old": {"engine": "lightdock"}})
        write_sidecar(scores_file, {"new": {"engine": "haddock"}})

        loaded = read_sidecar(scores_file)
        assert "old" in loaded
        assert "new" in loaded

    def test_read_missing(self, tmp_path):
        loaded = read_sidecar(str(tmp_path / "nonexistent.tsv"))
        assert loaded == {}
