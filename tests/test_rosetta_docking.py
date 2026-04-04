"""Tests for the Rosetta docking pipeline.

These tests require PyRosetta. If it is not installed, or if the
``ROSETTA_AVAILABLE`` environment variable is not set to ``1``, the
entire module is skipped automatically.

To run these tests locally with PyRosetta installed::

    ROSETTA_AVAILABLE=1 pytest tests/test_rosetta_docking.py
"""

import pytest

pytestmark = pytest.mark.requires_rosetta

pyrosetta = pytest.importorskip("pyrosetta", reason="PyRosetta not installed")

from ppinsight.docking import DockingPipeline


# ── edge test (no real files needed) ──────────────────────────────

def test_edge_cases():
    """n_runs < 1 must raise ValueError."""
    with pytest.raises(ValueError, match="n_runs must be at least 1"):
        DockingPipeline("protein1.pdb", "protein2.pdb", n_runs=0)


# ── smoke test: structures load and pipeline starts ───────────────

def test_smoke(pdb_rec, pdb_lig):
    """Pipeline initialises and prepare() loads real PDBs."""
    pipeline = DockingPipeline(pdb_rec, pdb_lig, n_runs=1, verbose=False)
    pose = pipeline.prepare()
    assert pose is not None


# ── one-shot test: a single docking run completes ────────────────

def test_oneshot(pdb_rec, pdb_lig):
    """A single-run pipeline executes end-to-end."""
    pipeline = DockingPipeline(pdb_rec, pdb_lig, n_runs=1, top_n=1, verbose=False)
    result = pipeline.run()
    assert "final_score" in result
    # Score should be a finite number
    assert isinstance(result["final_score"], (int, float))


# ── pattern test: prepare works for various n_runs ───────────────

def test_pattern(pdb_rec, pdb_lig):
    """prepare() succeeds regardless of the n_runs value."""
    for n_runs in [1, 2, 3]:
        pipeline = DockingPipeline(pdb_rec, pdb_lig, n_runs=n_runs, verbose=False)
        pose = pipeline.prepare()
        assert pose is not None