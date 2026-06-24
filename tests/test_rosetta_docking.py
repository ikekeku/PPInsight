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

from ppinsight.docking import DockingPipeline  # noqa: E402
from ppinsight.rosetta.prepare_structure import prepare_structures  # noqa: E402


def _make_atom_line(serial, atom_name, residue_name, chain_id, residue_number, x):
    return (
        f"ATOM  {serial:5d} {atom_name:>4} {residue_name:>3} "
        f"{chain_id}{residue_number:4d}    "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.00:6.2f}"
        f"{20.00:6.2f}          {atom_name[0]:>2}\n"
    )


def _write_single_residue_chain(
    handle,
    start_serial,
    residue_name,
    chain_id,
    residue_number,
    x,
):
    handle.write(
        _make_atom_line(start_serial, "N", residue_name, chain_id, residue_number, x)
        + _make_atom_line(
            start_serial + 1,
            "CA",
            residue_name,
            chain_id,
            residue_number,
            x + 1.5,
        )
        + _make_atom_line(
            start_serial + 2,
            "C",
            residue_name,
            chain_id,
            residue_number,
            x + 2.9,
        )
        + _make_atom_line(
            start_serial + 3,
            "O",
            residue_name,
            chain_id,
            residue_number,
            x + 3.9,
        )
        + (
            f"TER   {start_serial + 4:5d}      {residue_name:>3} "
            f"{chain_id}{residue_number:4d}\n"
        )
    )

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


def test_prepare_structures_filters_dbref_mixed_complex_for_rosetta(tmp_path):
    receptor = tmp_path / "P11111.pdb"
    ligand = tmp_path / "P22222.pdb"

    receptor.write_text(
        "DBREF  1ABC A    1     1  UNP    P11111   REC_HUMAN       1      1\n"
        + _make_atom_line(1, "N", "GLY", "A", 1, 0.0)
        + _make_atom_line(2, "CA", "GLY", "A", 1, 1.5)
        + _make_atom_line(3, "C", "GLY", "A", 1, 2.9)
        + _make_atom_line(4, "O", "GLY", "A", 1, 3.9)
        + "TER       5      GLY A   1\n"
        + "END\n",
        encoding="utf-8",
    )

    with ligand.open("w", encoding="utf-8") as handle:
        handle.write(
            "DBREF  2ABC V    1     1  UNP    P22222   LIG_HUMAN       1      1\n"
            "DBREF  2ABC W    1     1  UNP    P22222   LIG_HUMAN       2      2\n"
            "DBREF  2ABC X    1     1  UNP    Q99999   OTHER_HUMAN     1      1\n"
        )
        _write_single_residue_chain(handle, 1, "SER", "V", 1, 10.0)
        _write_single_residue_chain(handle, 6, "TYR", "W", 2, 20.0)
        _write_single_residue_chain(handle, 11, "ALA", "X", 3, 30.0)
        handle.write(
            "HETATM   16  S   SO4 Y   4      40.000   0.000   0.000"
            "  1.00 20.00           S\n"
            "END\n"
        )

    pose = prepare_structures(
        str(receptor),
        str(ligand),
        relax=False,
        verbose=False,
    )

    assert pose.total_residue() == 3
    assert pose.num_chains() == 3
    assert [
        pose.pdb_info().chain(pose.chain_begin(chain_index))
        for chain_index in range(1, pose.num_chains() + 1)
    ] == ["A", "B", "C"]


def test_prepare_structures_can_disable_dbref_auto_filter(tmp_path):
    receptor = tmp_path / "P11111.pdb"
    ligand = tmp_path / "P22222.pdb"

    receptor.write_text(
        "DBREF  1ABC A    1     1  UNP    P11111   REC_HUMAN       1      1\n"
        + _make_atom_line(1, "N", "GLY", "A", 1, 0.0)
        + _make_atom_line(2, "CA", "GLY", "A", 1, 1.5)
        + _make_atom_line(3, "C", "GLY", "A", 1, 2.9)
        + _make_atom_line(4, "O", "GLY", "A", 1, 3.9)
        + "TER       5      GLY A   1\n"
        + "END\n",
        encoding="utf-8",
    )

    with ligand.open("w", encoding="utf-8") as handle:
        handle.write(
            "DBREF  2ABC V    1     1  UNP    P22222   LIG_HUMAN       1      1\n"
            "DBREF  2ABC W    1     1  UNP    P22222   LIG_HUMAN       2      2\n"
            "DBREF  2ABC X    1     1  UNP    Q99999   OTHER_HUMAN     1      1\n"
        )
        _write_single_residue_chain(handle, 1, "SER", "V", 1, 10.0)
        _write_single_residue_chain(handle, 6, "TYR", "W", 2, 20.0)
        _write_single_residue_chain(handle, 11, "ALA", "X", 3, 30.0)
        handle.write(
            "HETATM   16  S   SO4 Y   4      40.000   0.000   0.000"
            "  1.00 20.00           S\n"
            "END\n"
        )

    pose = prepare_structures(
        str(receptor),
        str(ligand),
        relax=False,
        verbose=False,
        auto_filter=False,
    )

    assert pose.total_residue() == 4
    assert pose.num_chains() == 4
