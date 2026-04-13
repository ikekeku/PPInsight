import os
import sys

# Ensure the package under src/ is importable during tests
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)


import pytest  # noqa: E402


def pytest_collection_modifyitems(config, items):
    """Auto-skip tests that require unavailable external tools.

    - ``requires_rosetta``: skipped unless ``ROSETTA_AVAILABLE=1`` is set
      *and* PyRosetta can actually be imported.
    - ``integration``: skipped unless ``RUN_INTEGRATION_TESTS=1`` is set.
    """
    rosetta_available = os.environ.get("ROSETTA_AVAILABLE", "0") == "1"
    run_integration = os.environ.get("RUN_INTEGRATION_TESTS", "0") == "1"

    skip_rosetta = pytest.mark.skip(
        reason="PyRosetta not available; set ROSETTA_AVAILABLE=1 to run"
    )
    skip_integration = pytest.mark.skip(
        reason="Integration tests disabled; set RUN_INTEGRATION_TESTS=1 to run"
    )

    for item in items:
        if item.get_closest_marker("requires_rosetta") and not rosetta_available:
            item.add_marker(skip_rosetta)
        if item.get_closest_marker("integration") and not run_integration:
            item.add_marker(skip_integration)


@pytest.fixture
def sample_input_dirs(tmp_path):
    """Create a small input and output directory structure for tests.

    Returns a dict with keys: inp, rec, lig, work_root
    """
    inp = tmp_path / "input"
    inp.mkdir()
    rec = inp / "rec.pdb"
    lig = inp / "lig.pdb"
    rec.write_text("ATOM\n")
    lig.write_text("ATOM\n")

    work_root = tmp_path / "output"
    work_root.mkdir()

    return {
        "inp": inp,
        "rec": rec,
        "lig": lig,
        "work_root": work_root,
    }


# --------------- PDB fixtures for rosetta tests ---------------

EXAMPLE_PDBS = os.path.join(ROOT, "examples", "ppinsight_data", "input_files")


@pytest.fixture
def pdb_rec():
    """Path to the receptor PDB bundled in the repo."""
    p = os.path.join(EXAMPLE_PDBS, "2UUY_rec.pdb")
    if not os.path.isfile(p):
        pytest.skip("Example PDB 2UUY_rec.pdb not found")
    return p


@pytest.fixture
def pdb_lig():
    """Path to the ligand PDB bundled in the repo."""
    p = os.path.join(EXAMPLE_PDBS, "2UUY_lig.pdb")
    if not os.path.isfile(p):
        pytest.skip("Example PDB 2UUY_lig.pdb not found")
    return p


# --------------- Dummy score CSV fixtures for visualizer tests ---------------

@pytest.fixture
def dummy_score_csvs(tmp_path):
    """Create three tiny CSV files that mimic docking score output.

    Returns a list of three pathlib.Path objects.
    """
    import csv

    paths = []
    for i in range(1, 4):
        p = tmp_path / f"model_{i}.csv"
        with open(p, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["score_a", "score_b"])
            writer.writerow([-(10 + i), 0.5 * i])
        paths.append(str(p))
    return paths
