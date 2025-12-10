import sys
import os

# Ensure the package under src/ is importable during tests
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)


import pytest


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
