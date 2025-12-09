"""
FORMATTED FOR HW3!
Tests for `pdb_to_haddock.py`.

The tests monkeypatch `run_command` so haddock3 is not actually invoked.
"""
from pathlib import Path
from ppinsight import pdb_to_haddock
import pytest


class CmdRecorder:
    def __init__(self):
        self.calls = []

    def __call__(self, cmd, cwd=None):
        self.calls.append((list(cmd), cwd))


def test_smoke_haddock_pipeline_runs(tmp_path, monkeypatch):
    """
    author: ikekeku
    reviewer: fmclary
    category: smoke test

    Simple smoke test to ensure haddock_pipeline runs end-to-end without error,
    creating expected files/folders.
    """
    inp = tmp_path / "input"
    inp.mkdir()
    rec = inp / "rec.pdb"
    lig = inp / "lig.pdb"
    rec.write_text("ATOM\n")
    lig.write_text("ATOM\n")

    work_root = tmp_path / "output"
    work_root.mkdir()

    recorder = CmdRecorder()
    monkeypatch.setattr(pdb_to_haddock, "run_command", recorder)

    # Run pipeline (do not actually run haddock)
    run_dir, cfg_path, _ = pdb_to_haddock.haddock_pipeline(str(rec), str(lig), runname="smoke", run_haddock=False, base_root=str(work_root), method="haddock_runs")

    # Basic smoke assertions
    assert (run_dir / "data" / rec.name).exists()
    assert (run_dir / "data" / lig.name).exists()
    assert cfg_path.exists()


def test_oneshot_cfg_contains_expected_entries(tmp_path, monkeypatch):
    """
    author: ikekeku
    reviewer: fmclary
    category: one-shot test

    Verifies that a known input case produces a cfg file with expected entries.
    """
    inp = tmp_path / "input"
    inp.mkdir()
    rec = inp / "R1_rec.pdb"
    lig = inp / "L1_lig.pdb"
    rec.write_text("ATOM\n")
    lig.write_text("ATOM\n")

    work_root = tmp_path / "output"
    work_root.mkdir()

    monkeypatch.setattr(pdb_to_haddock, "run_command", lambda *a, **k: None)

    run_dir, cfg_path, _ = pdb_to_haddock.haddock_pipeline(str(rec), str(lig), runname="oneshot", mode="local", ncores=3, run_haddock=False, base_root=str(work_root), method="haddock_runs")

    text = cfg_path.read_text()
    # Assert runname and molecules appear in the cfg (one-shot known output)
    assert "run_dir = \"oneshot\"" in text
    assert f"\"data/{rec.name}\"" in text
    assert f"\"data/{lig.name}\"" in text
    assert "ncores = 3" in text


def test_edge_invalid_input_raises(tmp_path, monkeypatch):
    """
    author: ikekeku
    reviewer: fmclary
    category: edge test

    Verifies that providing a non-existent input file raises FileNotFoundError.
    """
    # Point to a non-existent receptor file -> copying should raise FileNotFoundError
    inp = tmp_path / "input"
    inp.mkdir()
    rec = inp / "missing_rec.pdb"  # does not exist
    lig = inp / "lig.pdb"
    lig.write_text("ATOM\n")

    work_root = tmp_path / "output"
    work_root.mkdir()

    monkeypatch.setattr(pdb_to_haddock, "run_command", lambda *a, **k: None)

    with pytest.raises(FileNotFoundError) as excinfo:
        pdb_to_haddock.haddock_pipeline(str(rec), str(lig), runname="edge", run_haddock=False, base_root=str(work_root), method="haddock_runs")

    # Ensure the error message includes the missing receptor filename for clarity
    assert "missing_rec.pdb" in str(excinfo.value)


def test_pattern_ncores_reflected_in_cfg(tmp_path, monkeypatch):
    """
    author: ikekeku
    reviewer: fmclary
    category: pattern test

    Verifies that the 'ncores' parameter is propagated into the generated
    cfg across multiple values.
    """
    inp = tmp_path / "input"
    inp.mkdir()
    rec = inp / "rec.pdb"
    lig = inp / "lig.pdb"
    rec.write_text("ATOM\n")
    lig.write_text("ATOM\n")

    work_root = tmp_path / "output"
    work_root.mkdir()

    monkeypatch.setattr(pdb_to_haddock, "run_command", lambda *a, **k: None)

    for cores in (1, 2, 4):
        run_dir, cfg_path, _ = pdb_to_haddock.haddock_pipeline(str(rec), str(lig), runname=f"p{cores}", ncores=cores, run_haddock=False, base_root=str(work_root), method="haddock_runs")
        text = cfg_path.read_text()
        assert f"ncores = {cores}" in text


def test_cfg_custom_runname_and_filename(tmp_path, monkeypatch):
    """
    author: ikekeku
    reviewer: fmclary
    category: functional test
    
    Verifies that a custom runname produces a cfg with the matching filename and
    run_dir entry.
    """
    inp = tmp_path / "input"
    inp.mkdir()
    rec = inp / "recX.pdb"
    lig = inp / "ligY.pdb"
    rec.write_text("ATOM\n")
    lig.write_text("ATOM\n")

    work_root = tmp_path / "output"
    work_root.mkdir()

    monkeypatch.setattr(pdb_to_haddock, "run_command", lambda *a, **k: None)

    run_dir, cfg_path, _ = pdb_to_haddock.haddock_pipeline(str(rec), str(lig), runname="custom_run", run_haddock=False, base_root=str(work_root), method="haddock_runs")

    assert cfg_path.exists()
    assert cfg_path.name == "custom_run.cfg"
    assert 'run_dir = "custom_run"' in cfg_path.read_text()


def test_cfg_includes_ambig_file(tmp_path, monkeypatch):
    """
    author: ikekeku
    reviewer: fmclary
    category: integration test
    
    Ensures that when an ambig restraints file is provided it is
    copied into data/ and referenced in the cfg.
    """
    inp = tmp_path / "input"
    inp.mkdir()
    rec = inp / "rec.pdb"
    lig = inp / "lig.pdb"
    ambig = inp / "restraints.tbl"
    rec.write_text("ATOM\n")
    lig.write_text("ATOM\n")
    ambig.write_text("{0} 1 2\n")

    work_root = tmp_path / "output"
    work_root.mkdir()

    monkeypatch.setattr(pdb_to_haddock, "run_command", lambda *a, **k: None)

    run_dir, cfg_path, _ = pdb_to_haddock.haddock_pipeline(str(rec), str(lig), runname="withambig", ambig=str(ambig), run_haddock=False, base_root=str(work_root), method="haddock_runs")

    # ambig file should be copied into data/
    assert (run_dir / "data" / ambig.name).exists()
    text = cfg_path.read_text()
    assert f'"data/{rec.name}"' in text
    assert f'"data/{lig.name}"' in text
    assert "ambig_fname = " in text


def test_cfg_out_root_and_method_respected(tmp_path, monkeypatch):
    """
    author: ikekeku
    reviewer: fmclary
    category: regression test
    
    Verifies that base_root and method params control where outputs are staged.
    """
    inp = tmp_path / "input"
    inp.mkdir()
    rec = inp / "rec.pdb"
    lig = inp / "lig.pdb"
    rec.write_text("ATOM\n")
    lig.write_text("ATOM\n")

    work_root = tmp_path / "myroot"
    work_root.mkdir()

    monkeypatch.setattr(pdb_to_haddock, "run_command", lambda *a, **k: None)

    run_dir, cfg_path, _ = pdb_to_haddock.haddock_pipeline(str(rec), str(lig), runname="o1", run_haddock=False, base_root=str(work_root), method="custom_method")

    # run_dir should be inside work_root/custom_method
    assert str(work_root.resolve()) in str(run_dir.resolve())
    assert "custom_method" in str(run_dir)
    assert cfg_path.exists()