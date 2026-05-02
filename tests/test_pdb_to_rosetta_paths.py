"""Path-level tests for pdb_to_rosetta output directory behavior."""

import os
import sys
import types

from ppinsight import pdb_to_rosetta


def test_make_output_dir_auto_increments_when_existing(tmp_path):
    rec = tmp_path / "rec.pdb"
    lig = tmp_path / "lig.pdb"
    rec.write_text("ATOM\n", encoding="utf-8")
    lig.write_text("ATOM\n", encoding="utf-8")

    run1 = pdb_to_rosetta._make_output_dir(
        str(rec),
        str(lig),
        base_root=str(tmp_path / "out"),
        method="rosetta_runs",
    )
    run2 = pdb_to_rosetta._make_output_dir(
        str(rec),
        str(lig),
        base_root=str(tmp_path / "out"),
        method="rosetta_runs",
    )

    assert run1 != run2
    assert run2.endswith("_1")


def test_pyrosetta_installer_env_prefers_active_env(monkeypatch):
    monkeypatch.setenv("PATH", "/usr/bin")
    monkeypatch.setenv("PYTHONPATH", "/tmp/project-src")
    monkeypatch.setenv("PIP_USER", "1")
    monkeypatch.setenv("PYTHONNOUSERSITE", "0")

    env_bin_dir = os.path.dirname(sys.executable)

    with pdb_to_rosetta._pyrosetta_installer_env():
        assert os.environ["PATH"].split(os.pathsep)[0] == env_bin_dir
        assert os.environ["PYTHONNOUSERSITE"] == "1"
        assert "PYTHONPATH" not in os.environ
        assert "PIP_USER" not in os.environ

    assert os.environ["PATH"] == "/usr/bin"
    assert os.environ["PYTHONPATH"] == "/tmp/project-src"
    assert os.environ["PIP_USER"] == "1"
    assert os.environ["PYTHONNOUSERSITE"] == "0"


def test_cli_can_disable_dbref_auto_filter(monkeypatch):
    class _FakePipeline:
        kwargs = None

        def __init__(self, *args, **kwargs):
            type(self).kwargs = kwargs

        def run(self):
            return {"final_score": 0.0}

        def save_scores(self, output_path):
            return None

    monkeypatch.setattr(
        pdb_to_rosetta,
        "resolve_input_path",
        lambda path, search_root=None: f"/tmp/{path}.pdb",
    )
    monkeypatch.setattr(pdb_to_rosetta, "_ensure_pyrosetta", lambda: None)
    monkeypatch.setattr(pdb_to_rosetta, "_make_output_dir", lambda *a, **k: "/tmp/out")

    monkeypatch.setitem(
        sys.modules,
        "ppinsight.docking",
        types.SimpleNamespace(DockingPipeline=_FakePipeline),
    )

    try:
        pdb_to_rosetta.main(["rec", "lig", "--no-auto-filter", "--quiet"])
    finally:
        sys.modules.pop("ppinsight.docking", None)

    assert _FakePipeline.kwargs["auto_filter"] is False


def test_cli_creates_explicit_output_dir(monkeypatch, tmp_path):
    output_dir = tmp_path / "rosetta_output"

    class _FakePipeline:
        def __init__(self, *args, **kwargs):
            pass

        def run(self):
            return {"final_score": 0.0}

        def save_scores(self, output_path):
            assert os.path.isdir(output_dir)
            assert output_path == str(output_dir / "docking_scores.csv")
            with open(output_path, "w", encoding="utf-8") as handle:
                handle.write("run,total_score,i_sc\n1,-1.0,-2.0\n")

    monkeypatch.setattr(
        pdb_to_rosetta,
        "resolve_input_path",
        lambda path, search_root=None: f"/tmp/{path}.pdb",
    )
    monkeypatch.setattr(pdb_to_rosetta, "_ensure_pyrosetta", lambda: None)

    monkeypatch.setitem(
        sys.modules,
        "ppinsight.docking",
        types.SimpleNamespace(DockingPipeline=_FakePipeline),
    )

    try:
        pdb_to_rosetta.main([
            "rec",
            "lig",
            "--output-dir",
            str(output_dir),
            "--quiet",
            "--no-cluster",
        ])
    finally:
        sys.modules.pop("ppinsight.docking", None)

    assert output_dir.is_dir()
    assert (output_dir / "docking_scores.csv").is_file()


def test_rosetta_outputs_skip_clustering_without_decoy_metadata(tmp_path):
    csv_path = tmp_path / "docking_scores.csv"
    csv_path.write_text("run,total_score,i_sc\n1,-100.0,-10.0\n2,-90.0,-8.0\n")

    ready, reason = pdb_to_rosetta._rosetta_outputs_support_clustering(
        str(csv_path),
        str(tmp_path),
    )

    assert ready is False
    assert "per-decoy identifiers" in reason


def test_rosetta_outputs_allow_clustering_with_decoys(tmp_path):
    csv_path = tmp_path / "docking_scores.csv"
    csv_path.write_text(
        "run,description,total_score,i_sc\n"
        "1,decoy_1,-100.0,-10.0\n"
        "2,decoy_2,-90.0,-8.0\n"
    )
    (tmp_path / "decoy_1.pdb").write_text("ATOM\n", encoding="utf-8")
    (tmp_path / "decoy_2.pdb").write_text("ATOM\n", encoding="utf-8")

    ready, reason = pdb_to_rosetta._rosetta_outputs_support_clustering(
        str(csv_path),
        str(tmp_path),
    )

    assert ready is True
    assert reason is None
