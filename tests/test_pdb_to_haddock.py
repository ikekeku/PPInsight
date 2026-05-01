"""
FORMATTED FOR HW3!
Tests for `pdb_to_haddock.py`.

The tests monkeypatch `run_command` so haddock3 is not actually invoked.
"""
import pytest

from ppinsight import pdb_to_haddock


class CmdRecorder:
    """Record calls to run_command(cmd, cwd=...).

    A tiny helper used by tests to capture subprocess calls.
    """

    def __init__(self):
        self.calls = []

    def __call__(self, cmd, cwd=None):
        self.calls.append((list(cmd), cwd))

    def reset(self):
        """Clear recorded calls."""
        self.calls = []


def test_smoke_haddock_pipeline_runs(tmp_path, monkeypatch):
    """Smoke test for haddock_pipeline: stages are created and cfg written."""
    # Metadata: author=ikekeku, reviewer=fmclary, category=smoke test

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
    run_dir, cfg_path, _ = pdb_to_haddock.haddock_pipeline(
        str(rec),
        str(lig),
        runname="smoke",
        run_haddock=False,
        base_root=str(work_root),
        method="haddock_runs",
    )

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

    _, cfg_path, _ = pdb_to_haddock.haddock_pipeline(
        str(rec),
        str(lig),
        runname="oneshot",
        mode="local",
        ncores=3,
        run_haddock=False,
        base_root=str(work_root),
        method="haddock_runs",
    )

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
        pdb_to_haddock.haddock_pipeline(
            str(rec),
            str(lig),
            runname="edge",
            run_haddock=False,
            base_root=str(work_root),
            method="haddock_runs",
        )

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
        _, cfg_path, _ = pdb_to_haddock.haddock_pipeline(
            str(rec),
            str(lig),
            runname=f"p{cores}",
            ncores=cores,
            run_haddock=False,
            base_root=str(work_root),
            method="haddock_runs",
        )
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

    _, cfg_path, _ = pdb_to_haddock.haddock_pipeline(
        str(rec),
        str(lig),
        runname="custom_run",
        run_haddock=False,
        base_root=str(work_root),
        method="haddock_runs",
    )

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

    run_dir, cfg_path, _ = pdb_to_haddock.haddock_pipeline(
        str(rec),
        str(lig),
        runname="withambig",
        ambig=str(ambig),
        run_haddock=False,
        base_root=str(work_root),
        method="haddock_runs",
    )

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

    run_dir, cfg_path, _ = pdb_to_haddock.haddock_pipeline(
        str(rec),
        str(lig),
        runname="o1",
        run_haddock=False,
        base_root=str(work_root),
        method="custom_method",
    )

    # run_dir should be inside work_root/custom_method
    assert str(work_root.resolve()) in str(run_dir.resolve())
    assert "custom_method" in str(run_dir)
    assert cfg_path.exists()


def test_execute_haddock_run_uses_container_before_local_cns_check(
    tmp_path,
    monkeypatch,
):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    cfg_path = run_dir / "test.cfg"
    cfg_path.write_text("run_dir = \"test\"\n")

    container_calls = []

    monkeypatch.setattr(
        pdb_to_haddock,
        "_detect_container_runtime",
        lambda: "docker",
    )
    monkeypatch.setattr(
        pdb_to_haddock,
        "_docker_runtime_usable",
        lambda docker_executable=None: True,
    )
    monkeypatch.setattr(
        pdb_to_haddock,
        "_check_cns_compatibility",
        lambda: pytest.fail("local CNS check should be skipped in container mode"),
    )
    monkeypatch.setattr(
        pdb_to_haddock,
        "_run_in_container",
        lambda runtime, image, host_workspace, run_dir_path, cfg_name: (
            container_calls.append(
                (runtime, image, host_workspace, run_dir_path, cfg_name)
            )
        ),
    )

    executed = pdb_to_haddock._execute_haddock_run(
        run_dir,
        cfg_path,
        run_haddock=True,
        haddock_cmd="haddock3",
        container="auto",
        container_image="ghcr.io/haddocking/haddock3:latest",
        workspace_root=str(tmp_path),
    )

    assert container_calls == [
        (
            "docker",
            "ghcr.io/haddocking/haddock3:latest",
            str(tmp_path),
            run_dir,
            cfg_path.name,
        )
    ]
    assert executed.startswith("container:docker image=")


def test_execute_haddock_run_checks_cns_compatibility_for_local_runs(
    tmp_path,
    monkeypatch,
):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    cfg_path = run_dir / "test.cfg"
    cfg_path.write_text("run_dir = \"test\"\n")

    monkeypatch.setattr(pdb_to_haddock, "_detect_container_runtime", lambda: None)
    monkeypatch.setattr(
        pdb_to_haddock,
        "_check_cns_compatibility",
        lambda: (_ for _ in ()).throw(RuntimeError("arch mismatch")),
    )
    monkeypatch.setattr(
        pdb_to_haddock,
        "run_command",
        lambda *args, **kwargs: pytest.fail("local run should not start on mismatch"),
    )

    with pytest.raises(RuntimeError, match="arch mismatch"):
        pdb_to_haddock._execute_haddock_run(
            run_dir,
            cfg_path,
            run_haddock=True,
            haddock_cmd="haddock3",
            container="auto",
            container_image="ghcr.io/haddocking/haddock3:latest",
            workspace_root=str(tmp_path),
        )


def test_detect_container_runtime_skips_unusable_docker(monkeypatch):
    def fake_which(name):
        if name == "docker":
            return "/usr/local/bin/docker"
        if name == "apptainer":
            return "/usr/bin/apptainer"
        return None

    monkeypatch.setattr(pdb_to_haddock.shutil, "which", fake_which)
    monkeypatch.setattr(
        pdb_to_haddock,
        "_docker_runtime_usable",
        lambda docker_executable=None: False,
    )

    assert pdb_to_haddock._detect_container_runtime() == "apptainer"


def test_execute_haddock_run_explicit_docker_requires_daemon(
    tmp_path,
    monkeypatch,
):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    cfg_path = run_dir / "test.cfg"
    cfg_path.write_text("run_dir = \"test\"\n")

    monkeypatch.setattr(
        pdb_to_haddock,
        "_docker_runtime_usable",
        lambda docker_executable=None: False,
    )
    monkeypatch.setattr(
        pdb_to_haddock,
        "_run_in_container",
        lambda *args, **kwargs: pytest.fail(
            "container run should not start when Docker daemon is unavailable"
        ),
    )

    with pytest.raises(RuntimeError, match="Docker is installed but the Docker daemon"):
        pdb_to_haddock._execute_haddock_run(
            run_dir,
            cfg_path,
            run_haddock=True,
            haddock_cmd="haddock3",
            container="docker",
            container_image="ghcr.io/haddocking/haddock3:latest",
            workspace_root=str(tmp_path),
        )


def test_main_reports_runtime_error_cleanly(tmp_path, monkeypatch, capsys):
    inp = tmp_path / "input"
    inp.mkdir()
    rec = inp / "rec.pdb"
    lig = inp / "lig.pdb"
    rec.write_text("ATOM\n")
    lig.write_text("ATOM\n")

    monkeypatch.setattr(
        pdb_to_haddock,
        "haddock_pipeline",
        lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("arch mismatch")),
    )

    with pytest.raises(SystemExit) as exc_info:
        pdb_to_haddock.main([str(rec), str(lig), "--run"])

    captured = capsys.readouterr()
    assert exc_info.value.code == 1
    assert "ERROR: arch mismatch" in captured.err
