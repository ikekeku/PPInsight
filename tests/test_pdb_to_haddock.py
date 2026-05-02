"""
FORMATTED FOR HW3!
Tests for `pdb_to_haddock.py`.

The tests monkeypatch `run_command` so haddock3 is not actually invoked.
"""
import subprocess

import pytest

from ppinsight import pdb_to_haddock


def _make_atom_line(
    serial: int,
    atom: str,
    resname: str,
    chain: str,
    resseq: int,
    *,
    segid: str | None = None,
    element: str = "C",
) -> str:
    segid = chain if segid is None else segid
    return (
        f"ATOM  {serial:5d} {atom:>4} {resname:>3} {chain:1}{resseq:4d}    "
        f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}{1.00:6.2f}{20.00:6.2f}      "
        f"{segid:<4}{element:>2}\n"
    )


def _write_simple_pdb(path, *, chain: str, segid: str | None = None):
    path.write_text(
        _make_atom_line(1, "N", "GLY", chain, 1, segid=segid, element="N")
        + _make_atom_line(2, "CA", "GLY", chain, 1, segid=segid)
        + "END\n"
    )


def _atom_lines(path):
    return [
        line
        for line in path.read_text().splitlines()
        if line.startswith("ATOM")
    ]


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
    _write_simple_pdb(rec, chain="A")
    _write_simple_pdb(lig, chain="B")

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
    _write_simple_pdb(rec, chain="A")
    _write_simple_pdb(lig, chain="B")

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
    _write_simple_pdb(lig, chain="B")

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
    _write_simple_pdb(rec, chain="A")
    _write_simple_pdb(lig, chain="B")

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
    _write_simple_pdb(rec, chain="A")
    _write_simple_pdb(lig, chain="B")

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
    _write_simple_pdb(rec, chain="A")
    _write_simple_pdb(lig, chain="B")
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
    assert "cmrest = true" not in text


def test_cfg_carries_cmrest_into_flexref_for_ab_initio_runs(tmp_path, monkeypatch):
    inp = tmp_path / "input"
    inp.mkdir()
    rec = inp / "rec.pdb"
    lig = inp / "lig.pdb"
    _write_simple_pdb(rec, chain="A")
    _write_simple_pdb(lig, chain="B")

    work_root = tmp_path / "output"
    work_root.mkdir()

    monkeypatch.setattr(pdb_to_haddock, "run_command", lambda *a, **k: None)

    _, cfg_path, _ = pdb_to_haddock.haddock_pipeline(
        str(rec),
        str(lig),
        runname="abinitio_cmrest",
        run_haddock=False,
        base_root=str(work_root),
        method="haddock_runs",
    )

    text = cfg_path.read_text()
    assert text.count("cmrest = true") == 2


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
    _write_simple_pdb(rec, chain="A")
    _write_simple_pdb(lig, chain="B")

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


def test_copy_inputs_normalizes_multichain_partner_without_restraints(tmp_path):
    data_dir = tmp_path / "data"
    data_dir.mkdir()

    rec = tmp_path / "rec.pdb"
    lig = tmp_path / "lig.pdb"
    rec.write_text(
        _make_atom_line(1, "N", "GLY", "R", 10, segid="R", element="N")
        + _make_atom_line(2, "CA", "GLY", "R", 10, segid="R")
        + "END\n"
    )
    lig.write_text(
        _make_atom_line(1, "N", "ALA", "H", 5, segid="H", element="N")
        + _make_atom_line(2, "CA", "ALA", "H", 5, segid="H")
        + "TER       3      ALA H   5\n"
        + _make_atom_line(4, "N", "SER", "J", 5, segid="J", element="N")
        + _make_atom_line(5, "CA", "SER", "J", 5, segid="J")
        + "END\n"
    )

    rec_dst, lig_dst, ambig_dst = pdb_to_haddock.copy_inputs(
        data_dir,
        str(rec),
        str(lig),
    )

    assert ambig_dst is None

    rec_atoms = _atom_lines(rec_dst)
    lig_atoms = _atom_lines(lig_dst)

    assert {line[21] for line in rec_atoms} == {"A"}
    assert {line[21] for line in lig_atoms} == {"B"}
    assert {line[72:76].strip() for line in rec_atoms} == {"A"}
    assert {line[72:76].strip() for line in lig_atoms} == {"B"}

    ligand_residues = []
    last_residue = None
    for line in lig_atoms:
        residue_id = line[21:27]
        if residue_id != last_residue:
            ligand_residues.append(line[22:26].strip())
            last_residue = residue_id

    assert ligand_residues == ["1", "2"]


def test_copy_inputs_normalizes_shared_chain_ids_without_restraints(tmp_path):
    data_dir = tmp_path / "data"
    data_dir.mkdir()

    rec = tmp_path / "rec.pdb"
    lig = tmp_path / "lig.pdb"
    _write_simple_pdb(rec, chain="A")
    _write_simple_pdb(lig, chain="A")

    rec_dst, lig_dst, _ = pdb_to_haddock.copy_inputs(
        data_dir,
        str(rec),
        str(lig),
    )

    rec_atoms = _atom_lines(rec_dst)
    lig_atoms = _atom_lines(lig_dst)
    assert {line[21] for line in rec_atoms} == {"A"}
    assert {line[21] for line in lig_atoms} == {"B"}


def test_copy_inputs_prefers_dbref_chains_for_matching_accession(tmp_path):
    data_dir = tmp_path / "data"
    data_dir.mkdir()

    rec = tmp_path / "rec.pdb"
    lig = tmp_path / "P12345.pdb"
    _write_simple_pdb(rec, chain="A")
    lig.write_text(
        "DBREF  1ABC A    1     1  UNP    P12345   TEST_HUMAN      1      1\n"
        "DBREF  1ABC B    1     1  UNP    P12345   TEST_HUMAN      2      2\n"
        "DBREF  1ABC X    1     1  UNP    Q99999   OTHER_HUMAN     1      1\n"
        + _make_atom_line(1, "N", "GLY", "A", 1, segid="A", element="N")
        + _make_atom_line(2, "CA", "GLY", "A", 1, segid="A")
        + "TER       3      GLY A   1\n"
        + _make_atom_line(4, "N", "SER", "B", 5, segid="B", element="N")
        + _make_atom_line(5, "CA", "SER", "B", 5, segid="B")
        + "TER       6      SER B   5\n"
        + _make_atom_line(7, "N", "TYR", "X", 9, segid="X", element="N")
        + _make_atom_line(8, "CA", "TYR", "X", 9, segid="X")
        + "END\n"
    )

    _, lig_dst, _ = pdb_to_haddock.copy_inputs(
        data_dir,
        str(rec),
        str(lig),
    )

    lig_atoms = _atom_lines(lig_dst)
    assert len(lig_atoms) == 4
    assert {line[21] for line in lig_atoms} == {"B"}

    ligand_residues = []
    last_residue = None
    for line in lig_atoms:
        residue_id = line[21:27]
        if residue_id != last_residue:
            ligand_residues.append(line[22:26].strip())
            last_residue = residue_id

    assert ligand_residues == ["1", "2"]


def test_copy_inputs_rejects_multichain_inputs_with_ambig(tmp_path):
    data_dir = tmp_path / "data"
    data_dir.mkdir()

    rec = tmp_path / "rec.pdb"
    lig = tmp_path / "lig.pdb"
    ambig = tmp_path / "ambig.tbl"
    _write_simple_pdb(rec, chain="A")
    lig.write_text(
        _make_atom_line(1, "N", "ALA", "H", 5, segid="H", element="N")
        + _make_atom_line(2, "CA", "ALA", "H", 5, segid="H")
        + "TER       3      ALA H   5\n"
        + _make_atom_line(4, "N", "SER", "J", 5, segid="J", element="N")
        + _make_atom_line(5, "CA", "SER", "J", 5, segid="J")
        + "END\n"
    )
    ambig.write_text("assign (segid A) (segid B) 2.0 2.0 0.0\n")

    with pytest.raises(RuntimeError, match="before using --ambig"):
        pdb_to_haddock.copy_inputs(data_dir, str(rec), str(lig), str(ambig))


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


def test_run_in_container_builds_docker_command_with_haddock_entrypoint(
    tmp_path,
    monkeypatch,
):
    commands = []
    run_dir = tmp_path / "output" / "haddock_runs" / "run1"
    run_dir.mkdir(parents=True)

    monkeypatch.setattr(
        pdb_to_haddock,
        "run_command",
        lambda cmd, cwd=None: commands.append((list(cmd), cwd)),
    )

    pdb_to_haddock._run_in_container(
        "docker",
        "ghcr.io/haddocking/haddock3:latest",
        str(tmp_path),
        run_dir,
        "run1.cfg",
    )

    assert commands == [
        (
            [
                "docker",
                "run",
                "--rm",
                "-v",
                f"{tmp_path}:/workspace",
                "-w",
                "/workspace/output/haddock_runs/run1",
                "--entrypoint",
                "haddock3",
                "ghcr.io/haddocking/haddock3:latest",
                "run1.cfg",
            ],
            None,
        )
    ]


def test_run_in_container_builds_apptainer_command_in_run_dir(
    tmp_path,
    monkeypatch,
):
    commands = []
    run_dir = tmp_path / "output" / "haddock_runs" / "run1"
    run_dir.mkdir(parents=True)

    monkeypatch.setattr(
        pdb_to_haddock,
        "run_command",
        lambda cmd, cwd=None: commands.append((list(cmd), cwd)),
    )

    pdb_to_haddock._run_in_container(
        "apptainer",
        "ghcr.io/haddocking/haddock3:latest",
        str(tmp_path),
        run_dir,
        "run1.cfg",
    )

    assert commands == [
        (
            [
                "apptainer",
                "exec",
                "--bind",
                f"{tmp_path}:/workspace",
                "--pwd",
                "/workspace/output/haddock_runs/run1",
                "docker://ghcr.io/haddocking/haddock3:latest",
                "haddock3",
                "run1.cfg",
            ],
            None,
        )
    ]


def test_main_reports_runtime_error_cleanly(tmp_path, monkeypatch, capsys):
    inp = tmp_path / "input"
    inp.mkdir()
    rec = inp / "rec.pdb"
    lig = inp / "lig.pdb"
    _write_simple_pdb(rec, chain="A")
    _write_simple_pdb(lig, chain="B")

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


def test_main_reports_called_process_error_cleanly(tmp_path, monkeypatch, capsys):
    inp = tmp_path / "input"
    inp.mkdir()
    rec = inp / "rec.pdb"
    lig = inp / "lig.pdb"
    _write_simple_pdb(rec, chain="A")
    _write_simple_pdb(lig, chain="B")

    monkeypatch.setattr(
        pdb_to_haddock,
        "haddock_pipeline",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            subprocess.CalledProcessError(2, ["docker", "run", "... "])
        ),
    )

    with pytest.raises(SystemExit) as exc_info:
        pdb_to_haddock.main([str(rec), str(lig), "--run", "--container", "docker"])

    captured = capsys.readouterr()
    assert exc_info.value.code == 1
    assert "ERROR: command failed with exit code 2: docker run ... " in captured.err
