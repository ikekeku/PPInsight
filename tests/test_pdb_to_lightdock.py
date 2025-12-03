"""
NOT AFFILIATED WITH HW3
Tests for `pdb_to_lightdock.py` command construction.
This test monkeypatches `run_command` so LightDock binaries are not invoked.
"""
import os

from ppinsight import pdb_to_lightdock


class CmdRecorder:
    """Record calls to run_command(cmd, cwd=...)."""

    def __init__(self):
        self.calls = []

    def __call__(self, cmd, cwd=None):
        # store a copy for inspection
        self.calls.append((list(cmd), cwd))


def test_pipeline_constructs_expected_commands(tmp_path, monkeypatch):
    # Prepare fake receptor/ligand files in a temp input dir
    inp = tmp_path / "input"
    inp.mkdir()
    """
    Tests for `pdb_to_lightdock.py` command construction.
    This test monkeypatches `run_command` so LightDock binaries are not invoked.
    """
    import os

    from ppinsight import pdb_to_lightdock


    class CmdRecorder:
        """Record calls to run_command(cmd, cwd=...)."""

        def __init__(self):
            self.calls = []

        def __call__(self, cmd, cwd=None):
            # store a copy for inspection
            self.calls.append((list(cmd), cwd))


    def test_pipeline_constructs_expected_commands(tmp_path, monkeypatch):
        # Prepare fake receptor/ligand files in a temp input dir
        inp = tmp_path / "input"
        inp.mkdir()
        rec = inp / "rec.pdb"
        lig = inp / "lig.pdb"
        rec.write_text("ATOM\n")
        lig.write_text("ATOM\n")

        # Prepare working dir
        work_root = tmp_path / "output"
        work_root.mkdir()

        # Record run_command calls
        recorder = CmdRecorder()
        monkeypatch.setattr(pdb_to_lightdock, "run_command", recorder)

        # Call pipeline with some options
        run_dir = pdb_to_lightdock.make_output_dir(str(rec), str(lig), base_root=str(work_root), method="lightdock_runs")

        pdb_to_lightdock.lightdock_pipeline(
            str(rec),
            str(lig),
            working_dir=run_dir,
            swarms=50,
            glowworms=100,
            steps=10,
            swarm_list=[0],
            cores=2,
            generate_models=False,
        )

        # There should be two calls recorded: setup and simulation (generation skipped)
        assert len(recorder.calls) == 2

        setup_cmd, setup_cwd = recorder.calls[0]
        sim_cmd, sim_cwd = recorder.calls[1]

        # cwd should be the run directory
        assert setup_cwd == run_dir
        assert sim_cwd == run_dir

        # setup command should start with lightdock3_setup.py and have filenames (basenames)
        assert setup_cmd[0].endswith("lightdock3_setup.py")
        assert setup_cmd[1] == os.path.basename(str(rec))
        assert setup_cmd[2] == os.path.basename(str(lig))
        assert "-s" in setup_cmd and "50" in setup_cmd
        assert "-g" in setup_cmd and "100" in setup_cmd

        # simulation command should reflect steps and cores and use -l for swarm list
        assert sim_cmd[0].endswith("lightdock3.py")
        assert sim_cmd[1] == "setup.json"
        assert "10" in sim_cmd
        assert "-c" in sim_cmd and "2" in sim_cmd
        assert "-l" in sim_cmd and "0" in sim_cmd


    def test_generation_runs_when_requested(tmp_path, monkeypatch):
        inp = tmp_path / "input"
        inp.mkdir()
        rec = inp / "rec.pdb"
        lig = inp / "lig.pdb"
        rec.write_text("ATOM\n")
        lig.write_text("ATOM\n")

        work_root = tmp_path / "output"
        work_root.mkdir()

        recorder = CmdRecorder()
        monkeypatch.setattr(pdb_to_lightdock, "run_command", recorder)

        run_dir = pdb_to_lightdock.make_output_dir(str(rec), str(lig), base_root=str(work_root), method="lightdock_runs")

        pdb_to_lightdock.lightdock_pipeline(
            str(rec),
            str(lig),
            working_dir=run_dir,
            swarms=None,
            glowworms=None,
            steps=5,
            swarm_list=None,
            cores=1,
            generate_models=True,
        )

        # now we expect three calls: setup, simulation, generation
        assert len(recorder.calls) == 3
        gen_cmd, gen_cwd = recorder.calls[2]
        assert gen_cmd[0].endswith("lgd_generate_conformations.py")
        # receptor/ligand passed to generator should be basenames
        assert gen_cmd[1] == os.path.basename(str(rec))
        assert gen_cmd[2] == os.path.basename(str(lig))
        # output .out path should be relative (swarm_0/gso_5.out)
        assert gen_cmd[3].endswith(os.path.join("swarm_0", "gso_5.out"))