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


def test_pipeline_constructs_expected_commands(sample_input_dirs, monkeypatch):
    """Ensure setup and simulation commands are constructed correctly."""
    rec = sample_input_dirs["rec"]
    lig = sample_input_dirs["lig"]
    work_root = sample_input_dirs["work_root"]

    # Record run_command calls
    recorder = CmdRecorder()
    monkeypatch.setattr(pdb_to_lightdock, "run_command", recorder)

    # Call pipeline with skip_postprocess to only get setup + simulation
    run_dir = pdb_to_lightdock.make_output_dir(
        str(rec), str(lig), base_root=str(work_root), method="lightdock_runs"
    )

    pdb_to_lightdock.lightdock_pipeline(
        str(rec),
        str(lig),
        working_dir=run_dir,
        swarms=50,
        glowworms=100,
        steps=10,
        swarm_list=[0],
        cores=2,
        skip_postprocess=True,
    )

    # There should be two calls recorded: setup and simulation (post-processing skipped)
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
    assert "-g" in setup_cmd and "100" in setup_cmd
    # ANM should be enabled by default
    assert "-anm" in setup_cmd

    # simulation command should reflect steps and cores and use -l for swarm list
    assert sim_cmd[0].endswith("lightdock3.py")
    assert sim_cmd[1] == "setup.json"
    assert "10" in sim_cmd
    assert "-c" in sim_cmd and "2" in sim_cmd
    assert "-l" in sim_cmd and "0" in sim_cmd


def test_full_pipeline_all_stages(sample_input_dirs, monkeypatch, tmp_path):
    """When skip_postprocess is False (default), all 5 stages are invoked.

    We create fake swarm directories so generation/clustering/ranking proceed.
    """
    rec = sample_input_dirs["rec"]
    lig = sample_input_dirs["lig"]
    work_root = sample_input_dirs["work_root"]

    recorder = CmdRecorder()
    monkeypatch.setattr(pdb_to_lightdock, "run_command", recorder)

    run_dir = pdb_to_lightdock.make_output_dir(
        str(rec), str(lig), base_root=str(work_root), method="lightdock_runs"
    )

    # Create fake swarm dirs with gso files so generation + clustering run
    for i in range(3):
        swarm = os.path.join(run_dir, f"swarm_{i}")
        os.makedirs(swarm, exist_ok=True)
        gso = os.path.join(swarm, "gso_5.out")
        with open(gso, "w") as f:
            f.write("# fake\n")

    pdb_to_lightdock.lightdock_pipeline(
        str(rec),
        str(lig),
        working_dir=run_dir,
        swarms=None,
        glowworms=None,
        steps=5,
        swarm_list=None,
        cores=1,
        # skip_postprocess defaults to False → all 5 stages
    )

    # Expected calls:
    # 1. setup, 2. simulation, 3-5. generation (3 swarms),
    # 6-8. clustering (3 swarms), 9. ranking
    assert len(recorder.calls) == 9

    # Verify generation commands iterate all swarms
    gen_cmds = [c for c, _ in recorder.calls if "lgd_generate_conformations.py" in c[0]]
    assert len(gen_cmds) == 3
    gen_swarms = sorted([c[3] for c in gen_cmds])  # out_file arg
    assert "swarm_0/gso_5.out" in gen_swarms[0]
    assert "swarm_1/gso_5.out" in gen_swarms[1]
    assert "swarm_2/gso_5.out" in gen_swarms[2]

    # Verify clustering commands
    cluster_cmds = [c for c, _ in recorder.calls if "lgd_cluster_bsas.py" in c[0]]
    assert len(cluster_cmds) == 3

    # Verify ranking command
    rank_cmds = [c for c, _ in recorder.calls if "lgd_rank.py" in c[0]]
    assert len(rank_cmds) == 1


def test_no_anm_flag(sample_input_dirs, monkeypatch):
    """When anm=False, the -anm flag should NOT appear in setup."""
    rec = sample_input_dirs["rec"]
    lig = sample_input_dirs["lig"]
    work_root = sample_input_dirs["work_root"]

    recorder = CmdRecorder()
    monkeypatch.setattr(pdb_to_lightdock, "run_command", recorder)

    run_dir = pdb_to_lightdock.make_output_dir(
        str(rec), str(lig), base_root=str(work_root), method="lightdock_runs"
    )

    pdb_to_lightdock.lightdock_pipeline(
        str(rec),
        str(lig),
        working_dir=run_dir,
        anm=False,
        skip_postprocess=True,
    )

    setup_cmd, _ = recorder.calls[0]
    assert "-anm" not in setup_cmd


def test_scoring_function_flag(sample_input_dirs, monkeypatch):
    """Custom scoring function is passed to the simulation step."""
    rec = sample_input_dirs["rec"]
    lig = sample_input_dirs["lig"]
    work_root = sample_input_dirs["work_root"]

    recorder = CmdRecorder()
    monkeypatch.setattr(pdb_to_lightdock, "run_command", recorder)

    run_dir = pdb_to_lightdock.make_output_dir(
        str(rec), str(lig), base_root=str(work_root), method="lightdock_runs"
    )

    pdb_to_lightdock.lightdock_pipeline(
        str(rec),
        str(lig),
        working_dir=run_dir,
        scoring="cpydock",
        skip_postprocess=True,
    )

    # Scoring function goes to simulation (lightdock3.py), not setup
    setup_cmd, _ = recorder.calls[0]
    sim_cmd, _ = recorder.calls[1]
    # Setup should NOT have -s cpydock (setup -s is for swarms)
    assert "cpydock" not in setup_cmd
    # Simulation should have -s cpydock
    assert "-s" in sim_cmd
    idx = sim_cmd.index("-s")
    assert sim_cmd[idx + 1] == "cpydock"


def test_backward_compat_generate_models_false(sample_input_dirs, monkeypatch):
    """Legacy generate_models=False maps to skip_postprocess=True."""
    rec = sample_input_dirs["rec"]
    lig = sample_input_dirs["lig"]
    work_root = sample_input_dirs["work_root"]

    recorder = CmdRecorder()
    monkeypatch.setattr(pdb_to_lightdock, "run_command", recorder)

    run_dir = pdb_to_lightdock.make_output_dir(
        str(rec), str(lig), base_root=str(work_root), method="lightdock_runs"
    )

    pdb_to_lightdock.lightdock_pipeline(
        str(rec),
        str(lig),
        working_dir=run_dir,
        generate_models=False,
    )

    # Only setup + simulation (no generation/clustering/ranking)
    assert len(recorder.calls) == 2
