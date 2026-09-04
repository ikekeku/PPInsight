"""Tests for ppinsight.consrank."""

import gzip
import os
import stat

import pytest

from ppinsight import consrank

# ---------------------------------------------------------------------------
# CONTROL file format
# ---------------------------------------------------------------------------

def test_write_control_matches_upstream_format(tmp_path):
    control_path = tmp_path / "CONTROL"
    consrank.write_control(
        str(control_path), "C", "BA", ["a.pdb", "b.pdb"],
        distance=5.0, gen_mat=0,
    )
    text = control_path.read_text(encoding="utf-8")
    lines = text.splitlines()

    assert lines[0].startswith("PairwiseChain1")
    assert "1" in lines[0]
    assert lines[1] == "C"
    assert lines[2].startswith("PairwiseChain2")
    assert "2" in lines[2]
    assert lines[3] == "B"
    assert lines[4] == "A"
    assert lines[5].startswith("CUTOffDistance")
    assert "5.0" in lines[5]
    assert lines[6].startswith("GenMat")
    assert lines[7].startswith("NumberOfPDBFiles")
    assert "2" in lines[7]
    assert lines[8:] == ["a.pdb", "b.pdb"]


def test_control_header_lines_stop_before_pdb_list():
    header = consrank._control_header_lines("C", "BA", distance=5.0, gen_mat=1)
    assert header[-1] == "GenMat\t\t\t1"
    assert not any(line.startswith("NumberOfPDBFiles") for line in header)


# ---------------------------------------------------------------------------
# Score parsing + selection (cut.f port)
# ---------------------------------------------------------------------------

def test_parse_consrank_scores(tmp_path):
    p = tmp_path / "Consrank_score.txt"
    p.write_text("model_1.pdb 0.42\nmodel_2.pdb 0.91\n\nmodel_3.pdb 0.10\n")
    rows = consrank._parse_consrank_scores(str(p))
    assert rows == [("model_1.pdb", 0.42), ("model_2.pdb", 0.91), ("model_3.pdb", 0.10)]


def test_parse_consrank_scores_skips_unparseable_lines(tmp_path):
    p = tmp_path / "Consrank_score.txt"
    p.write_text("model_1.pdb 0.5\nnot_a_score_line\nmodel_2.pdb notafloat\n")
    rows = consrank._parse_consrank_scores(str(p))
    assert rows == [("model_1.pdb", 0.5)]


def test_parse_consrank_scores_raises_when_empty(tmp_path):
    p = tmp_path / "Consrank_score.txt"
    p.write_text("\n\n")
    with pytest.raises(consrank.ConsrankError):
        consrank._parse_consrank_scores(str(p))


def test_select_top_fraction_keeps_highest_scoring_tail():
    scored = [("a", 0.1), ("b", 0.9), ("c", 0.5), ("d", 0.7)]
    kept = consrank._select_top_fraction(scored, cutoff=0.5)
    kept_names = {name for name, _ in kept}
    assert kept_names == {"b", "d"}


def test_select_top_fraction_never_drops_below_one():
    scored = [("a", 0.1), ("b", 0.2)]
    kept = consrank._select_top_fraction(scored, cutoff=0.01)
    assert len(kept) == 1
    assert kept[0][0] == "b"


def test_select_top_fraction_cutoff_one_keeps_all():
    scored = [("a", 0.1), ("b", 0.2), ("c", 0.3)]
    kept = consrank._select_top_fraction(scored, cutoff=1.0)
    assert len(kept) == 3


# ---------------------------------------------------------------------------
# Binary resolution
# ---------------------------------------------------------------------------

def test_consrank_binary_available_false_when_missing(tmp_path):
    assert consrank.consrank_binary_available(str(tmp_path / "nope")) is False


def test_consrank_binary_available_requires_executable_bit(tmp_path):
    fake = tmp_path / "CONSRANK"
    fake.write_text("#!/bin/sh\n")
    assert consrank.consrank_binary_available(str(fake)) is False
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    assert consrank.consrank_binary_available(str(fake)) is True


def test_require_consrank_binary_raises_with_helpful_message(tmp_path):
    with pytest.raises(consrank.ConsrankError, match="setup.sh"):
        consrank._require_consrank_binary(str(tmp_path / "missing"))


# ---------------------------------------------------------------------------
# Iteration loop
# ---------------------------------------------------------------------------

def test_run_iterative_consrank_shrinks_pool_and_stops(tmp_path, monkeypatch):
    pose_dir = tmp_path / "poses"
    pose_dir.mkdir()
    filenames = [f"m{i}.pdb" for i in range(10)]
    for name in filenames:
        (pose_dir / name).write_text("ATOM\n")

    monkeypatch.setattr(consrank, "consrank_binary_available", lambda path=None: True)

    call_count = {"n": 0}

    def _fake_run_command(cmd, cwd=None):
        call_count["n"] += 1
        # Read the CONTROL file's pool to know which poses are "in play".
        control_text = (pose_dir / "CONTROL").read_text()
        pool = [
            line for line in control_text.splitlines()
            if line.endswith(".pdb")
        ]
        # Fake deterministic ascending scores by index in the pool.
        lines = [f"{name} {idx / 10:.2f}" for idx, name in enumerate(pool)]
        (pose_dir / "Consrank_score.txt").write_text("\n".join(lines) + "\n")
        return ""

    monkeypatch.setattr(consrank, "run_command", _fake_run_command)

    result = consrank.run_iterative_consrank(
        str(pose_dir), filenames, "A", "B",
        cutoff=0.5, iterations=10, binary="/fake/CONSRANK",
    )

    history = result["history"]
    assert history[0]["n_input"] == 10
    assert history[0]["n_kept"] == 5
    # Pool keeps halving until it stops shrinking or bottoms out at a
    # single pose, at which point the loop stops early rather than
    # looping to the full --iterations count.
    last = history[-1]
    assert last["n_kept"] == last["n_input"] or last["n_kept"] <= 1
    assert len(history) < 10
    assert call_count["n"] == len(history)
    assert os.path.isfile(pose_dir / f"Consrank_score-{len(history)}.txt")
    assert os.path.isfile(pose_dir / f"CONTROL-{len(history)}.txt")


def test_run_iterative_consrank_requires_at_least_two_poses(tmp_path, monkeypatch):
    monkeypatch.setattr(consrank, "consrank_binary_available", lambda path=None: True)
    with pytest.raises(consrank.ConsrankError, match="at least 2 poses"):
        consrank.run_iterative_consrank(
            str(tmp_path), ["only_one.pdb"], "A", "B", binary="/fake/CONSRANK",
        )


# ---------------------------------------------------------------------------
# Pose discovery
# ---------------------------------------------------------------------------

def test_discover_poses_lightdock(tmp_path):
    run_dir = tmp_path / "run"
    for swarm in ("swarm_0", "swarm_1"):
        d = run_dir / swarm
        d.mkdir(parents=True)
        (d / "lightdock_1.pdb").write_text("ATOM\n")
    found = consrank.discover_poses(str(run_dir), "lightdock")
    assert len(found) == 2


def test_discover_poses_rosetta(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "decoy_1.pdb").write_text("ATOM\n")
    (run_dir / "decoy_2.pdb").write_text("ATOM\n")
    (run_dir / "not_a_decoy.pdb").write_text("ATOM\n")
    found = consrank.discover_poses(str(run_dir), "rosetta")
    assert len(found) == 2


def test_discover_poses_haddock(tmp_path):
    run_dir = tmp_path / "run" / "7_seletopclusts"
    run_dir.mkdir(parents=True)
    with gzip.open(run_dir / "cluster_1_model_1.pdb.gz", "wb") as fh:
        fh.write(b"ATOM\n")
    found = consrank.discover_poses(str(tmp_path / "run"), "haddock")
    assert len(found) == 1


def test_discover_poses_raises_when_none_found(tmp_path):
    with pytest.raises(consrank.ConsrankError, match="No rosetta poses"):
        consrank.discover_poses(str(tmp_path), "rosetta")


def test_discover_poses_respects_max_poses(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    for i in range(5):
        (run_dir / f"decoy_{i}.pdb").write_text("ATOM\n")
    found = consrank.discover_poses(str(run_dir), "rosetta", max_poses=2)
    assert len(found) == 2


# ---------------------------------------------------------------------------
# Staging poses (copy + gunzip)
# ---------------------------------------------------------------------------

def test_stage_poses_gunzips_and_copies(tmp_path):
    src_dir = tmp_path / "src"
    src_dir.mkdir()
    plain = src_dir / "plain.pdb"
    plain.write_text("ATOM plain\n")
    gz_path = src_dir / "cluster_1.pdb.gz"
    with gzip.open(gz_path, "wb") as fh:
        fh.write(b"ATOM gz\n")

    pose_dir = tmp_path / "staged"
    staged = consrank.stage_poses([str(plain), str(gz_path)], str(pose_dir))

    assert set(staged) == {"plain.pdb", "cluster_1.pdb"}
    assert (pose_dir / "cluster_1.pdb").read_text() == "ATOM gz\n"
    assert (pose_dir / "plain.pdb").read_text() == "ATOM plain\n"


def test_stage_poses_disambiguates_name_collisions(tmp_path):
    swarm0 = tmp_path / "swarm_0"
    swarm1 = tmp_path / "swarm_1"
    swarm0.mkdir()
    swarm1.mkdir()
    (swarm0 / "lightdock_1.pdb").write_text("ATOM 0\n")
    (swarm1 / "lightdock_1.pdb").write_text("ATOM 1\n")

    pose_dir = tmp_path / "staged"
    staged = consrank.stage_poses(
        [str(swarm0 / "lightdock_1.pdb"), str(swarm1 / "lightdock_1.pdb")],
        str(pose_dir),
    )
    assert len(staged) == 2
    assert len(set(staged)) == 2


def test_stage_poses_applies_engine_prefix(tmp_path):
    src_dir = tmp_path / "src"
    src_dir.mkdir()
    (src_dir / "decoy_1.pdb").write_text("ATOM\n")

    pose_dir = tmp_path / "staged"
    staged = consrank.stage_poses(
        [str(src_dir / "decoy_1.pdb")], str(pose_dir), prefix="rosetta",
    )
    assert staged == ["rosetta_decoy_1.pdb"]
    assert (pose_dir / "rosetta_decoy_1.pdb").is_file()


# ---------------------------------------------------------------------------
# Chain-ID collision relabeling (LightDock-specific)
# ---------------------------------------------------------------------------

def _atom_line(chain, serial=1):
    return (
        f"ATOM  {serial:>5}  CA  ALA {chain}{serial:>4}    "
        "0.000   0.000   0.000  1.00 20.00           C\n"
    )


def test_relabel_pose_chains_by_block_order(tmp_path):
    # Mirrors the real collision found in P00648_vs_P11540: receptor's
    # chains A,B then ligand's chains A,B, reusing letters across the
    # rec/lig boundary but staying contiguous per block.
    pose = tmp_path / "pose.pdb"
    lines = (
        [_atom_line("A", 1), _atom_line("A", 2)]  # receptor chain A (2 atoms)
        + [_atom_line("B", 3)]                     # receptor chain B (1 atom)
        + [_atom_line("A", 4)]                     # ligand chain A (1 atom)
        + [_atom_line("B", 5), _atom_line("B", 6)]  # ligand chain B (2 atoms)
    )
    pose.write_text("".join(lines))

    consrank._relabel_pose_chains_by_block_order(str(pose), ["A", "B", "C", "D"])

    from collections import Counter
    counts = Counter()
    with open(pose) as fh:
        for line in fh:
            counts[line[21]] += 1
    assert counts == {"A": 2, "B": 1, "C": 1, "D": 2}


def test_relabel_pose_chains_raises_on_unexpected_extra_block(tmp_path):
    pose = tmp_path / "pose.pdb"
    # 3 contiguous blocks but caller only expects 2 total (rec=1, lig=1).
    pose.write_text(_atom_line("A") + _atom_line("B") + _atom_line("C"))
    with pytest.raises(consrank.ConsrankError, match="more contiguous chain blocks"):
        consrank._relabel_pose_chains_by_block_order(str(pose), ["A", "B"])


def test_relabel_colliding_lightdock_chains_noop_when_disjoint(tmp_path):
    pose_dir = tmp_path / "poses"
    pose_dir.mkdir()
    (pose_dir / "m1.pdb").write_text(_atom_line("A") + _atom_line("B"))
    rec, lig = consrank.relabel_colliding_lightdock_chains(
        str(pose_dir), ["m1.pdb"], "A", "B",
    )
    assert (rec, lig) == ("A", "B")
    # File untouched.
    assert (pose_dir / "m1.pdb").read_text() == _atom_line("A") + _atom_line("B")


def test_relabel_colliding_lightdock_chains_end_to_end(tmp_path):
    pose_dir = tmp_path / "poses"
    pose_dir.mkdir()
    lines = (
        [_atom_line("A", 1), _atom_line("A", 2)]
        + [_atom_line("B", 3)]
        + [_atom_line("A", 4)]
        + [_atom_line("B", 5), _atom_line("B", 6)]
    )
    (pose_dir / "m1.pdb").write_text("".join(lines))

    rec, lig = consrank.relabel_colliding_lightdock_chains(
        str(pose_dir), ["m1.pdb"], "AB", "AB",
    )
    assert (rec, lig) == ("AB", "CD")
    from collections import Counter
    counts = Counter()
    with open(pose_dir / "m1.pdb") as fh:
        for line in fh:
            counts[line[21]] += 1
    assert counts == {"A": 2, "B": 1, "C": 1, "D": 2}


# ---------------------------------------------------------------------------
# Chain resolution
# ---------------------------------------------------------------------------

def _write_atom_lines(path, chains):
    lines = []
    for chain in chains:
        lines.append(
            f"ATOM      1  CA  ALA {chain}   1      11.104  13.207   2.010  1.00 20.00"
            "           C\n"
        )
    path.write_text("".join(lines))


def test_resolve_chains_haddock_is_fixed():
    rec, lig = consrank.resolve_chains("/irrelevant", "haddock", [])
    assert (rec, lig) == ("A", "B")


def test_resolve_chains_rosetta_reads_ppinsight_partners_comment(tmp_path):
    decoy = tmp_path / "decoy_1.pdb"
    decoy.write_text(
        _atom_line("A") + "TER\n" + "##Begin comments##\n"
        "ppinsight_partners ABC_DE\n" + "##End comments##\n"
    )
    rec, lig = consrank.resolve_chains("/irrelevant", "rosetta", [str(decoy)])
    assert (rec, lig) == ("ABC", "DE")


def test_resolve_chains_rosetta_raises_without_comment(tmp_path):
    decoy = tmp_path / "decoy_1.pdb"
    decoy.write_text(_atom_line("A"))
    with pytest.raises(consrank.ConsrankError, match="ppinsight_partners"):
        consrank.resolve_chains("/irrelevant", "rosetta", [str(decoy)])


def test_resolve_chains_lightdock_reads_original_inputs(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    rec = run_dir / "P00648.pdb"
    lig = run_dir / "P11540.pdb"
    _write_atom_lines(rec, ["A", "C"])
    _write_atom_lines(lig, ["B"])
    # A generated pose (prefixed lightdock_) must be ignored when locating
    # the two original top-level inputs.
    (run_dir / "lightdock_P00648.pdb").write_text("ATOM\n")

    rec_chains, lig_chains = consrank.resolve_chains(str(run_dir), "lightdock", [])
    assert rec_chains == "AC"
    assert lig_chains == "B"


def test_resolve_chains_lightdock_raises_without_exactly_two_inputs(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "only_one.pdb").write_text("ATOM\n")
    with pytest.raises(consrank.ConsrankError, match="Expected exactly 2"):
        consrank.resolve_chains(str(run_dir), "lightdock", [])


def test_resolve_chains_override_wins(tmp_path):
    rec, lig = consrank.resolve_chains(
        str(tmp_path), "lightdock", [], rec_chains="X", lig_chains="Y",
    )
    assert (rec, lig) == ("X", "Y")


def test_resolve_chains_requires_both_override_args(tmp_path):
    with pytest.raises(consrank.ConsrankError, match="must be given together"):
        consrank.resolve_chains(str(tmp_path), "lightdock", [], rec_chains="X")


# ---------------------------------------------------------------------------
# Output directory + scores DataFrame
# ---------------------------------------------------------------------------

def test_make_consrank_output_dir_auto_increments(tmp_path):
    base_root = str(tmp_path / "out")
    run1 = consrank.make_consrank_output_dir("pairA", base_root=base_root)
    run2 = consrank.make_consrank_output_dir("pairA", base_root=base_root)
    assert run1 != run2
    assert run2.endswith("_1")


def test_scores_to_dataframe_ranks_descending(tmp_path):
    (tmp_path / "m1.pdb").write_text("ATOM\n")
    (tmp_path / "m2.pdb").write_text("ATOM\n")
    df = consrank.scores_to_dataframe(
        [("m1.pdb", 0.3), ("m2.pdb", 0.9)],
        str(tmp_path),
        pair=("P1", "P2"),
    )
    assert list(df["pose_id"]) == ["m2", "m1"]
    assert list(df["pose_rank"]) == [1, 2]
    assert list(df["proteinA"]) == ["P1", "P1"]
    assert list(df["score_type"]) == ["consrank_score", "consrank_score"]


def test_scores_to_dataframe_uses_pose_engine_map(tmp_path):
    (tmp_path / "rosetta_m1.pdb").write_text("ATOM\n")
    (tmp_path / "haddock_m2.pdb").write_text("ATOM\n")
    df = consrank.scores_to_dataframe(
        [("rosetta_m1.pdb", 0.3), ("haddock_m2.pdb", 0.9)],
        str(tmp_path),
        pose_engine_map={"rosetta_m1.pdb": "rosetta", "haddock_m2.pdb": "haddock"},
    )
    assert dict(zip(df["pose_id"], df["model"], strict=True)) == {
        "haddock_m2": "haddock",
        "rosetta_m1": "rosetta",
    }


# ---------------------------------------------------------------------------
# Cross-engine pooling
# ---------------------------------------------------------------------------

def test_parse_pool_spec():
    assert consrank.parse_pool_spec("rosetta=/some/dir") == ("rosetta", "/some/dir")


def test_parse_pool_spec_rejects_bad_format():
    with pytest.raises(consrank.ConsrankError, match="ENGINE=RUN_DIR"):
        consrank.parse_pool_spec("no-equals-sign")


def test_parse_pool_spec_rejects_unknown_engine():
    with pytest.raises(consrank.ConsrankError, match="Unknown engine"):
        consrank.parse_pool_spec("swissdock=/some/dir")


def _make_source(engine, rec_chains, lig_chains, pose_dir, blocks):
    """Build a PoseSource with one staged pose matching *blocks*.

    *blocks* is a list of (chain_letter, n_atoms) contiguous runs, mirroring
    a real engine's pose layout (receptor's chains, then the ligand's).
    """
    filename = f"{engine}_pose.pdb"
    serial = 1
    lines = []
    for chain, n_atoms in blocks:
        for _ in range(n_atoms):
            lines.append(_atom_line(chain, serial))
            serial += 1
    (pose_dir / filename).write_text("".join(lines))
    source = consrank.PoseSource(
        engine=engine,
        run_dir=f"/fake/{engine}",
        pose_paths=[],
        rec_chains=rec_chains,
        lig_chains=lig_chains,
    )
    source.staged_filenames = [filename]
    return source


def test_harmonize_pool_chains_relabels_all_sources(tmp_path):
    # HADDOCK: rec=A (1 chain), lig=B (1 chain) -- already matches the
    # target scheme, but still gets rewritten for consistency.
    haddock_source = _make_source(
        "haddock", "A", "B", tmp_path, [("A", 2), ("B", 2)],
    )
    # LightDock: rec=C (1 chain), lig=D (1 chain) -- different original
    # letters than HADDOCK's, same chain *counts*, so poolable.
    lightdock_source = _make_source(
        "lightdock", "C", "D", tmp_path, [("C", 3), ("D", 1)],
    )

    rec_chains, lig_chains = consrank.harmonize_pool_chains(
        str(tmp_path), [haddock_source, lightdock_source],
    )

    assert (rec_chains, lig_chains) == ("A", "B")
    from collections import Counter
    haddock_counts = Counter()
    with open(tmp_path / "haddock_pose.pdb") as fh:
        for line in fh:
            haddock_counts[line[21]] += 1
    assert haddock_counts == {"A": 2, "B": 2}

    lightdock_counts = Counter()
    with open(tmp_path / "lightdock_pose.pdb") as fh:
        for line in fh:
            lightdock_counts[line[21]] += 1
    assert lightdock_counts == {"A": 3, "B": 1}


def test_harmonize_pool_chains_rejects_mismatched_chain_counts(tmp_path):
    # Rosetta: receptor kept 2 chains; HADDOCK: receptor normalized to 1.
    rosetta_source = _make_source(
        "rosetta", "AB", "C", tmp_path, [("A", 1), ("B", 1), ("C", 1)],
    )
    haddock_source = _make_source(
        "haddock", "A", "B", tmp_path, [("A", 1), ("B", 1)],
    )

    with pytest.raises(consrank.ConsrankError, match="mismatched chain topology"):
        consrank.harmonize_pool_chains(str(tmp_path), [rosetta_source, haddock_source])


# ---------------------------------------------------------------------------
# CLI: --pool end-to-end wiring
# ---------------------------------------------------------------------------

def test_cli_pool_ranks_across_engines(tmp_path, monkeypatch):
    # Rosetta pool: rec=A (1 chain), lig=B (1 chain), 2 decoys, matching
    # ppinsight_partners comment so the resolver can find the split.
    rosetta_dir = tmp_path / "rosetta_run"
    rosetta_dir.mkdir()
    for i in (1, 2):
        (rosetta_dir / f"decoy_{i}.pdb").write_text(
            _atom_line("A", 1) + _atom_line("B", 2) + "TER\n"
            "##Begin comments##\nppinsight_partners A_B\n##End comments##\n"
        )

    # HADDOCK pool: rec=A, lig=B (its fixed convention), 1 clustered model.
    haddock_dir = tmp_path / "haddock_run" / "7_seletopclusts"
    haddock_dir.mkdir(parents=True)
    with gzip.open(haddock_dir / "cluster_1_model_1.pdb.gz", "wb") as fh:
        fh.write((_atom_line("A", 1) + _atom_line("B", 2)).encode())

    monkeypatch.setattr(consrank, "consrank_binary_available", lambda path=None: True)
    monkeypatch.setattr(consrank, "_project_root", lambda: str(tmp_path))

    def _fake_run_command(cmd, cwd=None):
        pool = [
            name for name in os.listdir(cwd)
            if name.endswith(".pdb")
        ]
        lines = [f"{name} 0.5" for name in pool]
        with open(os.path.join(cwd, "Consrank_score.txt"), "w") as fh:
            fh.write("\n".join(lines) + "\n")
        return ""

    monkeypatch.setattr(consrank, "run_command", _fake_run_command)

    output_path = tmp_path / "combined.tsv"
    consrank.main([
        "--pool", f"rosetta={rosetta_dir}",
        "--pool", f"haddock={haddock_dir.parent}",
        "--iterations", "1",
        "-o", str(output_path),
    ])

    import pandas as pd
    df = pd.read_csv(output_path, sep="\t")
    assert set(df["model"]) == {"rosetta", "haddock"}
    assert len(df) == 3
