"""Tests for parse_pairs and batch_dock modules."""


from types import SimpleNamespace

import pandas as pd
import pytest

from ppinsight import parse_pairs

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_table(tmp_path):
    """Create a minimal RTK-interactome-style TSV."""
    tsv = tmp_path / "table.tsv"
    # Mimics the real format: leading tab, header rows, data rows
    tsv.write_text(
        "\tProteinA\t\tProtein B\t\t\t\t\t\t\t\t\n"
        "\tFamily\tMember\tInteractions\t\t\t\t\tNon-Interactions\t\t\t\n"
        "\tEphrin\tEPHA1\tERBB2 [6]\tFGFR2 [7]\t\t\t\tEGFR [7]\t\t\t\n"
        "\t\tEPHA2\tEGFR [5][8]\tERBB2 [5]\t\t\t\tMERTK [7]\t\t\t\n"
        "\tErbB\tEGFR\tFGFR1 [5]\t\t\t\t\t\t\t\t\n"
    )
    return str(tsv)


# ---------------------------------------------------------------------------
# parse_pairs
# ---------------------------------------------------------------------------

class TestParsePairs:
    def test_basic(self, sample_table):
        df = parse_pairs.parse_interaction_table(sample_table)
        assert "proteinA" in df.columns
        assert "proteinB" in df.columns
        assert "label" in df.columns
        assert "family" in df.columns
        assert len(df) > 0

    def test_interactions_and_non_interactions(self, sample_table):
        df = parse_pairs.parse_interaction_table(sample_table)
        labels = set(df["label"])
        assert "interaction" in labels
        assert "non-interaction" in labels

    def test_strips_references(self, sample_table):
        df = parse_pairs.parse_interaction_table(sample_table)
        # No brackets should remain in proteinB
        for val in df["proteinB"]:
            assert "[" not in val
            assert "]" not in val

    def test_references_extracted(self, sample_table):
        df = parse_pairs.parse_interaction_table(sample_table)
        # EGFR [5][8] should have references "5,8"
        egfr_row = df[(df["proteinA"] == "EPHA2") & (df["proteinB"] == "EGFR")]
        assert len(egfr_row) == 1
        assert "5" in egfr_row.iloc[0]["references"]
        assert "8" in egfr_row.iloc[0]["references"]

    def test_family_propagation(self, sample_table):
        df = parse_pairs.parse_interaction_table(sample_table)
        # EPHA2 is under Ephrin family (continuation row with no family)
        epha2_rows = df[df["proteinA"] == "EPHA2"]
        assert (epha2_rows["family"] == "Ephrin").all()

    def test_real_table(self):
        """Test against the real RTK interactome table if present."""
        import os
        real_path = os.path.join(
            os.path.dirname(__file__), "..",
            "RTK Interactome Tomasz Draft Table 1. Discovered cross-family "
            "interactions for PPInsight - Sheet1.tsv"
        )
        if not os.path.isfile(real_path):
            pytest.skip("Real RTK table not found")
        df = parse_pairs.parse_interaction_table(real_path)
        assert len(df) > 100  # should have 180+ pairs
        assert (df["label"] == "interaction").sum() > 50
        assert (df["label"] == "non-interaction").sum() > 20

    def test_empty_table(self, tmp_path):
        tsv = tmp_path / "empty.tsv"
        # Table with headers but no data rows at all
        tsv.write_text(
            "\tProteinA\t\tProtein B\t\t\t\t\t\t\t\t\n"
            "\tFamily\tMember\tInteractions\t\t\t\t\tNon-Interactions\t\t\t\n"
        )
        with pytest.raises(ValueError, match="No protein pairs"):
            parse_pairs.parse_interaction_table(str(tsv))


class TestParsePairsCLI:
    def test_basic(self, sample_table, tmp_path):
        out = str(tmp_path / "pairs.csv")
        parse_pairs.main([sample_table, "-o", out])
        df = pd.read_csv(out)
        assert len(df) > 0
        assert "label" in df.columns

    def test_tsv_output(self, sample_table, tmp_path):
        out = str(tmp_path / "pairs.tsv")
        parse_pairs.main([sample_table, "-o", out])
        df = pd.read_csv(out, sep="\t")
        assert len(df) > 0

    def test_stats_flag(self, sample_table, tmp_path, capsys):
        out = str(tmp_path / "pairs.csv")
        parse_pairs.main([sample_table, "-o", out, "--stats"])
        captured = capsys.readouterr()
        assert "Interactions:" in captured.out
        assert "Non-interactions:" in captured.out


# ---------------------------------------------------------------------------
# batch_dock (unit tests — no actual docking)
# ---------------------------------------------------------------------------

class TestBatchDock:
    def test_dry_run(self, tmp_path):
        from ppinsight import batch_dock

        # Create a minimal pairs DataFrame
        pairs_df = pd.DataFrame({
            "proteinA": ["EPHA1", "EPHA2"],
            "proteinB": ["ERBB2", "EGFR"],
            "label": ["interaction", "interaction"],
        })

        results = batch_dock.batch_dock(
            pairs_df,
            engines=["lightdock"],
            pdb_dir=str(tmp_path),  # no PDBs here
            dry_run=True,
        )
        assert len(results) == 2
        assert (results["status"] == "dry_run").all()

    def test_missing_pdbs(self, tmp_path):
        from ppinsight import batch_dock

        pairs_df = pd.DataFrame({
            "proteinA": ["FAKE_A"],
            "proteinB": ["FAKE_B"],
            "label": ["interaction"],
        })

        results = batch_dock.batch_dock(
            pairs_df,
            engines=["lightdock"],
            pdb_dir=str(tmp_path),
        )
        assert len(results) == 1
        assert results.iloc[0]["status"] == "pdb_missing_A"

    def test_limit(self, tmp_path):
        from ppinsight import batch_dock

        pairs_df = pd.DataFrame({
            "proteinA": ["A", "B", "C"],
            "proteinB": ["X", "Y", "Z"],
            "label": ["interaction"] * 3,
        })

        results = batch_dock.batch_dock(
            pairs_df,
            engines=["lightdock"],
            pdb_dir=str(tmp_path),
            limit=1,
            dry_run=True,
        )
        assert len(results) == 1

    def test_engine_kwargs_forwarded_to_runner(self, tmp_path, monkeypatch):
        from ppinsight import batch_dock

        (tmp_path / "A.pdb").write_text("END\n", encoding="utf-8")
        (tmp_path / "B.pdb").write_text("END\n", encoding="utf-8")

        pairs_df = pd.DataFrame({
            "proteinA": ["A"],
            "proteinB": ["B"],
            "label": ["interaction"],
        })

        seen = {}

        def fake_runner(rec_pdb, lig_pdb, output_root, pair_label, **kwargs):
            seen["output_root"] = output_root
            seen["kwargs"] = kwargs
            return str(tmp_path / "fake_run")

        plugin = SimpleNamespace(runner=fake_runner)
        monkeypatch.setattr(batch_dock.registry, "get", lambda _name: plugin)

        results = batch_dock.batch_dock(
            pairs_df,
            engines=["lightdock"],
            pdb_dir=str(tmp_path),
            output_root=str(tmp_path / "out"),
            engine_kwargs={"lightdock": {"cores": 8, "steps": 25}},
        )

        assert len(results) == 1
        assert results.iloc[0]["status"] == "ok"
        assert seen["kwargs"]["cores"] == 8
        assert seen["kwargs"]["steps"] == 25
        assert seen["output_root"] == str((tmp_path / "out").resolve())

    def test_resume_skips_successful_engine(self, tmp_path, monkeypatch):
        from ppinsight import batch_dock

        (tmp_path / "A.pdb").write_text("END\n", encoding="utf-8")
        (tmp_path / "B.pdb").write_text("END\n", encoding="utf-8")
        pairs_df = pd.DataFrame({
            "proteinA": ["A"],
            "proteinB": ["B"],
            "label": ["interaction"],
        })

        def fail_if_called(*_args, **_kwargs):
            raise AssertionError("resumed engine must not run")

        plugin = SimpleNamespace(runner=fail_if_called)
        monkeypatch.setattr(batch_dock.registry, "get", lambda _name: plugin)

        results = batch_dock.batch_dock(
            pairs_df,
            engines=["lightdock"],
            pdb_dir=str(tmp_path),
            completed_runs={("A", "B", "lightdock")},
        )

        assert results.empty

    def test_runner_exception_records_structured_failure(self, tmp_path, monkeypatch):
        from ppinsight import batch_dock, registry

        (tmp_path / "A.pdb").write_text(
            "ATOM      1  CA  GLY A   1\n", encoding="utf-8"
        )
        (tmp_path / "B.pdb").write_text(
            "ATOM      1  CA  GLY B   1\n", encoding="utf-8"
        )
        failed_dir = tmp_path / "output" / "lightdock_runs" / "A_vs_B"
        failed_dir.mkdir(parents=True)
        log_path = failed_dir / "lightdock.log"
        log_path.write_text("failure details\n", encoding="utf-8")
        pairs_df = pd.DataFrame({"proteinA": ["A"], "proteinB": ["B"]})

        def failed_runner(*_args, **_kwargs):
            raise registry.EngineRunError(
                "LightDock did not produce outputs.",
                output_dir=str(failed_dir),
                log_path=str(log_path),
                error_type="LightDockSimulationError",
            )

        monkeypatch.setattr(
            batch_dock.registry,
            "get",
            lambda _name: SimpleNamespace(runner=failed_runner),
        )

        results = batch_dock.batch_dock(
            pairs_df,
            engines=["lightdock"],
            pdb_dir=str(tmp_path),
            output_root=str(tmp_path / "output"),
        )

        row = results.iloc[0]
        assert row["status"] == "failed"
        assert row["error_type"] == "LightDockSimulationError"
        assert row["error_message"] == "LightDock did not produce outputs."
        assert row["output_dir"] == str(failed_dir)
        assert row["log_path"] == str(log_path)

    def test_clean_failed_removes_recorded_directory_before_retry(
        self, tmp_path, monkeypatch
    ):
        from ppinsight import batch_dock

        (tmp_path / "A.pdb").write_text(
            "ATOM      1  CA  GLY A   1\n", encoding="utf-8"
        )
        (tmp_path / "B.pdb").write_text(
            "ATOM      1  CA  GLY B   1\n", encoding="utf-8"
        )
        output_root = tmp_path / "output"
        failed_dir = output_root / "lightdock_runs" / "A_vs_B"
        failed_dir.mkdir(parents=True)
        (failed_dir / "partial.txt").write_text("partial\n", encoding="utf-8")
        pairs_df = pd.DataFrame({"proteinA": ["A"], "proteinB": ["B"]})

        def successful_runner(*_args, **_kwargs):
            assert not failed_dir.exists()
            return str(output_root / "lightdock_runs" / "A_vs_B_retry")

        monkeypatch.setattr(
            batch_dock.registry,
            "get",
            lambda _name: SimpleNamespace(runner=successful_runner),
        )

        results = batch_dock.batch_dock(
            pairs_df,
            engines=["lightdock"],
            pdb_dir=str(tmp_path),
            output_root=str(output_root),
            failed_run_dirs={("A", "B", "lightdock"): str(failed_dir)},
            clean_failed=True,
        )

        assert results.iloc[0]["status"] == "ok"
        assert not failed_dir.exists()

    def test_upsert_results_replaces_failed_row_after_retry(self):
        from ppinsight import batch_dock

        existing = pd.DataFrame([
            batch_dock._result_row(
                "A", "B", "", "", "haddock", "/tmp/failed", "failed",
                error_type="RuntimeError",
            )
        ])
        retry = pd.DataFrame([
            batch_dock._result_row(
                "A", "B", "", "", "haddock", "/tmp/success", "ok"
            )
        ])

        result = batch_dock._upsert_results(existing, retry)

        assert len(result) == 1
        assert result.iloc[0]["status"] == "ok"
        assert result.iloc[0]["output_dir"] == "/tmp/success"

    def test_cli_resume_upserts_successful_retry(self, tmp_path, monkeypatch):
        from ppinsight import batch_dock, registry

        pairs_path = tmp_path / "pairs.csv"
        pairs_path.write_text("proteinA,proteinB\nA,B\n", encoding="utf-8")
        (tmp_path / "A.pdb").write_text(
            "ATOM      1  CA  GLY A   1\n", encoding="utf-8"
        )
        (tmp_path / "B.pdb").write_text(
            "ATOM      1  CA  GLY B   1\n", encoding="utf-8"
        )
        output_root = tmp_path / "output"
        manifest = output_root / "scores" / "batch_results.csv"
        failed_dir = output_root / "lightdock_runs" / "A_vs_B"
        failed_dir.mkdir(parents=True)

        def failed_runner(*_args, **_kwargs):
            raise registry.EngineRunError(
                "initial failure",
                output_dir=str(failed_dir),
                error_type="RuntimeError",
            )

        plugin = SimpleNamespace(runner=failed_runner)
        monkeypatch.setattr(batch_dock.registry, "get", lambda _name: plugin)
        first_args = [
            str(pairs_path),
            "--engines",
            "lightdock",
            "--pdb-dir",
            str(tmp_path),
            "--output-root",
            str(output_root),
            "-o",
            str(manifest),
        ]
        batch_dock.main(first_args)

        plugin.runner = lambda *_args, **_kwargs: str(
            output_root / "lightdock_runs" / "A_vs_B_retry"
        )
        batch_dock.main(first_args + ["--resume"])

        results = pd.read_csv(manifest)
        assert len(results) == 1
        assert results.iloc[0]["status"] == "ok"
        assert results.iloc[0]["output_dir"].endswith("A_vs_B_retry")

    def test_preflight_reports_lightdock_heteroatom_warning(
        self, tmp_path, monkeypatch
    ):
        from ppinsight import batch_dock

        (tmp_path / "A.pdb").write_text(
            "ATOM      1  CA  GLY A   1\nHETATM    2  S   SO4 A   2\n",
            encoding="utf-8",
        )
        (tmp_path / "B.pdb").write_text(
            "ATOM      1  CA  GLY B   1\n", encoding="utf-8"
        )
        monkeypatch.setattr(batch_dock.shutil, "which", lambda _name: "/mock/bin")
        pairs_df = pd.DataFrame({"proteinA": ["A"], "proteinB": ["B"]})

        results = batch_dock.preflight_batch(
            pairs_df,
            engines=["lightdock"],
            pdb_dir=str(tmp_path),
            engine_kwargs={"lightdock": {"auto_clean_pdb": True}},
        )

        row = results.iloc[0]
        assert row["status"] == "preflight_ok"
        assert "HETATM" in row["preflight_warnings"]
        assert "Auto-clean retry is enabled" in row["preflight_warnings"]

    def test_cli_defaults_results_to_output_root(self, tmp_path, capsys):
        from ppinsight import batch_dock

        pairs_path = tmp_path / "pairs.csv"
        pairs_path.write_text(
            "proteinA,proteinB\nA,B\n", encoding="utf-8"
        )
        output_root = tmp_path / "custom_output"

        batch_dock.main([
            str(pairs_path),
            "--dry-run",
            "--output-root",
            str(output_root),
        ])

        assert (output_root / "scores" / "batch_results.csv").is_file()
        assert "Total batch runtime:" in capsys.readouterr().out

    def test_haddock_kwargs_include_refine_controls(self):
        from ppinsight import batch_dock

        args = SimpleNamespace(
            lightdock_steps=100,
            lightdock_swarms=400,
            lightdock_glowworms=200,
            cores=4,
            lightdock_anm=False,
            lightdock_scoring=None,
            lightdock_auto_clean_pdb=False,
            haddock_sampling=200,
            haddock_select_top=50,
            haddock_tolerance=25,
            haddock_skip_refinement=True,
            haddock_skip_flexref=False,
            haddock_skip_emref=False,
            rosetta_n_runs=5,
            rosetta_top_n=20,
            rosetta_relax=False,
            rosetta_no_cluster=False,
            rosetta_cluster_top_n=200,
            rosetta_rmsd_cutoff=4.0,
            rosetta_no_auto_filter=False,
            rosetta_debug_pyrosetta=False,
            rosetta_save_top=0,
        )

        kw = batch_dock._engine_kwargs_from_args(args)["haddock"]
        assert kw["sampling"] == 200
        assert kw["select_top"] == 50
        assert kw["tolerance"] == 25
        assert kw["skip_flexref"] is True
        assert kw["skip_emref"] is True

    def test_screening_preset_uses_reduced_end_to_end_values(self):
        from ppinsight import batch_dock

        args = SimpleNamespace(
            screening=True,
            lightdock_steps=100,
            lightdock_swarms=400,
            lightdock_glowworms=200,
            lightdock_anm=True,
            lightdock_auto_clean_pdb=False,
            haddock_sampling=10000,
            haddock_select_top=400,
            haddock_skip_refinement=False,
            haddock_skip_flexref=False,
            haddock_skip_emref=False,
            rosetta_n_runs=5000,
            rosetta_top_n=20,
            rosetta_cluster_top_n=200,
            rosetta_relax=True,
        )

        batch_dock._apply_screening_preset(args, ["--screening"])

        assert args.lightdock_steps == 50
        assert args.lightdock_swarms == 50
        assert args.lightdock_glowworms == 50
        assert args.lightdock_anm is False
        assert args.lightdock_auto_clean_pdb is True
        assert args.haddock_sampling == 1000
        assert args.haddock_select_top == 100
        assert args.haddock_skip_refinement is True
        assert args.rosetta_n_runs == 100
        assert args.rosetta_top_n == 20
        assert args.rosetta_cluster_top_n == 100
        assert args.rosetta_relax is False

    def test_screening_preset_preserves_explicit_overrides(self):
        from ppinsight import batch_dock

        args = SimpleNamespace(
            screening=True,
            lightdock_steps=300,
            lightdock_swarms=400,
            lightdock_glowworms=200,
            lightdock_anm=True,
            lightdock_auto_clean_pdb=False,
            haddock_sampling=3000,
            haddock_select_top=400,
            haddock_skip_refinement=False,
            haddock_skip_flexref=False,
            haddock_skip_emref=False,
            rosetta_n_runs=500,
            rosetta_top_n=20,
            rosetta_cluster_top_n=200,
            rosetta_relax=True,
        )

        batch_dock._apply_screening_preset(
            args,
            [
                "--screening",
                "--lightdock-steps",
                "300",
                "--haddock-sampling=3000",
                "--rosetta-n-runs",
                "500",
                "--lightdock-anm",
                "--rosetta-relax",
            ],
        )

        assert args.lightdock_steps == 300
        assert args.haddock_sampling == 3000
        assert args.rosetta_n_runs == 500
        assert args.lightdock_anm is True
        assert args.rosetta_relax is True

    def test_haddock_skip_flexref_implies_skip_emref(self):
        from ppinsight import batch_dock

        args = SimpleNamespace(
            lightdock_steps=100,
            lightdock_swarms=400,
            lightdock_glowworms=200,
            cores=4,
            lightdock_anm=False,
            lightdock_scoring=None,
            lightdock_auto_clean_pdb=False,
            haddock_sampling=200,
            haddock_select_top=50,
            haddock_tolerance=5,
            haddock_skip_refinement=False,
            haddock_skip_flexref=True,
            haddock_skip_emref=False,
            rosetta_n_runs=5,
            rosetta_top_n=20,
            rosetta_relax=False,
            rosetta_no_cluster=False,
            rosetta_cluster_top_n=200,
            rosetta_rmsd_cutoff=4.0,
            rosetta_no_auto_filter=False,
            rosetta_debug_pyrosetta=False,
            rosetta_save_top=0,
        )

        kw = batch_dock._engine_kwargs_from_args(args)["haddock"]
        assert kw["skip_flexref"] is True
        assert kw["skip_emref"] is True


# ---------------------------------------------------------------------------
# Visualizer label-aware features
# ---------------------------------------------------------------------------

class TestCompareScoresByLabel:
    @pytest.fixture
    def labeled_scores(self):
        return pd.DataFrame({
            "model": ["lightdock"] * 4 + ["haddock"] * 4,
            "score_type": ["dockq"] * 8,
            "score_value": [0.8, 0.7, 0.3, 0.2, 0.9, 0.85, 0.4, 0.25],
            "proteinA": ["A", "A", "B", "B"] * 2,
            "proteinB": ["X", "X", "Y", "Y"] * 2,
            "label": ["interaction", "interaction",
                       "non-interaction", "non-interaction"] * 2,
        })

    def test_basic(self, labeled_scores):
        import matplotlib

        from ppinsight.visualizer import compare_scores_by_label
        matplotlib.use("Agg")
        fig = compare_scores_by_label(labeled_scores, "dockq", output=None)
        # Should produce a figure without error
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_save_to_file(self, labeled_scores, tmp_path):
        import matplotlib

        from ppinsight.visualizer import compare_scores_by_label
        matplotlib.use("Agg")
        out = str(tmp_path / "labeled.png")
        fig = compare_scores_by_label(labeled_scores, "dockq", output=out)
        assert os.path.isfile(out)
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_no_label_column(self):
        import matplotlib

        from ppinsight.visualizer import compare_scores_by_label
        matplotlib.use("Agg")
        df = pd.DataFrame({
            "model": ["lightdock"],
            "score_type": ["dockq"],
            "score_value": [0.8],
        })
        with pytest.raises(ValueError, match="label"):
            compare_scores_by_label(df, "dockq")


class TestClassificationSummary:
    @pytest.fixture
    def labeled_scores(self):
        return pd.DataFrame({
            "model": ["lightdock"] * 4 + ["haddock"] * 4,
            "score_type": ["dockq"] * 8,
            "score_value": [0.8, 0.7, 0.3, 0.2, 0.9, 0.85, 0.4, 0.25],
            "proteinA": ["A", "A", "B", "B"] * 2,
            "proteinB": ["X", "X", "Y", "Y"] * 2,
            "label": ["interaction", "interaction",
                       "non-interaction", "non-interaction"] * 2,
        })

    def test_basic(self, labeled_scores):
        from ppinsight.visualizer import classification_summary
        result = classification_summary(labeled_scores, "dockq", threshold=0.5)
        assert "model" in result.columns
        assert "TP" in result.columns
        assert "accuracy" in result.columns
        assert len(result) == 2  # two models

    def test_auto_threshold(self, labeled_scores):
        from ppinsight.visualizer import classification_summary
        result = classification_summary(labeled_scores, "dockq")
        assert result["threshold"].iloc[0] is not None

    def test_lower_is_better(self, labeled_scores):
        from ppinsight.visualizer import classification_summary
        result = classification_summary(
            labeled_scores, "dockq", threshold=0.5, higher_is_better=False
        )
        assert len(result) == 2

    def test_auto_direction_dockq(self, labeled_scores):
        """DockQ is higher-is-better; auto-detection should pick that up."""
        from ppinsight.visualizer import classification_summary
        # Don't pass higher_is_better — let it auto-detect.
        # With threshold=0.5 and dockq (higher=better), interaction pairs
        # (scores 0.8, 0.7 / 0.9, 0.85 > 0.5) should be TP,
        # non-interaction pairs (0.3, 0.2 / 0.4, 0.25 < 0.5) should be TN.
        result = classification_summary(labeled_scores, "dockq", threshold=0.5)
        for _, row in result.iterrows():
            assert row["TP"] == 1   # one interaction pair, correctly classified
            assert row["TN"] == 1   # one non-interaction pair, correctly classified
            assert row["FP"] == 0
            assert row["FN"] == 0

    def test_auto_direction_haddock_score(self):
        """HADDOCK 'score' is lower-is-better; auto-detect must flip direction."""
        from ppinsight.visualizer import classification_summary
        df = pd.DataFrame({
            "model": ["haddock"] * 4,
            "score_type": ["score"] * 4,
            # Lower score = better binding.  Interacting pair should have
            # the lower (more negative) score.
            "score_value": [-120.0, -110.0, -20.0, -10.0],
            "proteinA": ["A", "A", "B", "B"],
            "proteinB": ["X", "X", "Y", "Y"],
            "label": ["interaction", "interaction",
                       "non-interaction", "non-interaction"],
        })
        # threshold at -65 (median).  With lower-is-better auto-detected:
        # interaction pair mean = -115 ≤ -65 → predicted interacting (TP)
        # non-interaction pair mean = -15 > -65 → predicted non-interacting (TN)
        result = classification_summary(df, "score", threshold=-65.0)
        assert result["TP"].iloc[0] == 1
        assert result["TN"].iloc[0] == 1
        assert result["FP"].iloc[0] == 0
        assert result["FN"].iloc[0] == 0


import os  # noqa: E402
