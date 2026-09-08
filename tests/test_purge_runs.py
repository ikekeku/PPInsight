"""Tests for manifest-driven cleanup of failed batch run directories."""

import pandas as pd

from ppinsight.purge_runs import failed_run_candidates, purge_failed_runs


def test_purge_failed_runs_is_dry_by_default(tmp_path):
    output_root = tmp_path / "output"
    failed_dir = output_root / "lightdock_runs" / "A_vs_B"
    failed_dir.mkdir(parents=True)
    (failed_dir / "partial.txt").write_text("partial\n", encoding="utf-8")
    manifest = tmp_path / "results.csv"
    pd.DataFrame([
        {
            "proteinA": "A",
            "proteinB": "B",
            "engine": "lightdock",
            "output_dir": str(failed_dir),
            "status": "failed",
        }
    ]).to_csv(manifest, index=False)

    removed, skipped = purge_failed_runs(manifest, output_root)

    assert removed == 1
    assert skipped == 0
    assert failed_dir.is_dir()


def test_purge_failed_runs_removes_only_safe_paths(tmp_path):
    output_root = tmp_path / "output"
    failed_dir = output_root / "haddock_runs" / "A_vs_B"
    failed_dir.mkdir(parents=True)
    unsafe_dir = tmp_path / "outside"
    unsafe_dir.mkdir()
    manifest = tmp_path / "results.csv"
    pd.DataFrame([
        {
            "proteinA": "A",
            "proteinB": "B",
            "engine": "haddock",
            "output_dir": str(failed_dir),
            "status": "failed",
        },
        {
            "proteinA": "C",
            "proteinB": "D",
            "engine": "rosetta",
            "output_dir": str(unsafe_dir),
            "status": "failed",
        },
    ]).to_csv(manifest, index=False)

    candidates = failed_run_candidates(manifest, output_root)
    removed, skipped = purge_failed_runs(manifest, output_root, execute=True)

    assert len(candidates) == 2
    assert removed == 1
    assert skipped == 1
    assert not failed_dir.exists()
    assert unsafe_dir.exists()
