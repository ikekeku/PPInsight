"""Tests for the lightweight Rosetta score-analysis helpers."""

import sys
import types

import pandas as pd

from ppinsight.rosetta.analyze import (
    analyze_scores,
    cluster_and_rank,
    export_scores_to_csv,
)


def test_analyze_scores_prefers_i_sc_over_total_score():
    results = [
        {"run": 1, "score": -140.0, "total_score": -140.0, "i_sc": -9.0},
        {"run": 2, "score": -110.0, "total_score": -110.0, "i_sc": -15.0},
    ]

    analysis = analyze_scores(results, top_n=1)

    assert analysis["primary_metric"] == "i_sc"
    assert analysis["top_results"][0]["run"] == 2
    assert analysis["final_score"] == -15.0


def test_export_scores_to_csv_writes_explicit_rosetta_schema(tmp_path):
    results = [
        {"run": 1, "score": -140.0, "total_score": -140.0, "i_sc": -9.0},
        {"run": 2, "score": -110.0, "total_score": -110.0, "i_sc": -15.0},
    ]
    output_path = tmp_path / "docking_scores.csv"

    export_scores_to_csv(results, output_path)

    assert output_path.read_text(encoding="utf-8") == (
        "run,description,total_score,i_sc\n"
        "1,decoy_1,-140.0,-9.0\n"
        "2,decoy_2,-110.0,-15.0\n"
    )


def test_cluster_and_rank_defaults_to_i_sc(monkeypatch, tmp_path):
    captured = {}

    def fake_cluster_decoys(df, pdb_dir, score_col, top_n, rmsd_cutoff):
        captured["columns"] = list(df.columns)
        captured["pdb_dir"] = str(pdb_dir)
        captured["score_col"] = score_col
        captured["top_n"] = top_n
        captured["rmsd_cutoff"] = rmsd_cutoff
        return pd.DataFrame(
            {
                "description": ["decoy_1"],
                "total_score": [-100.0],
                "i_sc": [-12.0],
                "cluster": [0],
                "cluster_size": [1],
                "cluster_rank": [0],
                "is_top_of_cluster": [True],
            }
        )

    monkeypatch.setitem(
        sys.modules,
        "ppinsight.rosetta.cluster",
        types.SimpleNamespace(cluster_decoys=fake_cluster_decoys),
    )

    scores_csv = tmp_path / "scores.csv"
    scores_csv.write_text(
        "description,total_score,i_sc\n"
        "decoy_1,-100.0,-12.0\n"
        "decoy_2,-90.0,-8.0\n",
        encoding="utf-8",
    )

    result = cluster_and_rank(scores_csv, tmp_path)

    assert captured["score_col"] == "i_sc"
    assert result is not None
