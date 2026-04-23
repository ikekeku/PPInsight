#!/usr/bin/env python3
"""Generate ALL possible PPInsight plot examples with fabricated data.

Every dataset and plot produced by this script is EXPLICITLY FABRICATED
-- synthetic numbers created for demonstration purposes only.  No real
docking data was used.  Do not cite these plots as experimental results.

The synthetic score tables intentionally keep the same lightweight
traceability columns used by real collected scores (`run_id`, `pose_id`,
`output_path`) so plot development exercises the same row structure.

Usage
-----
    cd PPInsight
    python examples/generate_example_plots.py

Outputs
-------
    examples/example_plots/fabricated_scores.tsv
        Standalone long-format fabricated scores table used to generate
        the plot gallery.
    examples/example_plots/*.png
        Plot examples rendered from that fabricated scores table.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

from ppinsight.visualizer import (
    DEFAULT_THEME,
    apply_theme,
    cdf_plot,
    model_agreement_scatter,
    pairwise_difference_plot,
    quality_bar_chart,
    ridge_plot,
    roc_curve_plot,
    violin_plot,
)

ROOT = Path(__file__).resolve().parent.parent

OUT = ROOT / "examples" / "example_plots"
OUT_SCORES = OUT / "fabricated_scores.tsv"

RNG_SEED = 42

ENGINES = ["HADDOCK", "LightDock", "Rosetta"]
PAIRS = [
    ("2UUY_rec", "2UUY_lig"),
    ("3HFM_H", "3HFM_Y"),
    ("1BRS_barnase", "1BRS_barstar"),
    ("1PPE_E", "1PPE_I"),
]
N_POSES = 50

_QUALITY_PARAMS = {
    "dockq": [
        {"HADDOCK": (0.85, 0.08), "LightDock": (0.82, 0.10), "Rosetta": (0.83, 0.09)},
        {"HADDOCK": (0.55, 0.10), "LightDock": (0.50, 0.12), "Rosetta": (0.52, 0.09)},
        {"HADDOCK": (0.30, 0.06), "LightDock": (0.28, 0.08), "Rosetta": (0.32, 0.07)},
        {"HADDOCK": (0.10, 0.05), "LightDock": (0.08, 0.06), "Rosetta": (0.12, 0.05)},
    ],
    "fnat": [
        {"HADDOCK": (0.80, 0.10), "LightDock": (0.72, 0.12), "Rosetta": (0.75, 0.10)},
        {"HADDOCK": (0.48, 0.10), "LightDock": (0.40, 0.12), "Rosetta": (0.45, 0.09)},
        {"HADDOCK": (0.22, 0.08), "LightDock": (0.18, 0.09), "Rosetta": (0.25, 0.07)},
        {"HADDOCK": (0.06, 0.04), "LightDock": (0.05, 0.04), "Rosetta": (0.07, 0.04)},
    ],
    "irmsd": [
        {"HADDOCK": (1.2, 0.5), "LightDock": (1.5, 0.7), "Rosetta": (1.3, 0.6)},
        {"HADDOCK": (3.5, 1.0), "LightDock": (4.0, 1.3), "Rosetta": (3.8, 1.1)},
        {"HADDOCK": (6.5, 1.5), "LightDock": (7.5, 2.0), "Rosetta": (6.0, 1.4)},
        {"HADDOCK": (10.0, 2.5), "LightDock": (11.0, 3.0), "Rosetta": (9.5, 2.2)},
    ],
    "lrmsd": [
        {"HADDOCK": (2.5, 1.0), "LightDock": (3.0, 1.5), "Rosetta": (2.8, 1.2)},
        {"HADDOCK": (7.0, 2.0), "LightDock": (8.5, 2.5), "Rosetta": (7.5, 2.0)},
        {"HADDOCK": (13.0, 3.0), "LightDock": (15.0, 4.0), "Rosetta": (12.0, 3.0)},
        {"HADDOCK": (22.0, 5.0), "LightDock": (24.0, 6.0), "Rosetta": (20.0, 4.5)},
    ],
}

_ENGINE_SCORE_PARAMS = {
    "score": [
        {"HADDOCK": (-105.0, 20.0)},
        {"HADDOCK": (-55.0, 18.0)},
        {"HADDOCK": (-20.0, 12.0)},
        {"HADDOCK": (-8.0, 8.0)},
    ],
    "luciferin_score": [
        {"LightDock": (18.0, 2.5)},
        {"LightDock": (14.0, 3.0)},
        {"LightDock": (9.0, 6.0)},
        {"LightDock": (-24.0, 18.0)},
    ],
    "total_score": [
        {"Rosetta": (-193.0, 10.0)},
        {"Rosetta": (-175.0, 15.0)},
        {"Rosetta": (-150.0, 18.0)},
        {"Rosetta": (-120.0, 20.0)},
    ],
}


def _traceability_fields(engine: str, protA: str, protB: str, pose_idx: int):
    run_id = f"{engine.lower()}_{protA}_{protB}_demo"
    pose_id = f"{engine.lower()}_{protA}_{protB}_pose_{pose_idx:03d}"
    output_path = f"/demo/{engine.lower()}/{protA}_vs_{protB}/{pose_id}.pdb"
    return run_id, pose_id, output_path


def _subset_scores(df, engines=None, pairs=None, metrics=None, labels=None):
    out = df
    if engines is not None:
        out = out[out["model"].isin(engines)]
    if pairs is not None:
        pair_set = set(pairs)
        pair_mask = [
            (protein_a, protein_b) in pair_set
            for protein_a, protein_b in zip(
                out["proteinA"], out["proteinB"], strict=False
            )
        ]
        out = out[pair_mask]
    if metrics is not None:
        out = out[out["score_type"].isin(metrics)]
    if labels is not None:
        out = out[out["label"].isin(labels)]
    return out.copy()


def _build_quality_scores(engines, pairs, metrics, rng):
    rows = []
    for pair_idx, (protA, protB) in enumerate(pairs):
        for eng in engines:
            for metric in metrics:
                if metric not in _QUALITY_PARAMS:
                    continue
                param_list = _QUALITY_PARAMS[metric]
                if pair_idx >= len(param_list):
                    continue
                params = param_list[pair_idx]
                if eng not in params:
                    continue
                mu, sigma = params[eng]
                values = rng.normal(mu, sigma, N_POSES)
                if metric in {"dockq", "fnat"}:
                    values = np.clip(values, 0.0, 1.0)
                elif metric in {"irmsd", "lrmsd"}:
                    values = np.clip(values, 0.0, None)
                for pose_idx, v in enumerate(values, start=1):
                    run_id, pose_id, output_path = _traceability_fields(
                        eng, protA, protB, pose_idx,
                    )
                    rows.append({
                        "model": eng,
                        "score_type": metric,
                        "score_value": float(v),
                        "proteinA": protA,
                        "proteinB": protB,
                        "run_id": run_id,
                        "pose_id": pose_id,
                        "output_path": output_path,
                    })
    return pd.DataFrame(rows)


def _build_engine_scores(engines, pairs, rng):
    engine_metric = {
        "HADDOCK": "score",
        "LightDock": "luciferin_score",
        "Rosetta": "total_score",
    }
    rows = []
    for pair_idx, (protA, protB) in enumerate(pairs):
        for eng in engines:
            metric = engine_metric.get(eng)
            if not metric or metric not in _ENGINE_SCORE_PARAMS:
                continue
            param_list = _ENGINE_SCORE_PARAMS[metric]
            if pair_idx >= len(param_list):
                continue
            params = param_list[pair_idx]
            if eng not in params:
                continue
            mu, sigma = params[eng]
            values = rng.normal(mu, sigma, N_POSES)
            for pose_idx, v in enumerate(values, start=1):
                run_id, pose_id, output_path = _traceability_fields(
                    eng, protA, protB, pose_idx,
                )
                rows.append({
                    "model": eng,
                    "score_type": metric,
                    "score_value": float(v),
                    "proteinA": protA,
                    "proteinB": protB,
                    "run_id": run_id,
                    "pose_id": pose_id,
                    "output_path": output_path,
                })
    return pd.DataFrame(rows)


def _build_roc_data(engines, pairs_binder, pairs_nonbinder, rng, metric="dockq"):
    rows = []
    binder_params = {
        "dockq": {
            "HADDOCK": (0.45, 0.18), "LightDock": (0.42, 0.20),
            "Rosetta": (0.48, 0.19),
        },
    }
    nonbinder_params = {
        "dockq": {
            "HADDOCK": (0.28, 0.16), "LightDock": (0.26, 0.18),
            "Rosetta": (0.30, 0.17),
        },
    }
    for eng in engines:
        if metric not in binder_params or eng not in binder_params[metric]:
            continue
        mu_b, sig_b = binder_params[metric][eng]
        for pose_idx, (protA, protB) in enumerate(pairs_binder, start=1):
            v = float(np.clip(rng.normal(mu_b, sig_b), 0.0, 1.0))
            run_id, pose_id, output_path = _traceability_fields(
                eng, protA, protB, pose_idx,
            )
            rows.append({
                "model": eng, "score_type": metric,
                "score_value": v,
                "proteinA": protA, "proteinB": protB,
                "label": "interaction",
                "run_id": run_id,
                "pose_id": pose_id,
                "output_path": output_path,
            })
        mu_n, sig_n = nonbinder_params[metric][eng]
        for pose_idx, (protA, protB) in enumerate(pairs_nonbinder, start=1):
            v = float(np.clip(rng.normal(mu_n, sig_n), 0.0, 1.0))
            run_id, pose_id, output_path = _traceability_fields(
                eng, protA, protB, pose_idx,
            )
            rows.append({
                "model": eng, "score_type": metric,
                "score_value": v,
                "proteinA": protA, "proteinB": protB,
                "label": "non-interaction",
                "run_id": run_id,
                "pose_id": pose_id,
                "output_path": output_path,
            })
    return pd.DataFrame(rows)


def _build_prodigy_scores(engines, pairs, rng):
    rows = []
    prodigy_params = {
        "HADDOCK": (-9.5, 2.0),
        "LightDock": (-7.8, 2.5),
        "Rosetta": (-8.6, 1.8),
    }
    for eng in engines:
        mu, sigma = prodigy_params[eng]
        for protA, protB in pairs:
            for pose_idx, v in enumerate(rng.normal(mu, sigma, N_POSES), start=1):
                run_id, pose_id, output_path = _traceability_fields(
                    eng, protA, protB, pose_idx,
                )
                rows.append({
                    "model": eng,
                    "score_type": "prodigy_ddg",
                    "score_value": float(v),
                    "proteinA": protA,
                    "proteinB": protB,
                    "run_id": run_id,
                    "pose_id": pose_id,
                    "output_path": output_path,
                })
    return pd.DataFrame(rows)


def build_fabricated_scores():
    rng = np.random.default_rng(RNG_SEED)
    binder_pairs = [(f"bind_rec_{i:02d}", f"bind_lig_{i:02d}") for i in range(50)]
    nonbinder_pairs = [
        (f"decoy_rec_{i:02d}", f"decoy_lig_{i:02d}") for i in range(50)
    ]

    frames = [
        _build_quality_scores(
            ENGINES, PAIRS, ["dockq", "fnat", "irmsd", "lrmsd"], rng,
        ),
        _build_engine_scores(ENGINES, PAIRS, rng),
        _build_roc_data(ENGINES, binder_pairs, nonbinder_pairs, rng),
        _build_prodigy_scores(ENGINES, PAIRS[:2], rng),
    ]

    return pd.concat(frames, ignore_index=True, sort=False)


def _write_fabricated_scores(df):
    df.to_csv(OUT_SCORES, sep="\t", index=False)
    print(f"  + {OUT_SCORES.relative_to(ROOT)}")


def _save(fig, name):
    import matplotlib.pyplot as plt
    path = OUT / name
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  + {path.relative_to(ROOT)}")


def gen_violin(scores_df):
    print("\n-- Violin plots --")
    for n_eng in (1, 2, 3):
        for n_pair in (1, 2, 3):
            df = _subset_scores(
                scores_df,
                engines=ENGINES[:n_eng],
                pairs=PAIRS[:n_pair],
                metrics=["dockq"],
            )
            fig = violin_plot(df, metric="dockq", output=None)
            _save(fig, f"violin_{n_eng}eng_{n_pair}pair_dockq.png")
    for eng in ENGINES:
        edf = _subset_scores(scores_df, engines=[eng], pairs=PAIRS[:1])
        edf = edf[edf["label"].isna() if "label" in edf.columns else True]
        for metric in edf["score_type"].unique():
            fig = violin_plot(edf, metric=metric, output=None)
            _save(fig, f"violin_1eng_1pair_{metric}.png")


def gen_ridge(scores_df):
    print("\n-- Ridge plots --")
    for n_eng in (1, 2, 3):
        for n_pair in (1, 2):
            df = _subset_scores(
                scores_df,
                engines=ENGINES[:n_eng],
                pairs=PAIRS[:n_pair],
                metrics=["dockq"],
            )
            fig = ridge_plot(df, metric="dockq", output=None)
            _save(fig, f"ridge_{n_eng}eng_{n_pair}pair_dockq.png")


def gen_scatter(scores_df):
    print("\n-- Scatter plots --")
    for eng_a, eng_b in [
        ("HADDOCK", "LightDock"),
        ("HADDOCK", "Rosetta"),
        ("LightDock", "Rosetta"),
    ]:
        df = _subset_scores(
            scores_df,
            engines=[eng_a, eng_b],
            pairs=PAIRS[:3],
            metrics=["dockq"],
        )
        fig = model_agreement_scatter(
            df, metric="dockq", model_x=eng_a, model_y=eng_b, output=None,
        )
        _save(fig, f"scatter_{eng_a}_vs_{eng_b}_dockq.png")


def gen_difference(scores_df):
    print("\n-- Difference plots --")
    for eng_a, eng_b in [("HADDOCK", "LightDock"), ("HADDOCK", "Rosetta")]:
        df = _subset_scores(
            scores_df,
            engines=[eng_a, eng_b],
            pairs=PAIRS[:3],
            metrics=["dockq"],
        )
        fig = pairwise_difference_plot(
            df, metric="dockq", model_a=eng_a, model_b=eng_b, output=None,
        )
        _save(fig, f"difference_{eng_a}_vs_{eng_b}_dockq.png")


def gen_roc(scores_df):
    print("\n-- ROC plots --")
    for n_eng in (1, 2, 3):
        df = _subset_scores(
            scores_df,
            engines=ENGINES[:n_eng],
            metrics=["dockq"],
            labels=["interaction", "non-interaction"],
        )
        fig = roc_curve_plot(df, metric="dockq", output=None)
        _save(fig, f"roc_{n_eng}eng_dockq.png")


def gen_quality_bar(scores_df):
    print("\n-- Quality bar charts --")
    quality_metrics = ["dockq", "fnat", "irmsd", "lrmsd"]
    for n_eng in (1, 2, 3):
        for n_pair in (1, 2, 3, 4):
            df = _subset_scores(
                scores_df,
                engines=ENGINES[:n_eng],
                pairs=PAIRS[:n_pair],
                metrics=quality_metrics,
            )
            fig = quality_bar_chart(df, output=None)
            _save(fig, f"quality_bar_{n_eng}eng_{n_pair}pair.png")


def gen_cdf(scores_df):
    print("\n-- CDF plots --")
    for n_eng in (1, 2, 3):
        df = _subset_scores(
            scores_df,
            engines=ENGINES[:n_eng],
            pairs=PAIRS[:2],
            metrics=["dockq"],
        )
        fig = cdf_plot(df, metric="dockq", output=None)
        _save(fig, f"cdf_{n_eng}eng_dockq.png")

        pdf = _subset_scores(
            scores_df,
            engines=ENGINES[:n_eng],
            pairs=PAIRS[:2],
            metrics=["prodigy_ddg"],
        )
        fig = cdf_plot(pdf, metric="prodigy_ddg", output=None)
        _save(fig, f"cdf_{n_eng}eng_prodigy_ddg.png")


def main():
    matplotlib.use("Agg")
    apply_theme(DEFAULT_THEME)
    OUT.mkdir(parents=True, exist_ok=True)
    print(
        "================================================================\n"
        "  PPInsight -- Example Plot Generator (FABRICATED DATA)\n"
        "================================================================"
    )
    print(f"Output directory: {OUT.relative_to(ROOT)}/")
    scores_df = build_fabricated_scores()
    _write_fabricated_scores(scores_df)
    gen_violin(scores_df)
    gen_ridge(scores_df)
    gen_scatter(scores_df)
    gen_difference(scores_df)
    gen_roc(scores_df)
    gen_quality_bar(scores_df)
    gen_cdf(scores_df)
    n_files = len(list(OUT.glob("*.png")))
    print(f"Fabricated scores rows: {len(scores_df)}")
    print(f"\nDone -- {n_files} plots saved to {OUT.relative_to(ROOT)}/")
    print("\nWARNING: ALL DATA IS FABRICATED. For demonstration only.")


if __name__ == "__main__":
    main()
