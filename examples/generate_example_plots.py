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
    funnel_plot,
    model_agreement_scatter,
    quality_bar_chart,
    roc_curve_plot,
    strip_plot,
    topn_plot,
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

# Separate target for the funnel example: engine scores and native-referenced
# DockQ (quality_dockq) share pose IDs, as after `ppinsight quality`.
FUNNEL_PAIR = ("demo_rec", "demo_lig")
FUNNEL_N_POSES = 150

# Native-referenced DockQ per pair, shaped like real runs: most poses are
# incorrect decoys and a fraction forms a near-native cluster, so every CAPRI
# tier appears.  Pairs get harder down the list.
# (near-native fraction, cluster mean, cluster sd)
_DOCKQ_MIXTURE = [
    {"HADDOCK": (0.36, 0.70, 0.14), "LightDock": (0.26, 0.62, 0.16),
     "Rosetta": (0.30, 0.66, 0.15)},
    {"HADDOCK": (0.22, 0.55, 0.13), "LightDock": (0.16, 0.50, 0.14),
     "Rosetta": (0.24, 0.52, 0.12)},
    {"HADDOCK": (0.12, 0.42, 0.10), "LightDock": (0.14, 0.38, 0.12),
     "Rosetta": (0.10, 0.45, 0.10)},
    {"HADDOCK": (0.04, 0.32, 0.08), "LightDock": (0.02, 0.30, 0.08),
     "Rosetta": (0.06, 0.35, 0.08)},
]

_QUALITY_PARAMS = {
    "quality_fnat": [
        {"HADDOCK": (0.80, 0.10), "LightDock": (0.72, 0.12), "Rosetta": (0.75, 0.10)},
        {"HADDOCK": (0.48, 0.10), "LightDock": (0.40, 0.12), "Rosetta": (0.45, 0.09)},
        {"HADDOCK": (0.22, 0.08), "LightDock": (0.18, 0.09), "Rosetta": (0.25, 0.07)},
        {"HADDOCK": (0.06, 0.04), "LightDock": (0.05, 0.04), "Rosetta": (0.07, 0.04)},
    ],
    "quality_irmsd": [
        {"HADDOCK": (1.2, 0.5), "LightDock": (1.5, 0.7), "Rosetta": (1.3, 0.6)},
        {"HADDOCK": (3.5, 1.0), "LightDock": (4.0, 1.3), "Rosetta": (3.8, 1.1)},
        {"HADDOCK": (6.5, 1.5), "LightDock": (7.5, 2.0), "Rosetta": (6.0, 1.4)},
        {"HADDOCK": (10.0, 2.5), "LightDock": (11.0, 3.0), "Rosetta": (9.5, 2.2)},
    ],
    "quality_lrmsd": [
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


def _mixture_dockq(pair_idx, eng, rng):
    """Decoys plus a near-native cluster, in random pose order."""
    fraction, mu, sigma = _DOCKQ_MIXTURE[pair_idx][eng]
    n_near = round(fraction * N_POSES)
    values = np.concatenate([
        rng.beta(1.3, 9.0, N_POSES - n_near),
        rng.normal(mu, sigma, n_near),
    ])
    return rng.permutation(np.clip(values, 0.0, 1.0))


def _build_quality_scores(engines, pairs, metrics, rng):
    rows = []
    for pair_idx, (protA, protB) in enumerate(pairs):
        for eng in engines:
            for metric in metrics:
                if metric == "quality_dockq":
                    if (
                        pair_idx >= len(_DOCKQ_MIXTURE)
                        or eng not in _DOCKQ_MIXTURE[pair_idx]
                    ):
                        continue
                    values = _mixture_dockq(pair_idx, eng, rng)
                elif metric not in _QUALITY_PARAMS:
                    continue
                else:
                    param_list = _QUALITY_PARAMS[metric]
                    if pair_idx >= len(param_list):
                        continue
                    params = param_list[pair_idx]
                    if eng not in params:
                        continue
                    mu, sigma = params[eng]
                    values = rng.normal(mu, sigma, N_POSES)
                if metric in {"quality_dockq", "quality_fnat"}:
                    values = np.clip(values, 0.0, 1.0)
                elif metric in {"quality_irmsd", "quality_lrmsd"}:
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


# Each engine's own score for the labelled binder / non-binder pairs.  ROC
# needs a score that every pair has, so it cannot use DockQ (which needs a
# native complex that non-binders do not have).
# (score_type, binder mean, non-binder mean, sd)
_LABELLED_ENGINE_SCORES = {
    "HADDOCK": ("score", -60.0, -45.0, 20.0),
    "LightDock": ("luciferin_score", 15.0, 11.0, 6.0),
    "Rosetta": ("i_sc", -12.0, -9.0, 5.0),
}

# Native-referenced DockQ for the binders only (they have a native complex).
_BINDER_DOCKQ = {
    "HADDOCK": (0.45, 0.18), "LightDock": (0.42, 0.20), "Rosetta": (0.48, 0.19),
}


def _build_roc_data(engines, pairs_binder, pairs_nonbinder, rng):
    """Labelled pairs: each engine's own score for every pair, plus
    native-referenced DockQ (quality_dockq) for the binders."""
    rows = []
    for eng in engines:
        if eng not in _LABELLED_ENGINE_SCORES:
            continue
        metric, mu_bind, mu_non, sigma = _LABELLED_ENGINE_SCORES[eng]
        for label, pairs, mu in (
            ("interaction", pairs_binder, mu_bind),
            ("non-interaction", pairs_nonbinder, mu_non),
        ):
            for pose_idx, (protA, protB) in enumerate(pairs, start=1):
                run_id, pose_id, output_path = _traceability_fields(
                    eng, protA, protB, pose_idx,
                )
                values = [(metric, float(rng.normal(mu, sigma)))]
                if label == "interaction":
                    dq_mu, dq_sigma = _BINDER_DOCKQ[eng]
                    dockq = float(np.clip(rng.normal(dq_mu, dq_sigma), 0.0, 1.0))
                    values.append(("quality_dockq", dockq))
                for score_type, value in values:
                    rows.append({
                        "model": eng, "score_type": score_type,
                        "score_value": value,
                        "proteinA": protA, "proteinB": protB,
                        "label": label,
                        "run_id": run_id,
                        "pose_id": pose_id,
                        "output_path": output_path,
                    })
    return pd.DataFrame(rows)


def _build_funnel_scores(engines, rng):
    """One target where each engine's own score partly tracks DockQ."""
    # (native score_type, intercept, slope per DockQ unit, noise sd)
    engine_models = {
        "HADDOCK": ("score", -40.0, -80.0, 12.0),
        "LightDock": ("luciferin_score", 10.0, 12.0, 8.0),
        "Rosetta": ("i_sc", -5.0, -15.0, 4.0),
    }
    protA, protB = FUNNEL_PAIR
    rows = []
    for eng in engines:
        if eng not in engine_models:
            continue
        metric, intercept, slope, noise = engine_models[eng]
        n_native_like = FUNNEL_N_POSES // 6
        dockq = np.concatenate([
            rng.beta(1.2, 10.0, FUNNEL_N_POSES - n_native_like),
            rng.normal(0.62, 0.15, n_native_like),
        ])
        dockq = np.clip(dockq, 0.0, 1.0)
        scores = intercept + slope * dockq + rng.normal(0.0, noise, len(dockq))
        for pose_idx, (dq, score) in enumerate(zip(dockq, scores, strict=False), 1):
            run_id, pose_id, output_path = _traceability_fields(
                eng, protA, protB, pose_idx,
            )
            for score_type, value in ((metric, score), ("quality_dockq", dq)):
                rows.append({
                    "model": eng,
                    "score_type": score_type,
                    "score_value": float(value),
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
            ENGINES, PAIRS,
            ["quality_dockq", "quality_fnat", "quality_irmsd", "quality_lrmsd"],
            rng,
        ),
        _build_engine_scores(ENGINES, PAIRS, rng),
        _build_roc_data(ENGINES, binder_pairs, nonbinder_pairs, rng),
        # Drawn last so the rows above stay identical to earlier versions.
        _build_funnel_scores(ENGINES, rng),
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
    # DockQ for 1, 2 and 3 engines, with 1, 2 and 3 pairs respectively.
    for n in (1, 2, 3):
        df = _subset_scores(
            scores_df,
            engines=ENGINES[:n],
            pairs=PAIRS[:n],
            metrics=["quality_dockq"],
        )
        fig = violin_plot(df, metric="quality_dockq", output=None)
        _save(fig, f"violin_{n}eng_{n}pair_dockq.png")
    # Two engine-native scores, which need no native structure.
    for eng, metric in (("HADDOCK", "score"), ("LightDock", "luciferin_score")):
        df = _subset_scores(
            scores_df, engines=[eng], pairs=PAIRS[:1], metrics=[metric],
        )
        fig = violin_plot(df, metric=metric, output=None)
        _save(fig, f"violin_{eng}_1pair_{metric}.png")


def gen_strip(scores_df):
    print("\n-- Strip plots --")
    for n in (1, 2, 3):
        df = _subset_scores(
            scores_df,
            engines=ENGINES[:n],
            pairs=PAIRS[:n],
            metrics=["quality_dockq"],
        )
        fig = strip_plot(df, metric="quality_dockq", output=None)
        _save(fig, f"strip_{n}eng_{n}pair_dockq.png")


def gen_scatter(scores_df):
    print("\n-- Scatter plots --")
    # The 50 labelled binders have a native complex, so their
    # native-referenced DockQ is meaningful, and 50 pairs give stable
    # correlations.  (Non-binders have no native, hence no DockQ.)
    df = _subset_scores(
        scores_df,
        engines=["HADDOCK", "Rosetta"],
        metrics=["quality_dockq"],
        labels=["interaction"],
    )
    fig = model_agreement_scatter(
        df, metric="quality_dockq", model_x="HADDOCK", model_y="Rosetta",
        output=None,
    )
    _save(fig, "scatter_HADDOCK_vs_Rosetta_dockq.png")


def gen_funnel(scores_df):
    print("\n-- Funnel plots (rank audit) --")
    df = _subset_scores(scores_df, engines=ENGINES, pairs=[FUNNEL_PAIR])
    fig = funnel_plot(df, output=None)
    _save(fig, "funnel_3eng_1pair.png")


def gen_topn(scores_df):
    print("\n-- Top-N plots (rank audit) --")
    df = _subset_scores(scores_df, engines=ENGINES, pairs=[FUNNEL_PAIR])
    fig = topn_plot(df, output=None)
    _save(fig, "topn_3eng_1pair.png")


def gen_roc(scores_df):
    print("\n-- ROC plots --")
    # ROC needs a score every pair has, so it uses HADDOCK's own score;
    # DockQ is not allowed because non-binders have no native structure.
    metric = _LABELLED_ENGINE_SCORES["HADDOCK"][0]
    df = _subset_scores(
        scores_df,
        engines=["HADDOCK"],
        metrics=[metric],
        labels=["interaction", "non-interaction"],
    )
    fig = roc_curve_plot(df, metric=metric, output=None)
    _save(fig, f"roc_HADDOCK_{metric}.png")


def gen_quality_bar(scores_df):
    print("\n-- Quality bar charts --")
    # One protein pair per plot; tiers come from DockQ.
    for n_eng in (1, 2, 3):
        df = _subset_scores(
            scores_df,
            engines=ENGINES[:n_eng],
            pairs=PAIRS[:1],
            metrics=["quality_dockq"],
        )
        fig = quality_bar_chart(df, output=None)
        _save(fig, f"quality_bar_{n_eng}eng_1pair.png")


def gen_cdf(scores_df):
    print("\n-- CDF plots --")
    for n_eng in (1, 2, 3):
        df = _subset_scores(
            scores_df,
            engines=ENGINES[:n_eng],
            pairs=PAIRS[:1],
            metrics=["quality_dockq"],
        )
        fig = cdf_plot(df, metric="quality_dockq", output=None)
        _save(fig, f"cdf_{n_eng}eng_dockq.png")


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
    gen_strip(scores_df)
    gen_scatter(scores_df)
    gen_roc(scores_df)
    gen_quality_bar(scores_df)
    gen_cdf(scores_df)
    gen_funnel(scores_df)
    gen_topn(scores_df)
    n_files = len(list(OUT.glob("*.png")))
    print(f"Fabricated scores rows: {len(scores_df)}")
    print(f"\nDone -- {n_files} plots saved to {OUT.relative_to(ROOT)}/")
    print("\nWARNING: ALL DATA IS FABRICATED. For demonstration only.")


if __name__ == "__main__":
    main()
