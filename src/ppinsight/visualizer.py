"""
visualizer – compare docking scores across PPI prediction models.

Supports three input modes (from simplest to most flexible):

1. **Unified scores file** (``scores.tsv`` / ``scores.csv``)
   Columns: ``proteinA, proteinB, model, score_type, score_value``
   Provenance: companion ``.provenance.json`` sidecar (written by ``collect_scores``)
   → :func:`load_scores` → :func:`compare_scores`

2. **Per-model score files** (one CSV/TSV per model, each with a numeric
   column the user wants to compare)
   → :func:`to_plot` → :func:`compare_scores`

3. **Native tool outputs** read directly
   - HADDOCK ``capri_ss.tsv``  → columns like *score*, *dockq*, *irmsd*, *fnat*, …
    - Rosetta ``docking_scores.csv`` → columns *run*, *total_score*, *i_sc*
   → :func:`to_plot` handles both CSV and TSV automatically.

CLI usage::

    compare_scores scores.tsv --metric dockq
    compare_scores scores.tsv --metric score --pair 2UUY_rec:2UUY_lig
    compare_scores model1.csv model2.csv --metric score --names HADDOCK Rosetta
    compare_scores scores.tsv --metric dockq --output plot.png
"""

import argparse
import os
import sys

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# Theme configuration
# ---------------------------------------------------------------------------

# Named themes are internal; PPInsight uses the default theme.
THEMES: dict[str, dict] = {
    "white-deep":         {"style": "white",      "palette": "deep"},
    "whitegrid-Set2":     {"style": "whitegrid",  "palette": "Set2"},
    "whitegrid-husl":     {"style": "whitegrid",  "palette": "husl"},
    "whitegrid-deep":     {"style": "whitegrid",  "palette": "deep"},
    "whitegrid-muted":    {"style": "whitegrid",  "palette": "muted"},
    "white-pastel":       {"style": "white",      "palette": "pastel"},
    "ticks-colorblind":   {"style": "ticks",      "palette": "colorblind"},
    "darkgrid-Set2":      {"style": "darkgrid",   "palette": "Set2"},
    "white-rocket":       {"style": "white",      "palette": "rocket_r"},
    "ticks-viridis":      {"style": "ticks",      "palette": "viridis"},
}

DEFAULT_THEME = "white-deep"

PLOT_TITLE_SIZE = 14
FACET_TITLE_SIZE = 11
AXIS_LABEL_SIZE = 12
TICK_LABEL_SIZE = 10
LEGEND_FONT_SIZE = 8
LEGEND_TITLE_SIZE = 9
ANNOTATION_FONT_SIZE = 9
N_LABEL_SIZE = 9
TITLE_PAD = 10
QUALITY_BAR_TITLE_PAD = 12
RIDGE_TITLE_OFFSET_IN = 0.29
RIDGE_LEGEND_OFFSET_IN = 0.024

DOCKQ_CAPRI_THRESHOLDS = [
    (0.23, "acceptable", "#f0c040"),
    (0.49, "medium", "#4090e0"),
    (0.80, "high", "#40c040"),
]


def apply_theme(name: str = DEFAULT_THEME) -> None:
    """Apply a named PPInsight theme (seaborn style + palette).

    Parameters
    ----------
    name : str
        Key in :data:`THEMES`.  Falls back to :data:`DEFAULT_THEME` if
        *name* is not recognised.
    """
    cfg = THEMES.get(name, THEMES[DEFAULT_THEME])
    sns.set_theme(style=cfg["style"], palette=cfg["palette"],
                  font_scale=1.0,
                  rc={"figure.dpi": 150, "savefig.dpi": 180,
                      "axes.edgecolor": ".3", "axes.linewidth": 0.8,
                      "axes.titlesize": PLOT_TITLE_SIZE,
                      "axes.labelsize": AXIS_LABEL_SIZE,
                      "axes.titlepad": TITLE_PAD,
                      "xtick.labelsize": TICK_LABEL_SIZE,
                      "ytick.labelsize": TICK_LABEL_SIZE,
                      "legend.fontsize": LEGEND_FONT_SIZE,
                      "legend.title_fontsize": LEGEND_TITLE_SIZE})


def _model_palette(models: list[str]) -> dict[str, tuple[float, float, float]]:
    """Return a stable palette mapping for model names.

    Known engines keep the first three colors of the active seaborn
    palette so HADDOCK / LightDock / Rosetta remain blue / orange /
    green across single-engine and multi-engine plots. Unknown engines
    receive the remaining colors in sorted order.
    """
    unique_models = sorted(set(models))
    base = sns.color_palette(n_colors=max(10, len(unique_models) + 3))
    preferred_positions = {"HADDOCK": 0, "LightDock": 1, "Rosetta": 2}

    palette: dict[str, tuple[float, float, float]] = {}
    for model, idx in preferred_positions.items():
        if model in unique_models:
            palette[model] = base[idx]

    extra_colors = iter(base[3:])
    for model in unique_models:
        if model not in palette:
            palette[model] = next(extra_colors)

    return palette


# ---------------------------------------------------------------------------
# Loading helpers
# ---------------------------------------------------------------------------

def to_plot(models: list[os.PathLike]) -> list[pd.DataFrame]:
    """Read a list of per-model CSV/TSV files into DataFrames.

    Accepts ``.csv`` and ``.tsv`` (tab-separated) files.  The separator is
    inferred from the file extension and falls back to :func:`csv.Sniffer`
    when the extension is ambiguous.

    Raises
    ------
    TypeError
        If *models* is not a list.
    """
    if not isinstance(models, list):
        raise TypeError("Format your input as a list of file paths")

    frames: list[pd.DataFrame] = []
    for m in models:
        sep = _guess_sep(str(m))
        df = pd.read_csv(m, sep=sep, header=0)
        frames.append(df)
    return frames


def load_scores(path: os.PathLike) -> pd.DataFrame:
    """Load a unified scores file (CSV or TSV).

    The file **must** contain at least the columns
    ``model`` and ``score_value``.  ``score_type``, ``proteinA``, and
    ``proteinB`` are optional but enable richer filtering.

    Returns
    -------
    pd.DataFrame
    """
    sep = _guess_sep(str(path))
    df = pd.read_csv(path, sep=sep, header=0)
    # Normalise column names (lowercase, strip whitespace)
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
    required = {"model", "score_value"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Unified scores file is missing required column(s): {missing}. "
            f"Found columns: {list(df.columns)}"
        )
    return df


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def compare_scores(
    frames: list[pd.DataFrame],
    models: list[str],
    score_type: str,
    plot_title: str | None = None,
    output: str | None = None,
) -> plt.Figure:
    """Bar chart comparing *score_type* across *models*.

    Each element of *frames* is a DataFrame that must contain a column named
    *score_type*.  The **mean** of that column is used as the bar height, with
    an error bar showing ±1 standard deviation when there are multiple rows.

    Parameters
    ----------
    frames : list[DataFrame]
        One DataFrame per model (from :func:`to_plot` or manual construction).
    models : list[str]
        Display name for each model (same length as *frames*).
    score_type : str
        Column name to plot (e.g. ``"score"``, ``"dockq"``, ``"irmsd"``).
    plot_title : str, optional
        Figure title.  Defaults to ``"<score_type> by model"``.
    output : str, optional
        If given, save the figure to this path (png/svg/pdf) instead of
        calling ``plt.show()``.

    Returns
    -------
    matplotlib.figure.Figure
        The figure object (useful for programmatic post-processing).

    Raises
    ------
    LookupError
        If *score_type* is not a column in one of the DataFrames.
    """
    if len(frames) != len(models):
        raise ValueError(
            f"Number of frames ({len(frames)}) must match number of "
            f"model names ({len(models)})"
        )

    means: list[float] = []
    stds: list[float] = []
    for name, df in zip(models, frames, strict=False):
        if score_type not in df.columns:
            raise LookupError(
                f"The score type '{score_type}' does not exist for the "
                f"'{name}' model. Available columns: {list(df.columns)}"
            )
        col = pd.to_numeric(df[score_type], errors="coerce").dropna()
        means.append(col.mean())
        stds.append(col.std() if len(col) > 1 else 0.0)

    # Warn if error bars dwarf the bars — bar chart is the wrong
    # visualisation for high-variance data.
    import warnings as _warnings
    for m_val, s_val, name in zip(means, stds, models, strict=False):
        if m_val != 0 and abs(s_val / m_val) > 1.5:
            _warnings.warn(
                f"'{name}' has std ({s_val:.2f}) much larger than the mean "
                f"({m_val:.2f}).  Consider --plot-type violin for the full "
                "distribution or --plot-type ridge for a compact density "
                "comparison.",
                stacklevel=2,
            )

    palette = sns.color_palette(n_colors=len(models))
    bar_width = max(0.45, min(0.8, 0.65 * len(models)))
    fig, ax = plt.subplots(figsize=(max(6, len(models) * 2.5), 5))
    ax.bar(models, means, width=bar_width, yerr=stds, capsize=6,
           color=palette, edgecolor="white", linewidth=0.8,
           error_kw={"linewidth": 1.5})
    ax.set_xlabel("Interaction model")
    ax.set_ylabel(f"{score_type}")
    ax.set_title(plot_title or f"{score_type} by model")
    sns.despine(ax=ax)

    if output:
        fig.savefig(output, bbox_inches="tight", dpi=150)
        print(f"Plot saved to {output}")
    else:
        plt.show()

    return fig


def compare_scores_unified(
    scores_df: pd.DataFrame,
    metric: str,
    pair: tuple[str, str] | None = None,
    plot_title: str | None = None,
    output: str | None = None,
) -> plt.Figure:
    """Bar chart from a **unified** scores DataFrame.

    When the data contains **≥ 2 protein pairs** (and no ``--pair`` filter
    is active), each pair gets its own side-by-side subplot so bars are
    never averaged across different interactions.

    Parameters
    ----------
    scores_df : DataFrame
        Must have columns ``model``, ``score_value`` and optionally
        ``score_type``, ``proteinA``, ``proteinB``.
    metric : str
        The ``score_type`` value to filter on (e.g. ``"dockq"``).
        If the DataFrame has no ``score_type`` column, all rows are used.
    pair : tuple[str, str], optional
        ``(proteinA, proteinB)`` filter.  Ignored when *None*.
    plot_title : str, optional
    output : str, optional
        Save path (png/svg/pdf).

    Returns
    -------
    matplotlib.figure.Figure
    """
    df = scores_df.copy()

    # Filter by score_type if the column exists
    if "score_type" in df.columns:
        df = df[df["score_type"].str.lower() == metric.lower()]

    # Filter by protein pair
    if pair and "proteina" in df.columns and "proteinb" in df.columns:
        a, b = pair
        df = df[
            (df["proteina"].str.lower() == a.lower())
            & (df["proteinb"].str.lower() == b.lower())
        ]

    if df.empty:
        raise ValueError(
            f"No rows match metric='{metric}'"
            + (f" pair={pair}" if pair else "")
        )

    pa, pb, pair_tuples = _detect_pairs(df)
    n_pairs = len(pair_tuples)

    # ── Facet by pair when ≥ 2 pairs and no explicit --pair filter ─
    if n_pairs >= 2 and pair is None:
        fig, axes = plt.subplots(1, n_pairs,
                                 figsize=(max(5, 3 * n_pairs) * 1.2, 5),
                                 sharey=True)
        if n_pairs == 1:
            axes = [axes]
        for idx, (a, b) in enumerate(pair_tuples):
            ax = axes[idx]
            sub = df[(df[pa] == a) & (df[pb] == b)]
            grouped = sub.groupby("model")
            model_names = sorted(grouped.groups.keys())
            means = [pd.to_numeric(grouped.get_group(m)["score_value"],
                                   errors="coerce").dropna().mean()
                     for m in model_names]
            stds = [pd.to_numeric(grouped.get_group(m)["score_value"],
                                  errors="coerce").dropna().std()
                    if len(grouped.get_group(m)) > 1 else 0.0
                    for m in model_names]
            palette = sns.color_palette(n_colors=len(model_names))
            bar_w = max(0.45, min(0.8, 0.65 * len(model_names)))
            ax.bar(model_names, means, width=bar_w, yerr=stds,
                   capsize=6, color=palette, edgecolor="white",
                   linewidth=0.8, error_kw={"linewidth": 1.5})
            ax.set_title(_pair_label(a, b), fontsize=10)
            ax.set_xlabel("")
            ax.set_ylabel(metric if idx == 0 else "")
            sns.despine(ax=ax)

        suptitle = plot_title or f"{metric} by model"
        fig.suptitle(suptitle, fontsize=13, y=1.02)
        fig.tight_layout()
        # Shared x-label centred across all facets
        fig.text(0.5, -0.02, "Docking Engine", ha="center", fontsize=11)

        if output:
            fig.savefig(output, bbox_inches="tight", dpi=150)
            print(f"Plot saved to {output}")
        else:
            plt.show()
        return fig

    # ── Single pair / explicit --pair — single-panel bar chart ────
    grouped = df.groupby("model")
    frames: list[pd.DataFrame] = []
    model_names: list[str] = []
    for name, group in grouped:
        frame = pd.DataFrame({
            metric: pd.to_numeric(
                group["score_value"], errors="coerce"
            )
        })
        frames.append(frame.reset_index(drop=True))
        model_names.append(str(name))

    title = plot_title or (
        f"{metric} by model"
        + (f" ({pair[0]} vs {pair[1]})" if pair else "")
    )
    return compare_scores(frames, model_names, metric, plot_title=title, output=output)


# ---------------------------------------------------------------------------
# Convenience: list available metrics and pairs
# ---------------------------------------------------------------------------

def available_metrics(scores_df: pd.DataFrame) -> list[str]:
    """Return sorted unique ``score_type`` values from a unified scores DataFrame."""
    if "score_type" not in scores_df.columns:
        return sorted(
            c for c in scores_df.columns
            if c not in (
                "model", "proteina", "proteinb", "output_path", "timestamp",
                "run_id", "pose_id", "pose_rank", "label", "family",
                "source", "source_file", "swarm", "cluster_id",
                "cluster_rank", "cluster_size", "cluster_pop",
            )
        )
    return sorted(scores_df["score_type"].dropna().unique())


def available_pairs(scores_df: pd.DataFrame) -> list[tuple[str, str]]:
    """Return unique ``(proteinA, proteinB)`` pairs."""
    if "proteina" not in scores_df.columns or "proteinb" not in scores_df.columns:
        return []
    pairs = scores_df[["proteina", "proteinb"]].drop_duplicates()
    return list(pairs.itertuples(index=False, name=None))


# ---------------------------------------------------------------------------
# Label-aware plotting (interaction vs non-interaction)
# ---------------------------------------------------------------------------

def compare_scores_by_label(
    scores_df: pd.DataFrame,
    metric: str,
    plot_title: str | None = None,
    output: str | None = None,
) -> plt.Figure:
    """Grouped bar chart comparing *metric* across models, split by label.

    Requires a ``label`` column (e.g. ``interaction`` / ``non-interaction``)
    in *scores_df*.  Each model gets two bars side-by-side.

    Parameters
    ----------
    scores_df : DataFrame
        Unified scores with ``model``, ``score_type``, ``score_value``,
        and ``label`` columns.
    metric : str
        The ``score_type`` value to plot.
    plot_title : str, optional
    output : str, optional

    Returns
    -------
    matplotlib.figure.Figure
    """
    df = scores_df.copy()
    if "score_type" in df.columns:
        df = df[df["score_type"].str.lower() == metric.lower()]
    if df.empty:
        raise ValueError(f"No rows match metric='{metric}'")

    if "label" not in df.columns:
        raise ValueError(
            "No 'label' column in scores — cannot split by interaction status. "
            "Use --pairs with collect_scores to annotate, or add a 'label' column."
        )

    df["score_value"] = pd.to_numeric(df["score_value"], errors="coerce")
    models = sorted(df["model"].unique())
    n_models = len(models)

    fig, ax = plt.subplots(figsize=(max(8, n_models * 2), 6))
    sns.barplot(
        data=df, x="model", y="score_value", hue="label",
        order=models, errorbar="sd", capsize=0.08,
        edgecolor="white", linewidth=0.8, ax=ax,
    )

    ax.set_xlabel("Docking Engine")
    ax.set_ylabel(metric)
    ax.set_title(plot_title or f"{metric} — interaction vs non-interaction")
    ax.legend(title="Label", frameon=True, fancybox=True, framealpha=0.9)
    sns.despine(ax=ax)
    fig.tight_layout()

    if output:
        fig.savefig(output, bbox_inches="tight", dpi=150)
        print(f"Plot saved to {output}")
    else:
        plt.show()

    return fig


# ---------------------------------------------------------------------------
# Additional plot types
# ---------------------------------------------------------------------------

def _detect_pairs(df: pd.DataFrame) -> tuple[str, str, list[tuple[str, str]]]:
    """Return (colA, colB, pair_tuples) from a DataFrame.

    Handles both lowercase and camelCase column names.  Returns empty
    list when pair columns are absent.
    """
    if "proteina" in df.columns and "proteinb" in df.columns:
        pa, pb = "proteina", "proteinb"
    elif "proteinA" in df.columns and "proteinB" in df.columns:
        pa, pb = "proteinA", "proteinB"
    else:
        return "", "", []
    tuples = [tuple(r) for r in df[[pa, pb]].drop_duplicates().values.tolist()]
    return pa, pb, tuples


def _pair_label(a: str, b: str) -> str:
    """Human-readable pair label for subplot titles."""
    return f"{a} vs {b}"


def violin_plot(
    scores_df: pd.DataFrame,
    metric: str,
    split_by_label: bool = False,
    plot_title: str | None = None,
    output: str | None = None,
) -> plt.Figure:
    """Violin plot showing the full score distribution per model.

    When the data contains **≥ 2 protein pairs**, each pair gets its own
    side-by-side subplot so distributions are never merged across
    different interactions.

    Uses seaborn for rich rendering: inner quartile lines, optional
    split by interaction label, and automatic pair annotation so the
    user always knows which protein pair(s) are being evaluated.

    Parameters
    ----------
    scores_df : DataFrame
        Unified scores file.
    metric : str
        Which ``score_type`` to plot.
    split_by_label : bool
        If True and a ``label`` column exists, draw separate halves per
        label (interaction vs non-interaction) inside each violin.
    plot_title : str | None
    output : str | None
        Save path.

    Returns
    -------
    matplotlib.figure.Figure
    """
    df = scores_df.copy()
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
    if "score_type" in df.columns:
        df = df[df["score_type"].str.lower() == metric.lower()]
    df["score_value"] = pd.to_numeric(df["score_value"], errors="coerce")
    df = df.dropna(subset=["score_value"])

    if df.empty:
        raise ValueError(f"No rows match metric='{metric}'")

    metric_name = get_metric_display_name(metric)
    models = sorted(df["model"].unique())
    model_palette = _model_palette(models)
    pa, pb, pair_tuples = _detect_pairs(df)
    has_label = "label" in df.columns and split_by_label and df["label"].notna().any()

    # ── Facet by pair when ≥ 2 pairs ──────────────────────────────
    n_pairs = len(pair_tuples)
    if n_pairs >= 2:
        n_models = len(models)
        panel_w = max(5, n_models * 2.0)
        fig, axes = plt.subplots(1, n_pairs,
                                 figsize=(panel_w * n_pairs, 6),
                                 sharey=True)
        if n_pairs == 1:
            axes = [axes]
        for idx, (a, b) in enumerate(pair_tuples):
            ax = axes[idx]
            sub = df[(df[pa] == a) & (df[pb] == b)]
            sub_models = sorted(sub["model"].unique())
            if has_label:
                sns.violinplot(
                    data=sub, x="model", y="score_value", hue="label",
                    split=True, inner="box", cut=0, ax=ax, linewidth=0.8,
                    density_norm="width", order=sub_models,
                )
                if idx == n_pairs - 1:
                    ax.legend(
                        title="Label",
                        fontsize=LEGEND_FONT_SIZE,
                        title_fontsize=LEGEND_TITLE_SIZE,
                    )
                else:
                    ax.legend_.remove() if ax.get_legend() else None
            else:
                sns.violinplot(
                    data=sub, x="model", y="score_value", hue="model",
                    inner="box", cut=0, ax=ax, linewidth=0.8,
                    density_norm="width", order=sub_models,
                    hue_order=sub_models, palette=model_palette,
                    dodge=False, legend=False,
                )
                _draw_deterministic_stripplot(ax=ax, data=sub, order=sub_models)
            ax.set_title(_pair_label(a, b), fontsize=FACET_TITLE_SIZE)
            ax.set_xlabel("")
            ax.set_ylabel(metric_name if idx == 0 else "")
            _add_faint_horizontal_gridlines(ax)
            sns.despine(ax=ax)

        suptitle = plot_title or f"{metric_name} distribution by model"
        fig.suptitle(suptitle, fontsize=PLOT_TITLE_SIZE, y=1.03)
        fig.tight_layout(rect=(0, 0, 1, 0.96), pad=1.0)
        # Shared x-label centred across all facets
        fig.text(0.5, -0.02, "Docking Engine", ha="center", fontsize=AXIS_LABEL_SIZE)

        if output:
            fig.savefig(output, bbox_inches="tight", dpi=180)
            print(f"Plot saved to {output}")
        else:
            plt.show()
        return fig

    # ── Single pair (or no pair columns) — single-panel plot ──────
    pairs = [_pair_label(a, b) for a, b in pair_tuples]
    n_models = len(models)
    fig_w = max(7, n_models * 2.2)
    fig, ax = plt.subplots(figsize=(fig_w, 6))

    if has_label:
        sns.violinplot(
            data=df, x="model", y="score_value", hue="label",
            split=True, inner="box", cut=0, ax=ax, linewidth=0.8,
            density_norm="width", order=models,
        )
        ax.legend(
            title="Label",
            fontsize=LEGEND_FONT_SIZE,
            title_fontsize=LEGEND_TITLE_SIZE,
        )
    else:
        sns.violinplot(
            data=df, x="model", y="score_value", hue="model",
            inner="box", cut=0, ax=ax, linewidth=0.8,
            density_norm="width", order=models,
            hue_order=models, palette=model_palette,
            dodge=False, legend=False,
        )
        _draw_deterministic_stripplot(ax=ax, data=df, order=models)

    ax.set_xlabel("Docking Engine")
    ax.set_ylabel(metric_name)
    _add_faint_horizontal_gridlines(ax)

    if plot_title:
        title = plot_title
    elif len(pairs) == 1:
        title = f"{metric_name} distribution by model — {pairs[0]}"
    else:
        title = f"{metric_name} distribution by model"
    ax.set_title(title, pad=TITLE_PAD)

    sns.despine(ax=ax)
    fig.tight_layout()

    if output:
        fig.savefig(output, bbox_inches="tight", dpi=180)
        print(f"Plot saved to {output}")
    else:
        plt.show()
    return fig


def roc_curve_plot(
    scores_df: pd.DataFrame,
    metric: str,
    higher_is_better: bool | None = None,
    output: str | None = None,
) -> plt.Figure:
    """ROC curve per model — sweep thresholds and plot TPR vs FPR.

    This is the gold-standard way to evaluate a binary classifier.  It
    shows how well each docking engine separates known interacting pairs
    from non-interacting ones at *every possible* threshold, rather than
    at just a single arbitrary cut-off.

    The area under each curve (AUC) is annotated in the legend.  A
    perfect engine has AUC = 1.0; random guessing gives AUC = 0.5.

    Requires a ``label`` column with values ``"interaction"`` and
    ``"non-interaction"``.

    Parameters
    ----------
    scores_df : DataFrame
    metric : str
        Which ``score_type`` to evaluate.
    higher_is_better : bool | None
        If None, looked up from :data:`METRIC_METADATA`.
    output : str | None
    """
    import numpy as np

    df = scores_df.copy()
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
    if "score_type" in df.columns:
        df = df[df["score_type"].str.lower() == metric.lower()]
    df["score_value"] = pd.to_numeric(df["score_value"], errors="coerce")
    df = df.dropna(subset=["score_value"])

    if "label" not in df.columns:
        raise ValueError("Need 'label' column for ROC curve")

    if higher_is_better is None:
        higher_is_better = get_metric_direction(metric)

    # Aggregate to one score per (model, pair)
    agg = df.groupby(
        ["model", "proteina", "proteinb", "label"], as_index=False,
    )["score_value"].mean()
    agg["actual"] = (agg["label"].str.lower() == "interaction").astype(int)

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot([0, 1], [0, 1], "--", color=".6", alpha=0.5, label="Random (AUC=0.50)")

    model_palette = _model_palette(list(agg["model"].unique()))
    metric_name = get_metric_display_name(metric)
    for model, grp in agg.groupby("model"):
        y_true = grp["actual"].values
        y_score = grp["score_value"].values
        if not higher_is_better:
            y_score = -y_score  # flip so higher = better for ROC

        # Simple ROC implementation (no sklearn dependency)
        order = np.argsort(-y_score)
        y_true_sorted = y_true[order]
        tpr_list = [0.0]
        fpr_list = [0.0]
        tp = fp = 0
        n_pos = y_true.sum()
        n_neg = len(y_true) - n_pos
        for yt in y_true_sorted:
            if yt == 1:
                tp += 1
            else:
                fp += 1
            tpr_list.append(tp / n_pos if n_pos else 0)
            fpr_list.append(fp / n_neg if n_neg else 0)

        # AUC via trapezoidal rule.
        # np.trapz was deprecated / removed in NumPy 2.0 in favour of
        # np.trapezoid.  We try the new name first for forward-compat.
        _trapz = getattr(np, "trapezoid", None) or np.trapz
        auc = _trapz(tpr_list, fpr_list)
        ax.plot(fpr_list, tpr_list, label=f"{model} (AUC={auc:.3f})",
                linewidth=2.2, color=model_palette[model])

    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title(f"ROC Curve — {metric_name}", fontsize=PLOT_TITLE_SIZE, pad=TITLE_PAD)
    ax.legend(
        loc="lower right",
        frameon=True,
        fancybox=True,
        framealpha=0.9,
        fontsize=LEGEND_FONT_SIZE,
        title_fontsize=LEGEND_TITLE_SIZE,
    )
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    sns.despine(ax=ax)
    fig.tight_layout(pad=1.0)

    if output:
        fig.savefig(output, bbox_inches="tight", dpi=150)
        print(f"Plot saved to {output}")
    else:
        plt.show()
    return fig


def model_agreement_scatter(
    scores_df: pd.DataFrame,
    metric: str,
    model_x: str,
    model_y: str,
    output: str | None = None,
) -> plt.Figure:
    """Scatter plot: same metric, same pairs, one model per axis.

    Each point is one protein pair.  The x-coordinate is model A's mean
    score for that pair; the y-coordinate is model B's.  If the two
    engines rank pairs the same way, points cluster along the diagonal
    — indicating harmonious scoring.

    Pearson *r* is annotated on the plot so users can immediately
    quantify the agreement.  Points are coloured by label when
    interaction annotations are available.

    Parameters
    ----------
    scores_df : DataFrame
        Unified scores (long format with ``score_type`` / ``score_value``).
    metric : str
        Quality measure present in both models (e.g. ``"dockq"``).
    model_x, model_y : str
        Names of the two models to compare.
    output : str | None
    """
    import numpy as np

    df = scores_df.copy()
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]

    if "score_type" in df.columns:
        df = df[df["score_type"].str.lower() == metric.lower()]
    df["score_value"] = pd.to_numeric(df["score_value"], errors="coerce")
    df = df.dropna(subset=["score_value"])

    # Need pair columns to align rows across models
    if "proteina" not in df.columns or "proteinb" not in df.columns:
        raise ValueError(
            "Need proteinA and proteinB columns to align pairs across "
            "models.  Re-run 'ppinsight collect' with --pair."
        )

    # Aggregate to one score per (model, pair)
    group_cols = ["model", "proteina", "proteinb"]
    extra = ["label"] if "label" in df.columns and df["label"].notna().any() else []
    agg = df.groupby(group_cols + extra, as_index=False)["score_value"].mean()

    dx = agg[agg["model"].str.lower() == model_x.lower()].set_index(
        ["proteina", "proteinb"]
    )
    dy = agg[agg["model"].str.lower() == model_y.lower()].set_index(
        ["proteina", "proteinb"]
    )

    if dx.empty:
        raise ValueError(f"No rows for model '{model_x}' with metric '{metric}'")
    if dy.empty:
        raise ValueError(f"No rows for model '{model_y}' with metric '{metric}'")

    common = dx.index.intersection(dy.index)
    if len(common) == 0:
        raise ValueError(
            f"No common pairs between {model_x} and {model_y} for "
            f"metric '{metric}'.  Both models must have scored the same "
            f"protein pairs."
        )

    x_vals = dx.loc[common, "score_value"].values
    y_vals = dy.loc[common, "score_value"].values
    metric_name = get_metric_display_name(metric)

    fig, ax = plt.subplots(figsize=(7, 7))

    # Colour by label if available
    if "label" in dx.columns:
        labels = dx.loc[common, "label"].values
        plot_df = pd.DataFrame({
            model_x: x_vals, model_y: y_vals, "label": labels,
        })
        sns.scatterplot(
            data=plot_df, x=model_x, y=model_y, hue="label",
            ax=ax, alpha=0.7, s=45, edgecolor="white", linewidth=0.5,
        )
        ax.legend(
            frameon=True,
            fancybox=True,
            framealpha=0.9,
            fontsize=LEGEND_FONT_SIZE,
            title_fontsize=LEGEND_TITLE_SIZE,
        )
    else:
        ax.scatter(
            x_vals, y_vals, alpha=0.7, s=45,
            edgecolors="white", linewidths=0.5,
        )

    # Diagonal reference line (perfect agreement)
    lo = min(x_vals.min(), y_vals.min())
    hi = max(x_vals.max(), y_vals.max())
    margin = (hi - lo) * 0.05 or 0.1
    ax.plot(
        [lo - margin, hi + margin], [lo - margin, hi + margin],
        "--", color=".6", alpha=0.5, zorder=0, label="y = x",
    )

    # Pearson r
    if len(common) >= 2:
        r = np.corrcoef(x_vals, y_vals)[0, 1]
        r_text = f"r = {r:.3f}" if not np.isnan(r) else "r = n/a"
    else:
        r_text = "r = n/a (< 2 points)"
    ax.annotate(
        f"{r_text}  (n = {len(common)} pairs)",
        xy=(0.05, 0.95), xycoords="axes fraction",
        fontsize=ANNOTATION_FONT_SIZE, va="top",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8),
    )

    ax.set_xlabel(f"{model_x} — {metric_name}")
    ax.set_ylabel(f"{model_y} — {metric_name}")
    ax.set_title(
        f"Model agreement: {model_x} vs {model_y} ({metric_name})",
        pad=TITLE_PAD,
    )
    sns.despine(ax=ax)
    fig.tight_layout(pad=1.0)

    if output:
        fig.savefig(output, bbox_inches="tight", dpi=150)
        print(f"Plot saved to {output}")
    else:
        plt.show()
    return fig


# ---------------------------------------------------------------------------
# Score metadata — direction and interpretation per engine / metric
# ---------------------------------------------------------------------------

# This table tells the normalizer which way "good" is for each metric.
# higher_is_better=True means a higher raw score indicates stronger binding.
# higher_is_better=False means a lower raw score is better (like an energy).
#
# When normalizing for cross-engine comparison we *flip* lower-is-better
# metrics so that **after normalization, higher always means "more likely
# to interact"** regardless of the original scoring convention.
METRIC_METADATA: dict[str, dict] = {
    # ── LightDock ──────────────────────────────────────────────────────
    "luciferin_score": {
        "higher_is_better": True,
        "display_name": "Luciferin Score",
        "description": (
            "LightDock luciferin/scoring (DFIRE by default)."
            " Higher = better fit."
        ),
    },
    # ── HADDOCK ────────────────────────────────────────────────────────
    "score": {
        "higher_is_better": False,
        "display_name": "Score",
        "description": "HADDOCK overall score (weighted energy). Lower = better.",
    },
    "dockq": {
        "higher_is_better": True,
        "display_name": "DockQ",
        "bounds": (0.0, 1.0),
        "description": "DockQ quality (0–1). Higher = closer to native complex.",
    },
    "irmsd": {
        "higher_is_better": False,
        "display_name": "iRMSD",
        "description": "Interface RMSD (Å). Lower = better.",
    },
    "lrmsd": {
        "higher_is_better": False,
        "display_name": "lRMSD",
        "description": "Ligand RMSD (Å). Lower = better.",
    },
    "fnat": {
        "higher_is_better": True,
        "display_name": "Fnat",
        "bounds": (0.0, 1.0),
        "description": "Fraction of native contacts recovered. Higher = better.",
    },
    # ── Rosetta ────────────────────────────────────────────────────────
    "interface_score": {
        "higher_is_better": False,
        "display_name": "Interface Score",
        "description": (
            "Legacy Rosetta interface-score label retained for backward "
            "compatibility with older collected files. Lower = better."
        ),
    },
    "i_sc": {
        "higher_is_better": False,
        "display_name": "I_sc",
        "description": (
            "Rosetta interface score (I_sc / dG_separated) — the primary "
            "RosettaDock quality metric.  Lower = stronger binding."
        ),
    },
    "total_score": {
        "higher_is_better": False,
        "display_name": "Total Score",
        "description": "Rosetta total energy (REU). Lower = better.",
    },
    "irms": {
        "higher_is_better": False,
        "display_name": "iRMS",
        "description": "Rosetta interface RMSD (Å). Lower = better.",
    },
    "rms": {
        "higher_is_better": False,
        "display_name": "RMS",
        "description": "Rosetta ligand RMSD (Å). Lower = better.",
    },
    "dg_separated": {
        "higher_is_better": False,
        "display_name": "ΔG Separated",
        "description": (
            "Rosetta dG_separated — binding energy after rigid-body "
            "separation.  Equivalent to I_sc.  Lower = stronger binding."
        ),
    },
    "cluster_size": {
        "higher_is_better": True,
        "display_name": "Cluster Size",
        "description": (
            "Number of decoys in the structural cluster.  Larger clusters "
            "indicate more frequently sampled (more confident) binding modes."
        ),
    },
    # ── DockQ quality metrics (from ppinsight.quality) ─────────────────
    # These use a ``quality_`` prefix to distinguish them from the
    # engine-native metrics above (e.g. HADDOCK's own ``dockq`` score).
    # They are produced by ``add_quality_to_scores()`` which compares
    # predicted complexes against a known native structure via DockQ.
    "quality_dockq": {
        "higher_is_better": True,
        "display_name": "DockQ",
        "bounds": (0.0, 1.0),
        "description": (
            "DockQ score (0–1) from comparison with native structure.  "
            "Higher = better.  Thresholds: acceptable ≥ 0.23, "
            "medium ≥ 0.49, high ≥ 0.80."
        ),
    },
    "quality_fnat": {
        "higher_is_better": True,
        "display_name": "Fnat",
        "bounds": (0.0, 1.0),
        "description": (
            "Fraction of native contacts (Fnat) recovered in the docked "
            "model.  Higher = better."
        ),
    },
    "quality_irmsd": {
        "higher_is_better": False,
        "display_name": "iRMSD",
        "description": (
            "Interface RMSD (Å) between docked model and native structure.  "
            "Lower = better."
        ),
    },
    "quality_lrmsd": {
        "higher_is_better": False,
        "display_name": "lRMSD",
        "description": (
            "Ligand RMSD (Å) between docked model and native structure.  "
            "Lower = better."
        ),
    },
    # ── PRODIGY binding affinity ────────────────────────────────────────
    # Produced by ``ppinsight.prodigy`` which wraps the ``prodigy-prot``
    # package.  Install with ``pip install ppinsight[prodigy]``.
    "prodigy_ddg": {
        "higher_is_better": False,
        "display_name": "PRODIGY ΔG",
        "description": (
            "PRODIGY predicted binding free energy (ΔG, kcal/mol).  "
            "More negative = stronger predicted binding.  "
            "Requires the optional prodigy-prot package."
        ),
    },
    "prodigy_kd": {
        "higher_is_better": False,
        "display_name": "PRODIGY Kd",
        "description": (
            "PRODIGY predicted dissociation constant (Kd, M) at 25 °C.  "
            "Lower = tighter predicted binding.  "
            "Requires the optional prodigy-prot package."
        ),
    },
}


def get_metric_direction(metric: str) -> bool:
    """Return *True* if higher values of *metric* mean better binding.

    Falls back to ``True`` (higher-is-better) for unknown metrics.
    """
    meta = METRIC_METADATA.get(metric.lower())
    if meta is not None:
        return meta["higher_is_better"]
    return True  # safe default


def get_metric_display_name(metric: str) -> str:
    """Return a user-facing display name for *metric*."""
    meta = METRIC_METADATA.get(metric.lower())
    if meta is not None and "display_name" in meta:
        return meta["display_name"]
    return metric.replace("_", " ").title()


def get_metric_bounds(metric: str) -> tuple[float, float] | None:
    """Return the valid numeric bounds for *metric*, if known."""
    meta = METRIC_METADATA.get(metric.lower())
    if meta is None or "bounds" not in meta:
        return None
    lo, hi = meta["bounds"]
    return float(lo), float(hi)


def _draw_deterministic_stripplot(
    ax,
    data: pd.DataFrame,
    order: list[str],
    color: str = ".2",
) -> None:
    """Draw a deterministic strip overlay for violin plots."""
    import numpy as np

    rng_state = np.random.get_state()
    np.random.seed(0)
    try:
        sns.stripplot(
            data=data,
            x="model",
            y="score_value",
            size=1.1,
            alpha=0.07,
            jitter=0.16,
            ax=ax,
            color=color,
            order=order,
            zorder=0,
        )
    finally:
        np.random.set_state(rng_state)


def _add_faint_horizontal_gridlines(ax) -> None:
    """Add subtle y-gridlines without competing with the data."""
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=".82", alpha=0.28, linewidth=0.6)


# ---------------------------------------------------------------------------
# CAPRI four-tier quality classification
# ---------------------------------------------------------------------------
#
# Standard thresholds from the CAPRI assessment protocol.
#
# Reference:
#   Lensink MF, Velankar S, Wodak SJ. (2016). "Prediction of homoprotein
#   and heteroprotein complexes by protein docking and template-based
#   modeling: A CASP-CAPRI experiment." Proteins, 84(Suppl 1):323-348.
#   DOI: 10.1002/prot.25007  (Table 3)
#
# The table defines four quality tiers based on f(nat), L-rms, and I-rms:
#
#   High (★★★):     f(nat) ≥ 0.5  AND (Lrms ≤ 1.0  OR irms ≤ 1.0)
#   Medium (★★):    f(nat) ≥ 0.3  AND (Lrms ≤ 5.0  OR irms ≤ 2.0)
#   Acceptable (★): f(nat) ≥ 0.1  AND (Lrms ≤ 10.0 OR irms ≤ 4.0)
#   Incorrect:      everything else
#
# PPInsight also supports DockQ-only classification when Lrms/irms are
# unavailable, using thresholds from:
#
#   Basu S, Wallner B. (2016). "DockQ: A Quality Measure for
#   Protein-Protein Docking Models." PLoS ONE 11(8):e0161879.
#   DOI: 10.1371/journal.pone.0161879
#
#   High:       DockQ ≥ 0.80
#   Medium:     DockQ ≥ 0.49
#   Acceptable: DockQ ≥ 0.23
#   Incorrect:  DockQ < 0.23

CAPRI_QUALITY_TIERS = ("high", "medium", "acceptable", "incorrect")


def _capri_pose_key_columns(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Return a copy of *df* plus the best available pose-level key columns.

    Preference order:
    1. explicit ``pose_id``
    2. explicit ``output_path``
    3. reconstructed ordinal within each run/pair/metric group
    """
    pose_cols = ["model"]
    for col in ("proteina", "proteinb", "run_id"):
        if col in df.columns:
            pose_cols.append(col)

    if (
        "pose_id" in df.columns
        and df["pose_id"].fillna("").astype(str).str.strip().ne("").any()
    ):
        return df, pose_cols + ["pose_id"]

    if (
        "output_path" in df.columns
        and df["output_path"].fillna("").astype(str).str.strip().ne("").any()
    ):
        return df, pose_cols + ["output_path"]

    reconstructed = df.copy()
    ordinal_group_cols = [*pose_cols, "score_type"]
    reconstructed["_pose_ordinal"] = reconstructed.groupby(
        ordinal_group_cols
    ).cumcount()
    return reconstructed, pose_cols + ["_pose_ordinal"]


def _capri_pose_table(scores_df: pd.DataFrame) -> pd.DataFrame:
    """Build one CAPRI classification row per pose.

    Uses explicit pose metadata when present and falls back to row-order
    reconstruction for older unified score files that lack ``pose_id``.
    """
    df = scores_df.copy()
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
    df["score_value"] = pd.to_numeric(df["score_value"], errors="coerce")

    if "score_type" not in df.columns:
        return pd.DataFrame(columns=["model", "capri_quality"])

    capri_metrics = {
        "fnat", "irmsd", "lrmsd", "dockq",
        "quality_fnat", "quality_irmsd", "quality_lrmsd", "quality_dockq",
    }
    df = df[df["score_type"].str.lower().isin(capri_metrics)].dropna(
        subset=["score_value"]
    )
    if df.empty:
        return pd.DataFrame(columns=["model", "capri_quality"])

    df, pose_cols = _capri_pose_key_columns(df)
    pivot = df.pivot_table(
        index=pose_cols,
        columns="score_type",
        values="score_value",
        aggfunc="first",
    )
    pivot.columns = [c.lower() for c in pivot.columns]

    def _row_quality(row):
        fnat = row.get("fnat") if "fnat" in row.index else row.get("quality_fnat")
        lrmsd = row.get("lrmsd") if "lrmsd" in row.index else row.get("quality_lrmsd")
        irmsd = row.get("irmsd") if "irmsd" in row.index else row.get("quality_irmsd")
        dockq = row.get("dockq") if "dockq" in row.index else row.get("quality_dockq")

        import math

        def _nan_to_none(val):
            if val is not None and isinstance(val, float) and math.isnan(val):
                return None
            return val

        return capri_quality(
            fnat=_nan_to_none(fnat),
            lrmsd=_nan_to_none(lrmsd),
            irmsd=_nan_to_none(irmsd),
            dockq=_nan_to_none(dockq),
        )

    pose_table = pivot.reset_index()
    pose_table["capri_quality"] = pivot.apply(_row_quality, axis=1).values
    return pose_table


def capri_quality(
    fnat: float | None = None,
    lrmsd: float | None = None,
    irmsd: float | None = None,
    dockq: float | None = None,
) -> str:
    """Classify a single docking model into a CAPRI quality tier.

    Uses the three-metric CAPRI protocol (f(nat) + Lrms + irms) when all
    three values are available.  Falls back to DockQ-only classification
    when RMSD values are missing.

    Parameters
    ----------
    fnat : float or None
        Fraction of native contacts (0–1).
    lrmsd : float or None
        Ligand RMSD in Å.
    irmsd : float or None
        Interface RMSD in Å.
    dockq : float or None
        DockQ score (0–1).  Used as fallback when RMSD values are None.

    Returns
    -------
    str
        One of ``"high"``, ``"medium"``, ``"acceptable"``, ``"incorrect"``.
    """
    # Full CAPRI protocol: fnat + at least one RMSD metric
    if fnat is not None and (lrmsd is not None or irmsd is not None):
        _lrmsd = lrmsd if lrmsd is not None else float("inf")
        _irmsd = irmsd if irmsd is not None else float("inf")

        if fnat >= 0.5 and (_lrmsd <= 1.0 or _irmsd <= 1.0):
            return "high"
        if fnat >= 0.3 and (_lrmsd <= 5.0 or _irmsd <= 2.0):
            return "medium"
        if fnat >= 0.1 and (_lrmsd <= 10.0 or _irmsd <= 4.0):
            return "acceptable"
        return "incorrect"

    # DockQ-only fallback
    if dockq is not None:
        if dockq >= 0.80:
            return "high"
        if dockq >= 0.49:
            return "medium"
        if dockq >= 0.23:
            return "acceptable"
        return "incorrect"

    return "incorrect"


def classify_capri(
    scores_df: pd.DataFrame,
) -> pd.DataFrame:
    """Add a ``capri_quality`` column to a unified scores DataFrame.

    The function pivots the long-format scores so that each
    (model, proteinA, proteinB) group has access to fnat, lrmsd, irmsd,
    and/or dockq values, then applies :func:`capri_quality` row-wise.

    Parameters
    ----------
    scores_df : DataFrame
        Unified scores (long format).

    Returns
    -------
    DataFrame
        A copy with an added ``capri_quality`` column containing one of
        ``"high"``, ``"medium"``, ``"acceptable"``, ``"incorrect"`` per
        row based on the model+pair combination.
    """
    df = scores_df.copy()
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
    df["score_value"] = pd.to_numeric(df["score_value"], errors="coerce")

    if "score_type" not in df.columns:
        df["capri_quality"] = "incorrect"
        return df

    pose_table = _capri_pose_table(df)
    if pose_table.empty:
        df["capri_quality"] = "incorrect"
        return df

    join_cols = ["model"]
    for col in ("proteina", "proteinb", "run_id"):
        if col in df.columns and col in pose_table.columns:
            join_cols.append(col)
    if "pose_id" in df.columns and "pose_id" in pose_table.columns:
        join_cols.append("pose_id")
    elif "output_path" in df.columns and "output_path" in pose_table.columns:
        join_cols.append("output_path")
    elif "_pose_ordinal" in pose_table.columns:
        df["_pose_ordinal"] = df.groupby([*join_cols, "score_type"]).cumcount()
        join_cols.append("_pose_ordinal")

    df = df.merge(
        pose_table[join_cols + ["capri_quality"]],
        on=join_cols,
        how="left",
    )
    if "_pose_ordinal" in df.columns:
        df = df.drop(columns=["_pose_ordinal"])
    df["capri_quality"] = df["capri_quality"].fillna("incorrect")
    return df


def capri_summary_table(
    scores_df: pd.DataFrame,
) -> pd.DataFrame:
    """Summary table: count of models in each CAPRI quality tier per engine.

    Analogous to the ranked summary tables in CAPRI assessments that
    report how many acceptable / medium / high models each group produced.

    Parameters
    ----------
    scores_df : DataFrame
        Unified scores (long format).  Must have metrics that
        :func:`classify_capri` can use (fnat, irmsd, lrmsd, or dockq).

    Returns
    -------
    DataFrame
        Columns: ``model``, ``high``, ``medium``, ``acceptable``,
        ``incorrect``, ``total``, ``pct_acceptable_plus``.
    """
    pose_table = _capri_pose_table(scores_df)
    if pose_table.empty or "capri_quality" not in pose_table.columns:
        raise ValueError(
            "capri_summary_table requires CAPRI metrics (fnat, irmsd, lrmsd, "
            "or dockq) but none were found.  Pass a scores DataFrame that "
            "includes at least one of: fnat, irmsd, lrmsd, dockq."
        )

    rows: list[dict] = []
    for model, grp in pose_table.groupby("model"):
        counts = grp["capri_quality"].value_counts()
        h = int(counts.get("high", 0))
        m = int(counts.get("medium", 0))
        a = int(counts.get("acceptable", 0))
        i = int(counts.get("incorrect", 0))
        total = h + m + a + i
        pct = round((h + m + a) / total * 100, 1) if total else 0.0
        rows.append({
            "model": model,
            "high": h,
            "medium": m,
            "acceptable": a,
            "incorrect": i,
            "total": total,
            "pct_acceptable+": pct,
        })

    result = pd.DataFrame(rows)
    # Sort by pct_acceptable+ descending, then by high count
    result = result.sort_values(
        ["pct_acceptable+", "high", "medium"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    return result


def quality_bar_chart(
    scores_df: pd.DataFrame,
    plot_title: str | None = None,
    output: str | None = None,
) -> plt.Figure:
    """100 % stacked bar chart: CAPRI quality tier distribution per engine.

    Each engine gets one bar that fills from 0 % to 100 %.  The bar is
    divided into four colour-coded tiers (incorrect, acceptable, medium,
    high).  A **N = X poses** label is shown to the right of each bar so
    the reader knows how many individual poses the percentages are based on.

    This function requires CAPRI-compatible metrics (fnat, irmsd, lrmsd,
    and/or dockq).  A :exc:`ValueError` is raised if none are present.

    Inspired by the cumulative bar charts in CAPRI assessments (Lensink
    MF *et al.*, *Proteins* 84(S1):323-348, 2016).

    Parameters
    ----------
    scores_df : DataFrame
        Unified scores with CAPRI-relevant metrics (fnat, irmsd, lrmsd,
        and/or dockq).  Engine-only metrics (e.g. HADDOCK ``score``) are
        ignored; only pose-level quality metrics drive the classification.
    plot_title : str or None
    output : str or None

    Returns
    -------
    matplotlib.figure.Figure

    Raises
    ------
    ValueError
        If *scores_df* contains no CAPRI metrics (fnat / irmsd / lrmsd /
        dockq).
    """
    import numpy as np

    # ── CAPRI guard ────────────────────────────────────────────────────
    _CAPRI_METRICS = {"fnat", "irmsd", "lrmsd", "dockq",
                      "quality_fnat", "quality_irmsd",
                      "quality_lrmsd", "quality_dockq"}
    present: set[str] = set()
    if "score_type" in scores_df.columns:
        present = _CAPRI_METRICS & set(
            scores_df["score_type"].str.lower().unique()
        )
    if not present:
        raise ValueError(
            "quality_bar_chart requires CAPRI metrics (fnat, irmsd, lrmsd, "
            "or dockq) but none were found in 'score_type'.  "
            "Pass a scores DataFrame that includes at least one of: "
            + ", ".join(sorted(_CAPRI_METRICS))
            + "."
        )

    # ── Classify per pose ──────────────────────────────────────────────
    classified = _capri_pose_table(scores_df)
    if classified.empty or "capri_quality" not in classified.columns:
        raise ValueError("CAPRI classification failed — no quality metrics found")

    models = sorted(classified["model"].unique())
    n = len(models)

    tier_colours = {
        "incorrect":  "#b0b0b0",
        "acceptable": "#31a354",
        "medium":     "#3182bd",
        "high":       "#6a51a3",
    }
    tier_order = ["incorrect", "acceptable", "medium", "high"]

    # Width: extra right margin for N= labels
    fig_w = max(6, n * 2.2) + 1.5
    fig, ax = plt.subplots(figsize=(fig_w, 5))
    x = np.arange(n)
    width = 0.55

    # Pre-compute percentages and totals per model
    pct_table: dict[str, dict[str, float]] = {}
    total_table: dict[str, int] = {}
    for model in models:
        sub = classified[classified["model"] == model]["capri_quality"]
        total = len(sub)
        total_table[model] = total
        counts = sub.value_counts()
        pct_table[model] = {
            tier: float(counts.get(tier, 0)) / total * 100 if total else 0.0
            for tier in tier_order
        }

    # Draw 100 % stacked bars
    bottom = np.zeros(n)
    handles = []
    for tier in tier_order:
        values = np.array([pct_table[m][tier] for m in models])
        bars = ax.bar(
            x, values, width, bottom=bottom,
            label=tier.capitalize(),
            color=tier_colours[tier],
            edgecolor="white", linewidth=0.8,
        )
        handles.append(bars[0])
        # Annotate segment if wide enough (≥ 4 %)
        for i, v in enumerate(values):
            if v >= 4.0:
                ax.text(
                    x[i], bottom[i] + v / 2,
                    f"{v:.0f}%",
                    ha="center", va="center",
                    fontsize=8, fontweight="bold", color="white",
                )
        bottom += values

    # N = X label to the right of each bar
    for i, model in enumerate(models):
        ax.text(
            x[i] + width / 2 + 0.07,
            100,
            f"N={total_table[model]}",
            ha="left", va="top",
            fontsize=N_LABEL_SIZE, color="#444444",
        )

    ax.set_xticks(x)
    ax.set_xticklabels(
        models,
        rotation=30 if n > 4 else 0,
        ha="right" if n > 4 else "center",
    )
    ax.set_ylim(0, 100)
    ax.set_xlabel("Docking Engine", labelpad=8)
    ax.set_ylabel("% of poses")
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0f}%"))
    ax.set_title(
        plot_title
        or "Pose-Level DockQ / CAPRI Quality Tier Distribution by Engine",
        pad=QUALITY_BAR_TITLE_PAD,
    )

    # Legend below x-axis title (labelpad above keeps space clear)
    ax.legend(
        handles=handles,
        labels=[t.capitalize() for t in tier_order],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.205),
        ncol=4,
        frameon=True, fancybox=True, framealpha=0.9,
        title="CAPRI tier",
        fontsize=LEGEND_FONT_SIZE,
        title_fontsize=LEGEND_TITLE_SIZE,
    )

    sns.despine(ax=ax)
    fig.tight_layout(rect=(0, 0.11, 1, 1), pad=1.0)

    if output:
        fig.savefig(output, bbox_inches="tight", dpi=150)
        print(f"Plot saved to {output}")
    else:
        plt.show()
    return fig


# ---------------------------------------------------------------------------
# Ridge plot
# ---------------------------------------------------------------------------

def ridge_plot(
    scores_df: pd.DataFrame,
    metric: str,
    plot_title: str | None = None,
    output: str | None = None,
) -> plt.Figure:
    """Ridgeline (joy) plot: one KDE density row per engine.

    Each engine occupies its own horizontally aligned KDE panel.  Panels
    share the same x-axis so distributions from different engines can be
    compared by eye.  This plot is especially useful when you have many
    engines and want to compare their score distributions without the
    visual clutter of overlapping violins.

    Parameters
    ----------
    scores_df : DataFrame
        Unified scores (long format).
    metric : str
        ``score_type`` value to visualise.
    plot_title : str or None
    output : str or None

    Returns
    -------
    matplotlib.figure.Figure
    """
    import numpy as np
    from scipy.stats import gaussian_kde

    df = scores_df.copy()
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
    sub = df[df["score_type"].str.lower() == metric.lower()].copy()
    sub["score_value"] = pd.to_numeric(sub["score_value"], errors="coerce")
    sub = sub.dropna(subset=["score_value"])
    if sub.empty:
        raise ValueError(f"No data found for metric '{metric}'")

    models = sorted(sub["model"].unique())
    n = len(models)
    metric_name = get_metric_display_name(metric)
    metric_bounds = get_metric_bounds(metric)
    model_palette = _model_palette(models)
    palette = [model_palette[m] for m in models]

    fig_h = max(3, n * 1.6)
    fig, axes = plt.subplots(n, 1, figsize=(8, fig_h), sharex=True)
    if n == 1:
        axes = [axes]

    all_vals = sub["score_value"].values
    if metric_bounds is not None:
        x_min, x_max = metric_bounds
    else:
        x_min, x_max = float(all_vals.min()), float(all_vals.max())
        pad = (x_max - x_min) * 0.05 or 0.5
        x_min -= pad
        x_max += pad
    xs = np.linspace(x_min, x_max, 400)
    threshold_guides = DOCKQ_CAPRI_THRESHOLDS if metric.lower() == "dockq" else []

    for ax, model, color in zip(axes, models, palette, strict=False):
        vals = sub[sub["model"] == model]["score_value"].values
        if len(vals) >= 2:
            ys = gaussian_kde(vals)(xs)
        else:
            ys = np.zeros_like(xs)
        for threshold, _, guide_color in threshold_guides:
            ax.axvline(
                threshold,
                color=guide_color,
                lw=0.9,
                ls="--",
                alpha=0.22,
                zorder=0,
            )
        ax.fill_between(xs, ys, alpha=0.55, color=color)
        ax.plot(xs, ys, color=color, lw=1.5)
        median_val = float(np.median(vals))
        median_y = float(np.interp(median_val, xs, ys)) if ys.size else 0.0
        ax.plot(
            [median_val],
            [median_y],
            color=color,
            marker="o",
            linestyle="None",
            markersize=5.5,
            markerfacecolor="white",
            markeredgecolor=color,
            markeredgewidth=1.3,
            zorder=3,
        )
        ax.set_ylabel(model, rotation=0, ha="right", va="center",
                      labelpad=6, fontsize=TICK_LABEL_SIZE)
        ax.set_yticks([])
        sns.despine(ax=ax, left=True, bottom=False)
        ax.tick_params(axis="x", which="both",
                       labelbottom=(ax is axes[-1]))

    if metric_bounds is not None:
        axes[-1].set_xlim(metric_bounds)
    axes[-1].set_xlabel(metric_name)
    legend_handles = [
        Line2D(
            [0], [0],
            marker="o",
            linestyle="None",
            markersize=5.5,
            markerfacecolor="white",
            markeredgecolor="#444444",
            markeredgewidth=1.2,
            label="Median",
        )
    ]
    for threshold, label, color in threshold_guides:
        legend_handles.append(
            Line2D(
                [0], [0],
                color=color,
                lw=1.2,
                ls="--",
                alpha=0.7,
                label=f"CAPRI {label} ({threshold:.2f})",
            )
        )
    fig.suptitle(
        plot_title or f"Ridge plot – {metric_name}",
        fontsize=PLOT_TITLE_SIZE,
        y=1 + (RIDGE_TITLE_OFFSET_IN / fig_h),
    )
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 1 + (RIDGE_LEGEND_OFFSET_IN / fig_h)),
        ncol=min(4, len(legend_handles)),
        frameon=True,
        fancybox=True,
        framealpha=0.9,
        fontsize=LEGEND_FONT_SIZE,
        handlelength=1.8,
        columnspacing=1.0,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90), pad=1.0)

    if output:
        fig.savefig(output, bbox_inches="tight", dpi=150)
        print(f"Plot saved to {output}")
    else:
        plt.show()
    return fig


# ---------------------------------------------------------------------------
# CDF plot
# ---------------------------------------------------------------------------

def cdf_plot(
    scores_df: pd.DataFrame,
    metric: str,
    plot_title: str | None = None,
    output: str | None = None,
) -> plt.Figure:
    """Empirical CDF plot: one step-function curve per engine.

    The CDF shows the fraction of poses scoring ≤ x.  For quality
    metrics where *higher is better* (e.g. DockQ, fnat) the CDF rises
    steeply at higher values for better engines.  CAPRI tier thresholds
    are drawn as vertical dashed lines when *metric* is ``dockq``.

    Parameters
    ----------
    scores_df : DataFrame
        Unified scores (long format).
    metric : str
        ``score_type`` value to visualise.
    plot_title : str or None
    output : str or None

    Returns
    -------
    matplotlib.figure.Figure
    """
    import numpy as np

    df = scores_df.copy()
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
    sub = df[df["score_type"].str.lower() == metric.lower()].copy()
    sub["score_value"] = pd.to_numeric(sub["score_value"], errors="coerce")
    sub = sub.dropna(subset=["score_value"])
    if sub.empty:
        raise ValueError(f"No data found for metric '{metric}'")

    models = sorted(sub["model"].unique())
    metric_name = get_metric_display_name(metric)
    metric_bounds = get_metric_bounds(metric)
    model_palette = _model_palette(models)

    fig, ax = plt.subplots(figsize=(7, 4.5))

    for model in models:
        color = model_palette[model]
        vals = np.sort(sub[sub["model"] == model]["score_value"].values)
        cdf = np.arange(1, len(vals) + 1) / len(vals)
        ax.step(vals, cdf, where="post", color=color, lw=2, label=model)

    # CAPRI threshold lines for dockq
    if metric.lower() == "dockq":
        for thr, label, color in DOCKQ_CAPRI_THRESHOLDS:
            ax.axvline(thr, color=color, lw=1.2, ls="--", alpha=0.8,
                       label=f"CAPRI {label} ({thr})")

    if metric_bounds is not None:
        ax.set_xlim(metric_bounds)
    ax.set_xlabel(metric_name)
    ax.set_ylabel("Cumulative fraction of poses")
    ax.set_ylim(0, 1.05)
    ax.set_title(plot_title or f"Empirical CDF – {metric_name}", pad=TITLE_PAD)
    ax.legend(
        frameon=True,
        fancybox=True,
        framealpha=0.9,
        fontsize=LEGEND_FONT_SIZE,
        title_fontsize=LEGEND_TITLE_SIZE,
    )
    sns.despine(ax=ax)
    fig.tight_layout(pad=1.0)

    if output:
        fig.savefig(output, bbox_inches="tight", dpi=150)
        print(f"Plot saved to {output}")
    else:
        plt.show()
    return fig


# ---------------------------------------------------------------------------
# Pairwise difference plot
# ---------------------------------------------------------------------------

def pairwise_difference_plot(
    scores_df: pd.DataFrame,
    metric: str,
    model_a: str,
    model_b: str,
    plot_title: str | None = None,
    output: str | None = None,
) -> plt.Figure:
    """Histogram of per-pair mean score differences (model A − model B).

    For each protein pair, this plot computes the mean score of *model_a*
    and the mean score of *model_b*, then histograms their differences.
    Bars to the right of zero indicate pairs where *model_a* scored
    higher; bars to the left indicate *model_b* scored higher.

    A vertical dashed line at zero is drawn for reference.  The median
    difference is annotated on the plot.

    Parameters
    ----------
    scores_df : DataFrame
        Unified scores (long format).  Must contain both *model_a* and
        *model_b* rows.
    metric : str
        ``score_type`` value to compute differences for.
    model_a : str
        First model name (positive axis).
    model_b : str
        Second model name (negative axis).
    plot_title : str or None
    output : str or None

    Returns
    -------
    matplotlib.figure.Figure

    Raises
    ------
    ValueError
        If either model is not found or there are no shared pairs.
    """
    import numpy as np

    df = scores_df.copy()
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
    sub = df[df["score_type"].str.lower() == metric.lower()].copy()
    sub["score_value"] = pd.to_numeric(sub["score_value"], errors="coerce")
    sub = sub.dropna(subset=["score_value"])
    metric_name = get_metric_display_name(metric)

    for name in (model_a, model_b):
        if name not in sub["model"].values:
            raise ValueError(
                f"Model '{name}' not found in scores for metric '{metric}'.  "
                f"Available: {sorted(sub['model'].unique())}"
            )

    group_cols = ["model"]
    if "proteina" in sub.columns and "proteinb" in sub.columns:
        group_cols = ["model", "proteina", "proteinb"]

    means = sub.groupby(group_cols)["score_value"].mean().reset_index()
    means.columns = [*group_cols, "mean_score"]

    pair_cols = [c for c in group_cols if c != "model"]

    if pair_cols:
        a_idx = means[means["model"] == model_a].set_index(pair_cols)["mean_score"]
        b_idx = means[means["model"] == model_b].set_index(pair_cols)["mean_score"]
    else:
        a_idx = means[means["model"] == model_a]["mean_score"].reset_index(drop=True)
        b_idx = means[means["model"] == model_b]["mean_score"].reset_index(drop=True)

    diff = (a_idx - b_idx).dropna()
    if diff.empty:
        raise ValueError(
            f"No shared protein pairs found between '{model_a}' and '{model_b}'."
        )

    median_diff = float(np.median(diff.values))

    fig, ax = plt.subplots(figsize=(7, 4.5))
    color = _model_palette([model_a, model_b])[model_a]
    ax.hist(diff.values, bins="auto", color=color, edgecolor="white",
            linewidth=0.6, alpha=0.85)
    ax.axvline(0, color="#444444", lw=1.5, ls="--", label="No difference")
    ax.axvline(median_diff, color="#e05050", lw=1.5, ls="-",
               label=f"Median Δ = {median_diff:+.3f}")
    ax.set_xlabel(f"{metric_name}  ({model_a} − {model_b})")
    ax.set_ylabel("Number of pairs")
    ax.set_title(
        plot_title or f"Per-pair {metric_name} difference: {model_a} vs {model_b}",
        pad=TITLE_PAD,
    )
    ax.legend(
        frameon=True,
        fancybox=True,
        framealpha=0.9,
        fontsize=LEGEND_FONT_SIZE,
        title_fontsize=LEGEND_TITLE_SIZE,
    )
    sns.despine(ax=ax)
    fig.tight_layout(pad=1.0)

    if output:
        fig.savefig(output, bbox_inches="tight", dpi=150)
        print(f"Plot saved to {output}")
    else:
        plt.show()
    return fig


def ranking_table(
    scores_df: pd.DataFrame,
    metric: str,
) -> pd.DataFrame:
    """Per-pair ranking table: which engine scores best on each pair.

    For each protein pair and the given metric, ranks engines from best
    to worst (respecting score direction from :data:`METRIC_METADATA`).
    Returns a tidy DataFrame with columns: ``pair``, ``rank``, ``model``,
    ``mean_score``.

    Parameters
    ----------
    scores_df : DataFrame
        Unified scores (long format).
    metric : str
        Which ``score_type`` to rank by.

    Returns
    -------
    DataFrame
        Columns: ``pair``, ``rank``, ``model``, ``mean_score``.
    """
    df = scores_df.copy()
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]

    if "score_type" in df.columns:
        df = df[df["score_type"].str.lower() == metric.lower()]
    df["score_value"] = pd.to_numeric(df["score_value"], errors="coerce")
    df = df.dropna(subset=["score_value"])

    if df.empty:
        raise ValueError(f"No rows match metric='{metric}'")

    higher_is_better = get_metric_direction(metric)

    group_cols = ["model"]
    if "proteina" in df.columns and "proteinb" in df.columns:
        group_cols = ["proteina", "proteinb", "model"]

    agg = df.groupby(group_cols, as_index=False)["score_value"].mean()

    if "proteina" in agg.columns and "proteinb" in agg.columns:
        agg["pair"] = agg["proteina"] + " vs " + agg["proteinb"]
    else:
        agg["pair"] = "all"

    rows: list[dict] = []
    for pair_name, grp in agg.groupby("pair"):
        sorted_grp = grp.sort_values(
            "score_value", ascending=not higher_is_better,
        ).reset_index(drop=True)
        for rank_idx, (_, row) in enumerate(sorted_grp.iterrows(), start=1):
            rows.append({
                "pair": pair_name,
                "rank": rank_idx,
                "model": row["model"],
                "mean_score": round(row["score_value"], 4),
            })

    return pd.DataFrame(rows)


def tabular_summary(
    scores_df: pd.DataFrame,
    metric: str,
) -> pd.DataFrame:
    """Aggregated summary table: per-model, per-pair statistics.

    Produces a table with mean, std, median, min, max, and count for
    the given metric, grouped by model and pair.  This is the tabular
    complement to the graphical plots — every CAPRI assessment paper
    leads with summary tables before showing figures.

    Parameters
    ----------
    scores_df : DataFrame
        Unified scores (long format).
    metric : str
        Which ``score_type`` to summarise.

    Returns
    -------
    DataFrame
        Columns: ``model``, ``pair``, ``count``, ``mean``, ``std``,
        ``median``, ``min``, ``max``.
    """
    df = scores_df.copy()
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]

    if "score_type" in df.columns:
        df = df[df["score_type"].str.lower() == metric.lower()]
    df["score_value"] = pd.to_numeric(df["score_value"], errors="coerce")
    df = df.dropna(subset=["score_value"])

    if df.empty:
        raise ValueError(f"No rows match metric='{metric}'")

    if "proteina" in df.columns and "proteinb" in df.columns:
        df["pair"] = df["proteina"] + " vs " + df["proteinb"]
        group_cols = ["model", "pair"]
    else:
        df["pair"] = "all"
        group_cols = ["model", "pair"]

    agg = df.groupby(group_cols)["score_value"].agg(
        ["count", "mean", "std", "median", "min", "max"]
    ).reset_index()

    # Round numeric columns
    for col in ["mean", "std", "median", "min", "max"]:
        agg[col] = agg[col].round(4)
    agg["count"] = agg["count"].astype(int)

    return agg.sort_values(["pair", "model"]).reset_index(drop=True)


def metric_overview_table(scores_df: pd.DataFrame) -> pd.DataFrame:
    """Metric-agnostic overview table grouped by model and pair.

    This summary does not require a metric filter. It is meant as a quick
    inventory of what data is present before choosing plot/metric options.

    Returns
    -------
    DataFrame
        Columns: ``model``, ``pair``, ``rows``, ``poses``, ``metrics``,
        ``metric_count``, ``labels_present``, ``paths_present``.
    """
    df = scores_df.copy()
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]

    if "model" not in df.columns:
        raise ValueError("overview table requires a 'model' column")

    if "proteina" in df.columns and "proteinb" in df.columns:
        df["pair"] = (
            df["proteina"].fillna("").astype(str)
            + " vs "
            + df["proteinb"].fillna("").astype(str)
        )
        df.loc[df["pair"].str.strip() == "vs", "pair"] = "all"
    else:
        df["pair"] = "all"

    def _joined_metrics(values: pd.Series) -> str:
        tokens = sorted({str(v).strip() for v in values if str(v).strip()})
        return ";".join(tokens)

    group_cols = ["model", "pair"]
    overview = df.groupby(group_cols, dropna=False).size().reset_index(name="rows")

    if "pose_id" in df.columns:
        pose_counts = (
            df.assign(_pose=df["pose_id"].fillna("").astype(str).str.strip())
            .groupby(group_cols)["_pose"]
            .apply(lambda s: s[s != ""].nunique())
            .reset_index(name="poses")
        )
        overview = overview.merge(pose_counts, on=group_cols, how="left")
        overview["poses"] = overview["poses"].fillna(0).astype(int)
    else:
        overview["poses"] = overview["rows"].astype(int)

    if "score_type" in df.columns:
        metrics = (
            df.groupby(group_cols)["score_type"]
            .agg(metrics=_joined_metrics, metric_count="nunique")
            .reset_index()
        )
        overview = overview.merge(metrics, on=group_cols, how="left")
    else:
        overview["metrics"] = ""
        overview["metric_count"] = 0

    if "label" in df.columns:
        labels_present = (
            df.assign(_label=df["label"].fillna("").astype(str).str.strip())
            .groupby(group_cols)["_label"]
            .apply(lambda s: any(v and v.lower() != "unknown" for v in s))
            .reset_index(name="labels_present")
        )
        overview = overview.merge(labels_present, on=group_cols, how="left")
    else:
        overview["labels_present"] = False

    if "output_path" in df.columns:
        paths_present = (
            df.assign(_path=df["output_path"].fillna("").astype(str).str.strip())
            .groupby(group_cols)["_path"]
            .apply(lambda s: any(v for v in s))
            .reset_index(name="paths_present")
        )
        overview = overview.merge(paths_present, on=group_cols, how="left")
    else:
        overview["paths_present"] = False

    overview["metric_count"] = overview["metric_count"].fillna(0).astype(int)
    overview["metrics"] = overview["metrics"].fillna("")
    overview["labels_present"] = overview["labels_present"].fillna(False)
    overview["paths_present"] = overview["paths_present"].fillna(False)

    return overview.sort_values(["pair", "model"]).reset_index(drop=True)


def _print_metric_overview_table(scores_df: pd.DataFrame) -> pd.DataFrame:
    """Print the metric-agnostic overview table for unified scores."""
    tbl = metric_overview_table(scores_df)
    print()
    print("── Metric-Agnostic Overview ──")
    print(tbl.to_string(index=False))
    print()
    return tbl


def _print_tabular_summary(scores_df: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Print the tabular summary to stdout with a clear header.

    Returns the summary DataFrame so it can be reused (e.g. for plot tags).
    """
    tbl = tabular_summary(scores_df, metric)
    direction = get_metric_direction(metric)
    direction_label = "higher is better" if direction else "lower is better"
    print()
    print(f"── Tabular Summary: {metric} ({direction_label}) ──")
    print(tbl.to_string(index=False))
    print()
    return tbl


def _annotate_plot_source(
    fig: plt.Figure,
    metric: str,
    summary_df: pd.DataFrame,
) -> None:
    """Add a small footnote to the figure referencing the tabular summary.

    Provides traceability: every plot documents its source data at a
    glance — metric name, model count, and total observations.
    """
    n_models = summary_df["model"].nunique()
    n_obs = int(summary_df["count"].sum())
    note = (
        f"metric: {metric}  |  {n_models} engine(s)  |  "
        f"{n_obs} observations  |  see tabular summary"
    )
    fig.text(
        0.5, -0.01, note,
        ha="center", va="top", fontsize=7, fontstyle="italic",
        color="grey",
    )


# ---------------------------------------------------------------------------
# Score normalization
# ---------------------------------------------------------------------------

def normalize_scores(
    scores_df: pd.DataFrame,
    method: str = "minmax",
    per_model: bool = True,
    direction_aware: bool = True,
) -> pd.DataFrame:
    """Normalize ``score_value`` so that different engines are comparable.

    Different docking engines produce scores on completely different scales
    **and in different directions**:

    * LightDock luciferin scores ∼0–100, **higher = better**
    * HADDOCK scores ∼ −150 to +50, **lower = better**
    * Rosetta interface energies ∼ −30 to +10, **lower = better**

    Simply scaling them to [0, 1] is not enough — you also need to flip
    the lower-is-better metrics so that a high normalised value always
    means "this engine thinks the pair interacts strongly".

    When *direction_aware* is True (default), the function automatically
    flips lower-is-better metrics using :data:`METRIC_METADATA`.  After
    normalisation, **0 = worst / least interacting, 1 = best / most
    interacting** regardless of the original scoring convention.

    Parameters
    ----------
    scores_df : DataFrame
        Unified scores (must have ``score_value``; ``model`` and
        ``score_type`` are optional but recommended).
    method : ``"minmax"`` | ``"zscore"`` | ``"rank"``
        * **minmax** – scales each group to [0, 1].
        * **zscore** – mean-centres and divides by std (μ=0, σ=1).
        * **rank** – replaces values with percentile ranks in [0, 1].
    per_model : bool
        If *True* (default), normalize within each ``(model, score_type)``
        group independently.  If *False*, normalize across all rows for
        each ``score_type``.
    direction_aware : bool
        If *True* (default), flip lower-is-better metrics **before**
        normalizing so that after normalization higher always means
        "better binding".  Set to *False* for raw normalisation without
        flipping.

    Returns
    -------
    DataFrame
        A copy with ``score_value`` replaced by normalised values and
        ``norm_method`` recording the strategy used.
    """
    METHODS = ("minmax", "zscore", "rank")
    if method not in METHODS:
        raise ValueError(f"method must be one of {METHODS}, got '{method}'")

    df = scores_df.copy()
    df["score_value"] = pd.to_numeric(df["score_value"], errors="coerce")

    # Flip lower-is-better metrics so that higher always = better binding.
    # Without this step, comparing HADDOCK (lower = better) against
    # LightDock (higher = better) on the same normalized scale would be
    # misleading — one engine's "good" would map to the other's "bad".
    if direction_aware and "score_type" in df.columns:
        for st in df["score_type"].dropna().unique():
            if not get_metric_direction(st):
                mask = df["score_type"] == st
                df.loc[mask, "score_value"] = -df.loc[mask, "score_value"]

    group_cols = ["score_type"]
    if per_model and "model" in df.columns:
        group_cols = ["model", "score_type"]

    def _norm(g: pd.Series) -> pd.Series:
        if method == "minmax":
            lo, hi = g.min(), g.max()
            rng = hi - lo
            return (g - lo) / rng if rng != 0 else g * 0 + 0.5
        elif method == "zscore":
            mu, sigma = g.mean(), g.std(ddof=0)
            return (g - mu) / sigma if sigma != 0 else g * 0
        else:  # rank
            return g.rank(pct=True)

    if "score_type" in df.columns:
        df["score_value"] = df.groupby(group_cols)["score_value"].transform(_norm)
    else:
        df["score_value"] = _norm(df["score_value"])

    df["norm_method"] = method
    return df


def classification_summary(
    scores_df: pd.DataFrame,
    metric: str,
    threshold: float | None = None,
    higher_is_better: bool | None = None,
) -> pd.DataFrame:
    """Compute per-model confusion-matrix-style statistics.

    For each model, splits pairs into those above and below *threshold*
    on the given *metric*, then cross-tabulates against ``label``.

    If *threshold* is ``None``, the overall median of *metric* across all
    models is used.

    Parameters
    ----------
    scores_df : DataFrame
        Unified scores with ``label`` column.
    metric : str
        The ``score_type`` to threshold.
    threshold : float | None
        Score threshold to call "predicted interacting".
    higher_is_better : bool | None
        If True, scores ≥ threshold → predicted interacting.
        If False, scores ≤ threshold → predicted interacting.
        If None (default), auto-detected from :data:`METRIC_METADATA`.

    Returns
    -------
    DataFrame with columns: model, TP, FP, TN, FN, accuracy, precision,
    recall, f1.
    """
    if higher_is_better is None:
        higher_is_better = get_metric_direction(metric)
    df = scores_df.copy()
    if "score_type" in df.columns:
        df = df[df["score_type"].str.lower() == metric.lower()]
    if "label" not in df.columns:
        raise ValueError("Need 'label' column for classification summary")

    df["score_value"] = pd.to_numeric(df["score_value"], errors="coerce")
    df = df.dropna(subset=["score_value"])

    # Aggregate to one score per (model, proteinA, proteinB) pair
    agg = df.groupby(
        ["model", "proteinA", "proteinB", "label"],
        as_index=False,
    )["score_value"].mean()

    if threshold is None:
        threshold = agg["score_value"].median()

    if higher_is_better:
        agg["predicted"] = agg["score_value"] >= threshold
    else:
        agg["predicted"] = agg["score_value"] <= threshold

    agg["actual"] = agg["label"].str.lower() == "interaction"

    rows: list[dict] = []
    for model, grp in agg.groupby("model"):
        tp = ((grp["predicted"]) & (grp["actual"])).sum()
        fp = ((grp["predicted"]) & (~grp["actual"])).sum()
        tn = ((~grp["predicted"]) & (~grp["actual"])).sum()
        fn = ((~grp["predicted"]) & (grp["actual"])).sum()
        total = tp + fp + tn + fn
        acc = (tp + tn) / total if total else 0
        prec = tp / (tp + fp) if (tp + fp) else 0
        rec = tp / (tp + fn) if (tp + fn) else 0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0
        rows.append({
            "model": model, "threshold": threshold,
            "TP": int(tp), "FP": int(fp), "TN": int(tn), "FN": int(fn),
            "accuracy": round(acc, 3), "precision": round(prec, 3),
            "recall": round(rec, 3), "f1": round(f1, 3),
            "n_pairs": int(total),
        })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# CLI helpers
# ---------------------------------------------------------------------------

def _print_plot_guide(scores_df: pd.DataFrame) -> None:
    """Inspect *scores_df* and print which plot types are applicable.

    Checks the data for the prerequisites each plot type needs (number
    of models, number of pairs, presence of labels, etc.) and prints a
    concise table so the user can plan their analysis.
    """
    models = (
        sorted(scores_df["model"].dropna().unique())
        if "model" in scores_df.columns
        else []
    )
    n_models = len(models)

    pairs = available_pairs(scores_df)
    n_pairs = len(pairs)

    has_labels = "label" in scores_df.columns and scores_df["label"].notna().any()

    metrics = available_metrics(scores_df)
    n_metrics = len(metrics)

    # For scatter: need ≥2 models sharing a metric on the same pairs
    metrics_per_model: dict[str, set[str]] = {}
    pairs_per_model: dict[str, set[tuple[str, str]]] = {}
    if "score_type" in scores_df.columns:
        for model in models:
            model_df = scores_df[scores_df["model"] == model]
            metrics_per_model[model] = set(
                model_df["score_type"].dropna().str.lower().unique()
            )
            if "proteina" in scores_df.columns and "proteinb" in scores_df.columns:
                pairs_per_model[model] = set(
                    zip(
                        model_df["proteina"].dropna(),
                        model_df["proteinb"].dropna(),
                        strict=False,
                    )
                )
    shared_metrics: set[str] = set()
    shared_pairs: set[tuple[str, str]] = set()
    if len(metrics_per_model) >= 2:
        model_list = list(metrics_per_model.values())
        shared_metrics = model_list[0]
        for s in model_list[1:]:
            shared_metrics = shared_metrics & s
    if len(pairs_per_model) >= 2:
        pair_list = list(pairs_per_model.values())
        shared_pairs = pair_list[0]
        for s in pair_list[1:]:
            shared_pairs = shared_pairs & s

    # ── Header ────────────────────────────────────────────────────
    print()
    print("═══ ppinsight compare — plot-type guide ═══")
    print()
    models_str = (
        ", ".join(models) if n_models <= 5
        else ", ".join(models[:5]) + ", …"
    )
    print(f"  Models   : {n_models}  ({models_str})")
    print(f"  Pairs    : {n_pairs}  ", end="")
    if n_pairs == 0:
        print("(no proteinA/proteinB columns)")
    elif n_pairs <= 4:
        print("(" + ", ".join(f"{a}:{b}" for a, b in pairs) + ")")
    else:
        print(f"({pairs[0][0]}:{pairs[0][1]}, … +{n_pairs - 1} more)")
    print(f"  Labels   : {'yes' if has_labels else 'no'}")
    metrics_str = ", ".join(metrics[:6])
    if n_metrics > 6:
        metrics_str += ", …"
    print(f"  Metrics  : {n_metrics}  ({metrics_str})")
    print()

    # ── Per-plot-type assessment ──────────────────────────────────
    rows: list[tuple[str, str, str]] = []

    # violin
    rows.append((
        "violin", "✅ yes",
        "Works with any setup — full distributions "
        "+ individual points.",
    ))

    rows.append((
        "ridge", "✅ yes",
        "Works with any setup — compact density comparison across engines.",
    ))

    # roc
    if has_labels:
        rows.append((
            "roc", "✅ yes",
            "Labels present — ROC with AUC per model.",
        ))
    else:
        rows.append((
            "roc", "❌ no",
            "Needs a 'label' column. Re-run "
            "'ppinsight collect' with --pairs <pairs_file>.",
        ))

    # scatter — needs ≥2 models sharing a metric on the same pairs
    if n_models < 2:
        rows.append((
            "scatter", "❌ no",
            f"Needs ≥ 2 models (you have {n_models}). "
            "Use --models MODEL_A MODEL_B.",
        ))
    elif shared_metrics and shared_pairs:
        rows.append((
            "scatter",
            "✅ yes",
            f"Shared metrics: {', '.join(sorted(shared_metrics)[:5])}.  "
            f"Shared pairs: {len(shared_pairs)}.  Use --models to pick two."
        ))
    elif not shared_pairs:
        rows.append((
            "scatter", "❌ no",
            "Models have no common pairs — both must "
            "score the same protein pairs.",
        ))
    elif not shared_metrics:
        rows.append((
            "scatter", "⚠️  limited",
            "Models share pairs but no common metric "
            "name. Try --normalize first.",
        ))
    else:
        rows.append((
            "scatter", "❌ no",
            "Needs proteinA/proteinB columns. Re-run "
            "'ppinsight collect' with --pair.",
        ))

    # quality_bar — needs fnat/irmsd/lrmsd or dockq metrics
    has_quality_metrics = False
    if "score_type" in scores_df.columns:
        st_lower = set(scores_df["score_type"].dropna().str.lower().unique())
        quality_metrics = {
            "fnat", "irmsd", "lrmsd", "dockq",
            "quality_fnat", "quality_irmsd",
            "quality_lrmsd", "quality_dockq",
        }
        has_quality_metrics = bool(st_lower & quality_metrics)
    if has_quality_metrics:
        rows.append((
            "quality_bar", "✅ yes",
            "Pose-level CAPRI quality tiers (fnat/RMSD/dockq "
            "detected). No --metric needed.",
        ))
    else:
        rows.append((
            "quality_bar", "❌ no",
            "Needs fnat + irmsd/lrmsd, or dockq metrics "
            "for CAPRI quality tiers.",
        ))

    rows.append((
        "cdf", "✅ yes",
        "Works with any numeric metric; DockQ adds CAPRI threshold lines.",
    ))

    if n_models < 2:
        rows.append((
            "difference", "❌ no",
            f"Needs ≥ 2 models (you have {n_models}) and shared pairs.",
        ))
    elif shared_metrics and shared_pairs:
        rows.append((
            "difference", "✅ yes",
            f"Shared metrics: {', '.join(sorted(shared_metrics)[:5])}.  "
            f"Shared pairs: {len(shared_pairs)}.  Use --models to pick two.",
        ))
    else:
        rows.append((
            "difference", "❌ no",
            "Needs two models that scored the same metric on the same pairs.",
        ))

    # ── Print table ───────────────────────────────────────────────
    col_w = [max(len(r[i]) for r in rows) for i in range(3)]
    col_w[0] = max(col_w[0], len("plot type"))
    col_w[1] = max(col_w[1], len("available?"))
    col_w[2] = max(col_w[2], len("details"))

    hdr = f"  {'plot type':<{col_w[0]}}  {'available?':<{col_w[1]}}  {'details'}"
    sep = f"  {'─' * col_w[0]}  {'─' * col_w[1]}  {'─' * col_w[2]}"
    print(hdr)
    print(sep)
    for name, status, detail in rows:
        print(f"  {name:<{col_w[0]}}  {status:<{col_w[1]}}  {detail}")

    # ── Extras ────────────────────────────────────────────────────
    print()
    extras: list[str] = []
    if has_labels:
        extras.append(
            "--classify       Print confusion-matrix "
            "summaries (uses labels)."
        )
        extras.append(
            "--split-label    Split violins by "
            "interaction label."
        )
    else:
        extras.append(
            "--classify       ❌ needs labels "
            "(collect with --pairs)."
        )
        extras.append(
            "--split-label    ❌ needs labels "
            "(collect with --pairs)."
        )
    extras.append(
        "--normalize      Normalise scores for "
        "cross-engine comparison (always available)."
    )
    extras.append(
        "--table          Tabular summary "
        "for a metric, or metric-agnostic overview when --metric is omitted."
    )
    extras.append(
        "--rank           Per-pair ranking: which "
        "engine scores best on each pair."
    )
    if has_quality_metrics:
        extras.append(
            "--capri-quality  CAPRI quality tier summary "
                "(high/medium/acceptable/incorrect) across poses."
        )
    else:
        extras.append(
            "--capri-quality  ❌ needs fnat/irmsd/lrmsd "
            "or dockq metrics."
        )
    print("  Additional flags:")
    for e in extras:
        print(f"    {e}")
    print()


def _cli_error_with_hint(msg: str, metric: str, scores_df) -> None:
    """Print an error message with a context-aware hint, then exit.

    Matches common ValueError messages from plotting / classification
    functions and suggests the CLI flag or upstream step that would
    resolve the problem.
    """
    print(f"ERROR: {msg}", file=sys.stderr)
    lower = msg.lower()

    if "label" in lower:
        # Missing interaction labels — need to annotate at collection time
        print(
            "Hint: your scores file has no 'label' column.  Re-run "
            "'ppinsight collect' with --pairs <pairs_file> to annotate "
            "rows with interaction / non-interaction labels.",
            file=sys.stderr,
        )
    elif "no rows match metric" in lower:
        # Metric not found — show available metrics
        avail = []
        if scores_df is not None and "score_type" in scores_df.columns:
            avail = sorted(scores_df["score_type"].dropna().unique())
        print(
            f"Hint: metric '{metric}' was not found.  Use --list-metrics "
            f"to see available values.",
            file=sys.stderr,
        )
        if avail:
            print(f"  Available metrics: {', '.join(avail)}", file=sys.stderr)
    elif "no common pairs" in lower:
        print(
            "Hint: the two models have no overlapping protein pairs.  "
            "Check that both were scored on the same pair set.",
            file=sys.stderr,
        )
    elif "proteinA" in lower or "proteinB" in lower:
        print(
            "Hint: this plot needs proteinA and proteinB columns.  "
            "Re-run 'ppinsight collect' with --pair or --pairs.",
            file=sys.stderr,
        )
    sys.exit(1)


def _parse_pair_arg(raw_pair: str | None) -> tuple[str, str] | None:
    """Parse ``proteinA:proteinB`` input from the CLI."""
    if not raw_pair:
        return None
    parts = [p.strip() for p in raw_pair.split(":")]
    if len(parts) != 2 or not parts[0] or not parts[1]:
        print("ERROR: --pair must be formatted as proteinA:proteinB", file=sys.stderr)
        sys.exit(2)
    return parts[0], parts[1]


def _slugify_output_token(text: str) -> str:
    """Return a filesystem-safe token for auto-generated output filenames."""
    token = "".join(ch if ch.isalnum() or ch in {"_", "-", "."} else "_" for ch in text)
    token = token.strip("._")
    return token or "plot"


def _next_available_output_path(path: str) -> str:
    """Return *path* if free, else append ``_N`` before file extension."""
    if not os.path.exists(path):
        return path

    stem, ext = os.path.splitext(path)
    idx = 1
    while True:
        candidate = f"{stem}_{idx}{ext}"
        if not os.path.exists(candidate):
            return candidate
        idx += 1


def _default_plot_output_path(
    plot_type: str,
    metric: str | None,
    pair: tuple[str, str] | None,
) -> str:
    """Build a default plot file path under ``data/output/plots``."""
    out_dir = os.path.join("data", "output", "plots")
    metric_token = _slugify_output_token(metric or "summary")
    name_parts = [_slugify_output_token(plot_type), metric_token]
    if pair is not None:
        pair_token = (
            f"{_slugify_output_token(pair[0])}_vs_{_slugify_output_token(pair[1])}"
        )
        name_parts.append(pair_token)
    base = "_".join(name_parts) + ".png"
    return _next_available_output_path(os.path.join(out_dir, base))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    """CLI entry-point for score comparison plots.

    Examples::

        ppinsight compare scores.tsv --metric dockq
        ppinsight compare model_a.csv model_b.csv --metric score --names HADDOCK Rosetta
    """
    parser = argparse.ArgumentParser(
        prog="ppinsight compare",
        description="Compare docking scores across PPI prediction models.",
    )
    parser.add_argument(
        "files",
        nargs="+",
        help=(
            "Path(s) to score file(s).  A single file is treated as a unified "
            "scores file (columns: model, score_type, score_value, …) "
            "produced by 'ppinsight collect'.  Multiple files are treated as "
            "one-per-model (legacy / quick-compare mode)."
        ),
    )
    parser.add_argument(
        "--metric", "-m",
        default=None,
        help=(
            "Score type to analyse (e.g. 'dockq', 'score', 'irmsd', "
            "'luciferin_score', 'i_sc').  Must match a value in the "
            "'score_type' column of the scores file.  Use --list-metrics "
            "to see what is available.  Required for plotting and "
            "classification; not needed for --list-metrics, --list-pairs, "
            "--guide, or metric-agnostic --table mode."
        ),
    )
    parser.add_argument(
        "--names", "-n",
        nargs="*",
        default=None,
        help=(
            "Model display names (one per file).  Required when passing "
            "multiple per-model files so the plot legend is meaningful.  "
            "Ignored in unified-file mode (names come from the 'model' "
            "column)."
        ),
    )
    parser.add_argument(
        "--pair", "-p",
        default=None,
        help=(
            "Filter to a single protein pair, formatted as "
            "'proteinA:proteinB'.  Use when the scores file contains many "
            "pairs and you want to drill into one specific interaction."
        ),
    )
    parser.add_argument(
        "--output", "-o",
        default=None,
        help=(
            "Save the plot to a file (png/svg/pdf).  If omitted, compare "
            "auto-saves a PNG under data/output/plots/ using plot type and "
            "metric in the filename.  Parent directories are created "
            "automatically."
        ),
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help=(
            "Display plots interactively instead of auto-saving when -o is "
            "not provided."
        ),
    )
    parser.add_argument(
        "--list-metrics",
        action="store_true",
        help=(
            "Print all available score_type values in the scores file and "
            "exit.  Use this to discover which metrics were collected."
        ),
    )
    parser.add_argument(
        "--list-pairs",
        action="store_true",
        help=(
            "Print all available proteinA:proteinB pairs in the scores file "
            "and exit.  Use this to discover pair names for --pair filtering."
        ),
    )
    parser.add_argument(
        "--guide",
        action="store_true",
        help=(
            "Inspect the scores file and print which plot types are "
            "applicable given your data (number of models, pairs, labels, "
            "metrics).  Useful before choosing --plot-type.  Exit after "
            "printing."
        ),
    )
    parser.add_argument(
        "--split-label",
        action="store_true",
        help=(
            "Split violins by interaction label (interaction vs "
            "non-interaction).  Requires a 'label' column (added by "
            "'ppinsight collect --pairs').  Lets you visually assess "
            "whether an engine separates known binders from non-binders."
        ),
    )
    parser.add_argument(
        "--classify",
        action="store_true",
        help=(
            "Print a confusion-matrix summary (TP, FP, TN, FN, accuracy, "
            "precision, recall, F1) for each model.  Requires a 'label' "
            "column.  The pipeline automatically determines whether higher "
            "or lower scores indicate binding for each metric.  "
            "Default threshold: median score (override with --threshold)."
        ),
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=None,
        help=(
            "Score threshold for --classify (default: median of the metric).  "
            "Set this when you have domain knowledge of a meaningful cutoff "
            "(e.g. DockQ ≥ 0.23 for 'acceptable' quality).  Omit to let "
            "the pipeline use the data-driven median."
        ),
    )
    parser.add_argument(
        "--normalize",
        choices=["minmax", "zscore", "rank"],
        default=None,
        help=(
            "Normalize score_value before plotting or classifying.  "
            "'minmax' scales to [0, 1], 'zscore' centres at μ=0 σ=1, "
            "'rank' converts to percentile ranks in [0, 1].  "
            "Normalization automatically flips lower-is-better metrics "
            "(HADDOCK score, Rosetta energies) so that higher always means "
            "stronger predicted binding — this makes cross-engine "
            "comparisons on the same axis meaningful."
        ),
    )
    parser.add_argument(
        "--plot-type",
        choices=["violin", "ridge", "roc", "scatter",
                 "quality_bar", "cdf", "difference"],
        default="violin",
        help=(
            "Plot style.  'violin' (default): full score distributions "
            "with individual data points — reveals bimodality or outliers.  "
            "'ridge': ridgeline KDE per engine — compact when comparing "
            "many engines.  'roc': ROC curves with AUC — needs labels.  "
            "'scatter': model-agreement plot — same metric, same pairs, "
            "one model per axis (use with --models).  "
            "'quality_bar': 100 %% stacked bar of CAPRI quality tiers "
            "per engine — requires fnat+RMSD or dockq.  "
            "'cdf': empirical CDF per engine — CAPRI thresholds shown "
            "for dockq.  "
            "'difference': histogram of per-pair mean score differences "
            "between two engines (use with --models)."
        ),
    )
    parser.add_argument(
        "--models",
        nargs=2,
        metavar=("MODEL_A", "MODEL_B"),
        default=None,
        help=(
            "Two model names for --plot-type scatter (e.g. --models "
            "haddock rosetta).  Each protein pair becomes one point: "
            "x = model A's score, y = model B's.  Points along the "
            "diagonal mean the engines agree.  Pearson r is annotated."
        ),
    )
    parser.add_argument(
        "--table",
        action="store_true",
        help=(
            "Print tables only (no plot).  With --metric, prints the "
            "per-metric summary (mean/std/median/min/max).  Without "
            "--metric, prints a metric-agnostic overview of models, "
            "pairs, and available metrics."
        ),
    )
    parser.add_argument(
        "--rank",
        action="store_true",
        help=(
            "Print a per-pair ranking table: which engine scores best on "
            "each protein pair for the given --metric.  Respects score "
            "direction (lower HADDOCK score = rank 1, higher DockQ = rank 1)."
        ),
    )
    parser.add_argument(
        "--capri-quality",
        action="store_true",
        help=(
            "Print a CAPRI quality tier summary: count of models classified "
            "as high / medium / acceptable / incorrect for each engine.  "
            "Uses the standard CAPRI protocol: f(nat) + Lrms + irms "
            "thresholds (Lensink et al. 2016, DOI: 10.1002/prot.25007, "
            "Table 3).  Falls back to DockQ thresholds (Basu & Wallner 2016, "
            "DOI: 10.1371/journal.pone.0161879) when RMSD values are absent.  "
            "Does not require --metric."
        ),
    )

    args = parser.parse_args(argv)
    pair_filter = _parse_pair_arg(args.pair)

    # Apply the default theme before any plotting
    apply_theme()

    # ── unified file mode (single file) ──────────────────────────
    if len(args.files) == 1:
        try:
            scores_df = load_scores(args.files[0])
        except ValueError:
            # Not a unified file — fall through to per-model mode with one file
            scores_df = None

        if scores_df is not None:
            # Apply normalization if requested.  Direction-aware flipping
            # is always enabled — the pipeline knows which metrics are
            # lower-is-better from METRIC_METADATA.  Per-model normalization
            # is the default (each engine's internal distribution is preserved).
            if args.normalize:
                scores_df = normalize_scores(
                    scores_df,
                    method=args.normalize,
                    per_model=True,
                    direction_aware=True,
                )

            if args.list_metrics:
                for m in available_metrics(scores_df):
                    print(m)
                return

            if args.list_pairs:
                for a, b in available_pairs(scores_df):
                    print(f"{a}:{b}")
                return

            if args.guide:
                _print_plot_guide(scores_df)
                return

            analysis_df = scores_df
            if pair_filter is not None:
                pa, pb, _pairs = _detect_pairs(analysis_df)
                if not pa or not pb:
                    print(
                        "ERROR: --pair requires proteinA/proteinB columns in "
                        "the scores file.",
                        file=sys.stderr,
                    )
                    print(
                        "Hint: re-run 'ppinsight collect' with --pair or --pairs.",
                        file=sys.stderr,
                    )
                    sys.exit(2)
                a, b = pair_filter
                analysis_df = analysis_df[
                    analysis_df[pa].fillna("").astype(str).str.lower().eq(a.lower())
                    & analysis_df[pb].fillna("").astype(str).str.lower().eq(b.lower())
                ]
                if analysis_df.empty:
                    _cli_error_with_hint(
                        f"No rows match pair={pair_filter}",
                        args.metric or "",
                        scores_df,
                    )

            # --capri-quality does NOT need --metric
            if args.capri_quality:
                try:
                    summary = capri_summary_table(analysis_df)
                except ValueError as exc:
                    print(f"ERROR: {exc}", file=sys.stderr)
                    sys.exit(1)
                print(summary.to_string(index=False))
                return

            if args.table and not args.metric:
                try:
                    _print_metric_overview_table(analysis_df)
                except ValueError as exc:
                    print(f"ERROR: {exc}", file=sys.stderr)
                    sys.exit(1)
                return

            # --metric is required for everything below this point
            # (unless --plot-type quality_bar, which also doesn't need it).
            if args.plot_type == "quality_bar":
                plot_out = args.output
                if plot_out is None and not args.show:
                    plot_out = _default_plot_output_path(
                        "quality_bar",
                        "capri_quality",
                        pair_filter,
                    )
                if plot_out is not None:
                    os.makedirs(os.path.dirname(plot_out) or ".", exist_ok=True)
                try:
                    quality_bar_chart(
                        analysis_df,
                        plot_title=None,
                        output=plot_out,
                    )
                except ValueError as exc:
                    print(f"ERROR: {exc}", file=sys.stderr)
                    sys.exit(1)
                return

            if not args.metric:
                print(
                    "ERROR: --metric is required for plotting and "
                    "classification.  Use --list-metrics to see available "
                    "values, or --guide to see which plot types work with "
                    "your data.  If you only need a metric-agnostic table, "
                    "use --table without --metric.",
                    file=sys.stderr,
                )
                sys.exit(2)

            # ── Mandatory tabular summary ─────────────────────────
            # Always printed before any metric-based output.  Every
            # plot is a *visualisation* of this summary — the table is
            # the canonical data reference.
            try:
                summary_tbl = _print_tabular_summary(analysis_df, args.metric)
            except ValueError as exc:
                _cli_error_with_hint(str(exc), args.metric, scores_df)

            if args.classify:
                try:
                    summary = classification_summary(
                        analysis_df,
                        metric=args.metric,
                        threshold=args.threshold,
                    )
                except ValueError as exc:
                    _cli_error_with_hint(str(exc), args.metric, scores_df)
                print(summary.to_string(index=False))
                return

            # --table: tabular summary only (already printed above)
            if args.table:
                return

            if args.rank:
                try:
                    rnk = ranking_table(analysis_df, metric=args.metric)
                except ValueError as exc:
                    _cli_error_with_hint(str(exc), args.metric, scores_df)
                print(rnk.to_string(index=False))
                return

            # Dispatch by plot type.
            # Plot functions receive the resolved output path so they avoid
            # interactive rendering when auto-save mode is active.
            pt = args.plot_type

            def _resolve_plot_output(plot_type: str) -> str | None:
                plot_out = args.output
                if plot_out is None and not args.show:
                    plot_out = _default_plot_output_path(
                        plot_type,
                        args.metric,
                        pair_filter,
                    )
                if plot_out is not None:
                    os.makedirs(os.path.dirname(plot_out) or ".", exist_ok=True)
                return plot_out

            def _save_or_show(fig, plot_out: str | None):
                """Annotate and re-save plots when writing to disk."""
                if not plot_out:
                    return
                _annotate_plot_source(fig, args.metric, summary_tbl)
                os.makedirs(os.path.dirname(plot_out) or ".", exist_ok=True)
                fig.savefig(plot_out, bbox_inches="tight", dpi=180)

            try:
                if pt == "violin":
                    plot_out = _resolve_plot_output(pt)
                    fig = violin_plot(
                        analysis_df,
                        metric=args.metric,
                        split_by_label=args.split_label,
                        plot_title=None,
                        output=plot_out,
                    )
                    _save_or_show(fig, plot_out)
                    return

                if pt == "roc":
                    plot_out = _resolve_plot_output(pt)
                    fig = roc_curve_plot(
                        analysis_df,
                        metric=args.metric,
                        output=plot_out,
                    )
                    _save_or_show(fig, plot_out)
                    return

                if pt == "scatter":
                    if not args.models or len(args.models) != 2:
                        print(
                            "ERROR: --plot-type scatter requires "
                            "--models MODEL_A MODEL_B",
                            file=sys.stderr,
                        )
                        print(
                            "Hint: pass two model names that both "
                            "scored the same pairs, e.g. "
                            "--models haddock rosetta",
                            file=sys.stderr,
                        )
                        sys.exit(2)
                    plot_out = _resolve_plot_output(pt)
                    fig = model_agreement_scatter(
                        analysis_df,
                        metric=args.metric,
                        model_x=args.models[0],
                        model_y=args.models[1],
                        output=plot_out,
                    )
                    _save_or_show(fig, plot_out)
                    return

                if pt == "ridge":
                    plot_out = _resolve_plot_output(pt)
                    fig = ridge_plot(
                        analysis_df,
                        metric=args.metric,
                        plot_title=None,
                        output=plot_out,
                    )
                    _save_or_show(fig, plot_out)
                    return

                if pt == "cdf":
                    plot_out = _resolve_plot_output(pt)
                    fig = cdf_plot(
                        analysis_df,
                        metric=args.metric,
                        plot_title=None,
                        output=plot_out,
                    )
                    _save_or_show(fig, plot_out)
                    return

                if pt == "difference":
                    if not args.models or len(args.models) != 2:
                        print(
                            "ERROR: --plot-type difference requires "
                            "--models MODEL_A MODEL_B",
                            file=sys.stderr,
                        )
                        sys.exit(2)
                    plot_out = _resolve_plot_output(pt)
                    fig = pairwise_difference_plot(
                        analysis_df,
                        metric=args.metric,
                        model_a=args.models[0],
                        model_b=args.models[1],
                        output=plot_out,
                    )
                    _save_or_show(fig, plot_out)
                    return

                raise ValueError(f"Unsupported plot type '{pt}'")
            except ValueError as exc:
                _cli_error_with_hint(str(exc), args.metric, scores_df)
            return

    # ── per-model file mode (multiple files) ─────────────────────
    if not args.metric:
        print(
            "ERROR: --metric is required for plotting.  Pass a column name "
            "from your score files.",
            file=sys.stderr,
        )
        sys.exit(2)
    frames = to_plot(args.files)
    names = args.names
    if not names:
        # Derive names from file stems
        names = [os.path.splitext(os.path.basename(f))[0] for f in args.files]
    if len(names) != len(frames):
        print(
            f"ERROR: {len(names)} names given for {len(frames)} files.",
            file=sys.stderr,
        )
        sys.exit(1)

    output = args.output
    if output is None and not args.show:
        output = _default_plot_output_path("bar", args.metric, None)
    compare_scores(frames, names, args.metric, plot_title=None, output=output)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _guess_sep(path: str) -> str:
    """Return '\\t' for .tsv files, ',' otherwise."""
    if path.lower().endswith(".tsv"):
        return "\t"
    return ","


if __name__ == "__main__":
    main()
