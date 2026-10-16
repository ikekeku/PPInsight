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
from contextlib import contextmanager

import pandas as pd
import seaborn as sns
from matplotlib import colors as mcolors
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# Theme configuration
# ---------------------------------------------------------------------------

# Okabe & Ito (2008) colorblind-safe palette — PPInsight's default for
# every plot.  Yellow is last because it has poor contrast on white.
OKABE_ITO = {
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "bluish_green": "#009E73",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "reddish_purple": "#CC79A7",
    "black": "#000000",
    "yellow": "#F0E442",
}

# PPInsight has one theme; there is no theme option.
THEMES: dict[str, dict] = {
    "white-okabe-ito": {"style": "white", "palette": list(OKABE_ITO.values())},
}

DEFAULT_THEME = "white-okabe-ito"

PLOT_TITLE_SIZE = 14
FACET_TITLE_SIZE = 11
AXIS_LABEL_SIZE = 12
TICK_LABEL_SIZE = 10
LEGEND_FONT_SIZE = 8
LEGEND_TITLE_SIZE = 9
# Figure-level legends placed below the plot (violin, strip, funnel,
# quality_bar) have room to be larger than legends inside the axes.
BELOW_LEGEND_FONT_SIZE = 10
BELOW_LEGEND_TITLE_SIZE = 11
VIOLIN_FILL_ALPHA = 0.55
ANNOTATION_FONT_SIZE = 9
TITLE_PAD = 10
BAR_ERROR_KW = {"linewidth": 1.0, "capthick": 1.0}

# Below this many shared pairs, correlation coefficients swing widely from
# one pair to the next, so the scatter still draws but warns.
MIN_SCATTER_PAIRS = 10

# Violin smoothing: 0.75 × seaborn's default (Scott's rule) bandwidth.  The
# default filled the gap between two-cluster DockQ distributions; 0.5 kept
# the gap but looked too bumpy.
VIOLIN_BW_ADJUST = 0.75

# Most points drawn per engine in strip and violin plots.  Above this a
# uniform random sample is drawn (so density stays honest) and rare good
# poses are overlaid; statistics still use all poses.
MAX_DRAWN_POSES = 5000
ALWAYS_DRAWN_BEST = 50

# Metrics on the CAPRI DockQ scale: HADDOCK's own ``dockq`` and the
# native-referenced ``quality_dockq`` written by ``ppinsight quality``.
DOCKQ_METRICS = {"dockq", "quality_dockq"}

# One colour per CAPRI tier, shared by the quality_bar segments and the
# DockQ threshold lines so a tier looks the same in every plot.
CAPRI_TIER_COLORS = {
    "incorrect": "#b0b0b0",
    "acceptable": "#31a354",
    "medium": "#3182bd",
    "high": "#6a51a3",
}

DOCKQ_CAPRI_THRESHOLDS = [
    (0.23, "acceptable", CAPRI_TIER_COLORS["acceptable"]),
    (0.49, "medium", CAPRI_TIER_COLORS["medium"]),
    (0.80, "high", CAPRI_TIER_COLORS["high"]),
]

ENGINE_COLORS = {
    "haddock": OKABE_ITO["vermillion"],
    "lightdock": OKABE_ITO["blue"],
    "rosetta": OKABE_ITO["reddish_purple"],
}

# Interaction labels keep the same colors in every plot that uses them.
LABEL_COLORS = {
    "interaction": OKABE_ITO["blue"],
    "non-interaction": OKABE_ITO["orange"],
}

ENGINE_DISPLAY_NAMES = {
    "haddock": "HADDOCK3",
    "lightdock": "LightDock",
    "rosetta": "RosettaDock",
}


def _engine_label(model: str) -> str:
    """Return the publication name for an engine (e.g. ``rosetta`` → RosettaDock)."""
    return ENGINE_DISPLAY_NAMES.get(str(model).lower(), str(model))


def _engine_tick_label(
    model: str, n: int | None = None, shown: int | None = None,
) -> str:
    """Engine name for an axis tick, with the pose count underneath."""
    label = _engine_label(model)
    if n is None:
        return label
    if shown is not None and shown < n:
        return f"{label}\n(n = {n:,}; {shown:,} sampled)"
    return f"{label}\n(n = {n:,})"


def _resolve_model_name(name: str, available) -> str:
    """Match a user-typed engine name to a ``model`` value in the data.

    Accepts the raw name (``haddock``) or the display name (``HADDOCK3``),
    case-insensitively.  Returns *name* unchanged when nothing matches so
    the caller's own "not found" error still fires.
    """
    wanted = str(name).lower()
    for model in available:
        if wanted in {str(model).lower(), _engine_label(model).lower()}:
            return model
    return name


def _with_engine_labels(df: pd.DataFrame) -> pd.DataFrame:
    """Copy of *df* with the ``model`` column shown as engine display names."""
    out = df.copy()
    if "model" in out.columns:
        out["model"] = out["model"].map(_engine_label)
    return out


def _point_alpha(n: int) -> float:
    """Point opacity that stays visible for small n and readable for large n."""
    return float(min(0.8, max(0.08, 6.0 / max(n, 1) ** 0.5)))


def _point_size(n: int) -> float:
    """Marker area that shrinks as the number of poses grows."""
    if n <= 300:
        return 16.0
    if n <= 2000:
        return 8.0
    return 4.0


def _capri_legend_handles() -> list[Line2D]:
    """Legend entries for the CAPRI DockQ threshold lines."""
    return [
        Line2D([0], [0], color=color, lw=2.0, ls="--",
               label=f"CAPRI {label} ({threshold:.2f})")
        for threshold, label, color in DOCKQ_CAPRI_THRESHOLDS
    ]


def _draw_capri_lines(
    ax, metric: str, orientation: str = "h", zorder: float = 1.5,
) -> bool:
    """Draw CAPRI DockQ thresholds when *metric* is on the DockQ scale.

    Returns True when lines were drawn so callers know to add a legend.
    """
    if metric.lower() not in DOCKQ_METRICS:
        return False
    draw = ax.axhline if orientation == "h" else ax.axvline
    for threshold, _, color in DOCKQ_CAPRI_THRESHOLDS:
        draw(threshold, color=color, lw=1.6, ls="--", alpha=0.75, zorder=zorder)
    return True


def _soften_violin_bodies(ax) -> None:
    """Make violin fills translucent while keeping their outlines solid.

    Threshold lines drawn behind the bodies then show through, as they did
    in the ridge plot.
    """
    from matplotlib.collections import PolyCollection

    for body in ax.collections:
        if isinstance(body, PolyCollection):
            faces = body.get_facecolors()
            if len(faces):
                faces[:, 3] = VIOLIN_FILL_ALPHA
                body.set_facecolors(faces)


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

    Known engines always get their fixed Okabe-Ito color from
    :data:`ENGINE_COLORS` (matched case-insensitively), regardless of
    which other engines are plotted or which theme is active.  Unknown
    engines receive the unreserved Okabe-Ito colors in sorted order,
    cycling if there are more engines than colors.
    """
    unique_models = sorted(set(models))
    reserved = set(ENGINE_COLORS.values())
    spare = [c for c in OKABE_ITO.values() if c not in reserved]

    palette: dict[str, tuple[float, float, float]] = {}
    n_unknown = 0
    for model in unique_models:
        fixed = ENGINE_COLORS.get(str(model).lower())
        if fixed is None:
            fixed = spare[n_unknown % len(spare)]
            n_unknown += 1
        palette[model] = mcolors.to_rgb(fixed)

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
                "distribution or --plot-type strip to see every pose.",
                stacklevel=2,
            )

    metric_name = get_metric_display_name(score_type)
    model_palette = _model_palette(models)
    palette = [model_palette[m] for m in models]
    bar_width = max(0.45, min(0.8, 0.65 * len(models)))
    fig, ax = plt.subplots(figsize=(max(6, len(models) * 2.5), 5))
    ax.bar([_engine_label(m) for m in models], means, width=bar_width,
           yerr=stds, capsize=6, color=palette, edgecolor="white",
           linewidth=0.8, error_kw=BAR_ERROR_KW)
    ax.set_xlabel("Interaction model")
    ax.set_ylabel(f"Mean {metric_name} (± 1 SD)")
    ax.set_title(plot_title or f"Mean {metric_name} by engine")
    _apply_metric_bounds(ax, score_type)
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
    metric_name = get_metric_display_name(metric)

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
            model_palette = _model_palette(model_names)
            palette = [model_palette[m] for m in model_names]
            bar_w = max(0.45, min(0.8, 0.65 * len(model_names)))
            ax.bar([_engine_label(m) for m in model_names], means, width=bar_w,
                   yerr=stds, capsize=6, color=palette, edgecolor="white",
                   linewidth=0.8, error_kw=BAR_ERROR_KW)
            ax.set_title(_pair_label(a, b), fontsize=FACET_TITLE_SIZE)
            ax.set_xlabel("")
            ax.set_ylabel(f"Mean {metric_name} (± 1 SD)" if idx == 0 else "")
            _apply_metric_bounds(ax, metric)
            sns.despine(ax=ax)

        _finish_figure(
            fig, axes,
            title=plot_title or f"Mean {metric_name} by engine",
            xlabel="Docking Engine",
        )

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

    if pair is None and len(pair_tuples) == 1:
        pair = pair_tuples[0]
    title = plot_title or (
        f"Mean {metric_name} by engine"
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
    user always knows which protein pair(s) are being evaluated.  The
    smoothing bandwidth is half of seaborn's default
    (:data:`VIOLIN_BW_ADJUST`) so two-cluster distributions keep their gap,
    and bounded metrics such as DockQ always span their full 0–1 range.

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
    if len(pair_tuples) >= 2:
        panels = [
            (_pair_label(a, b), df[(df[pa] == a) & (df[pb] == b)])
            for a, b in pair_tuples
        ]
        panel_w = max(4.5, len(models) * 1.9)
    else:
        panels = [(None, df)]
        panel_w = max(7, len(models) * 2.2)

    fig, axes = plt.subplots(
        1, len(panels), figsize=(panel_w * len(panels), 6),
        sharey=True, squeeze=False,
    )
    axes = list(axes[0])
    capri_drawn = False
    highlight_handle = None

    for idx, (ax, (panel_title, data)) in enumerate(zip(axes, panels, strict=False)):
        panel_models = sorted(data["model"].unique())
        shown = None
        if has_label:
            sns.violinplot(
                data=data, x="model", y="score_value", hue="label",
                split=True, inner="box", cut=0, bw_adjust=VIOLIN_BW_ADJUST,
                ax=ax, linewidth=0.8, density_norm="width", order=panel_models,
            )
            if idx == len(panels) - 1:
                ax.legend(
                    title="Label",
                    fontsize=LEGEND_FONT_SIZE,
                    title_fontsize=LEGEND_TITLE_SIZE,
                )
            elif ax.get_legend() is not None:
                ax.get_legend().remove()
        else:
            sns.violinplot(
                data=data, x="model", y="score_value", hue="model",
                inner="box", cut=0, bw_adjust=VIOLIN_BW_ADJUST,
                ax=ax, linewidth=0.8, density_norm="width", order=panel_models,
                hue_order=panel_models, palette=model_palette,
                dodge=False, legend=False,
            )
            shown, handle = _draw_deterministic_stripplot(
                ax=ax, data=data, order=panel_models, metric=metric,
            )
            highlight_handle = handle or highlight_handle
        _set_engine_ticks(ax, data, panel_models, shown)
        _soften_violin_bodies(ax)
        capri_drawn = _draw_capri_lines(ax, metric, zorder=0.5) or capri_drawn
        _apply_metric_bounds(ax, metric)
        if panel_title:
            ax.set_title(panel_title, fontsize=FACET_TITLE_SIZE)
        ax.set_xlabel("")
        ax.set_ylabel(metric_name if idx == 0 else "")
        _add_faint_horizontal_gridlines(ax)
        sns.despine(ax=ax)

    _finish_figure(
        fig, axes,
        title=plot_title or _engines_title(
            f"{metric_name} distribution", models, pair_tuples,
        ),
        xlabel="Docking Engine",
        legend_handles=(
            ([highlight_handle] if highlight_handle else [])
            + (_capri_legend_handles() if capri_drawn else [])
        ) or None,
        legend_ncol=2 if highlight_handle else 3,
    )

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

    _require_reference_free(metric, "A ROC curve")

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
        ax.plot(fpr_list, tpr_list,
                label=f"{_engine_label(model)} (AUC={auc:.3f})",
                linewidth=2.2, color=model_palette[model])

    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    n_pairs = len(agg[["proteina", "proteinb"]].drop_duplicates())
    ax.set_title(
        f"ROC Curve — {metric_name} ({n_pairs} protein pairs)",
        fontsize=PLOT_TITLE_SIZE, pad=TITLE_PAD,
    )
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

    model_x = _resolve_model_name(model_x, agg["model"].unique())
    model_y = _resolve_model_name(model_y, agg["model"].unique())
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
    label_x, label_y = _engine_label(model_x), _engine_label(model_y)

    few_pairs = len(common) < MIN_SCATTER_PAIRS
    if few_pairs:
        import warnings as _warnings

        _warnings.warn(
            f"Only {len(common)} protein pair(s) are shared by {label_x} and "
            f"{label_y}.  With fewer than {MIN_SCATTER_PAIRS} pairs, Pearson r "
            "and Spearman ρ can change a lot from a single pair, so treat the "
            "agreement shown as a rough impression, not a measurement.",
            stacklevel=2,
        )

    fig, ax = plt.subplots(figsize=(7, 7))

    # Shared, square axes so the y = x diagonal sits at 45°.
    bounds = get_metric_bounds(metric)
    if bounds is not None:
        lo, hi = bounds
    else:
        lo = float(min(x_vals.min(), y_vals.min()))
        hi = float(max(x_vals.max(), y_vals.max()))
        pad = (hi - lo) * 0.06 or 0.1
        lo, hi = lo - pad, hi + pad
    ax.plot([lo, hi], [lo, hi], "--", color=".55", lw=1.2, zorder=1,
            label="Perfect agreement (y = x)")

    point_style = dict(s=55, alpha=0.85, edgecolor="white", linewidth=0.7,
                       zorder=3, clip_on=False)
    if "label" in dx.columns:
        # Colour by label if available
        labels = dx.loc[common, "label"].values
        plot_df = pd.DataFrame({"x": x_vals, "y": y_vals, "label": labels})
        for label_value, grp in plot_df.groupby("label"):
            ax.scatter(
                grp["x"], grp["y"],
                color=LABEL_COLORS.get(str(label_value).lower(), ".3"),
                label=str(label_value), **point_style,
            )
    else:
        ax.scatter(x_vals, y_vals, color=".25", label="Protein pair",
                   **point_style)

    # Name the pairs the engines disagree on most.
    if len(common) >= 5:
        gaps = np.abs(y_vals - x_vals)
        for idx in np.argsort(gaps)[::-1][:3]:
            pair_a, pair_b = common[idx]
            # Point the label toward the centre so it stays inside the axes.
            mid = (lo + hi) / 2
            left_half = x_vals[idx] < mid
            lower_half = y_vals[idx] < mid
            ax.annotate(
                _pair_label(pair_a, pair_b),
                xy=(x_vals[idx], y_vals[idx]),
                xytext=(7 if left_half else -7, 6 if lower_half else -6),
                textcoords="offset points",
                ha="left" if left_half else "right",
                va="bottom" if lower_half else "top",
                fontsize=ANNOTATION_FONT_SIZE - 1, color=".2", zorder=4,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none",
                          alpha=0.8),
            )

    # Pearson r (linear agreement) and Spearman ρ (same pair ordering)
    if len(common) >= 2:
        r = np.corrcoef(x_vals, y_vals)[0, 1]
        rho = np.corrcoef(
            pd.Series(x_vals).rank().values, pd.Series(y_vals).rank().values,
        )[0, 1]
        stats_text = (
            f"Pearson r = {r:.2f}" if not np.isnan(r) else "Pearson r = n/a"
        ) + "\n" + (
            f"Spearman ρ = {rho:.2f}" if not np.isnan(rho) else "Spearman ρ = n/a"
        )
    else:
        stats_text = "r = n/a (< 2 points)"
    stats_text += f"\nn = {len(common)} pairs"
    if few_pairs:
        stats_text += (
            f"\nWarning: fewer than {MIN_SCATTER_PAIRS} pairs,\n"
            "correlations are unstable"
        )
    ax.annotate(
        stats_text,
        xy=(0.04, 0.96), xycoords="axes fraction",
        fontsize=ANNOTATION_FONT_SIZE, va="top",
        bbox=dict(
            boxstyle="round,pad=0.4", fc="white", alpha=0.9,
            ec=OKABE_ITO["vermillion"] if few_pairs else ".8",
        ),
    )

    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(f"{label_x} — mean {metric_name} per pair")
    ax.set_ylabel(f"{label_y} — mean {metric_name} per pair")
    ax.set_title(
        f"Model agreement: {label_x} vs {label_y} ({metric_name})",
        pad=TITLE_PAD,
    )
    ax.legend(
        loc="lower right",
        frameon=True,
        fancybox=True,
        framealpha=0.9,
        fontsize=LEGEND_FONT_SIZE,
        title_fontsize=LEGEND_TITLE_SIZE,
    )
    ax.grid(True, color=".9", linewidth=0.6)
    ax.set_axisbelow(True)
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
    # HADDOCK's caprieval values below are measured against the reference
    # HADDOCK was given.  PPInsight stages caprieval without one, and then
    # caprieval compares every model to HADDOCK's own best-scoring model —
    # which is why that model shows DockQ = 1.0.  They are not measured
    # against an experimental native; the quality_* metrics further down are.
    "dockq": {
        "higher_is_better": True,
        "display_name": "DockQ (HADDOCK caprieval)",
        "bounds": (0.0, 1.0),
        "reference_based": True,
        "description": (
            "HADDOCK's own DockQ from its caprieval step (0–1), measured "
            "against the reference HADDOCK was given — without one, against "
            "HADDOCK's own best-scoring model, so that model shows 1.0.  "
            "Not measured against your native complex: for that, run "
            "'ppinsight quality' and use quality_dockq."
        ),
    },
    "irmsd": {
        "higher_is_better": False,
        "display_name": "iRMSD (HADDOCK caprieval)",
        "reference_based": True,
        "description": (
            "HADDOCK's own interface RMSD (Å) from caprieval, against "
            "HADDOCK's reference, not your native.  Lower = better.  "
            "Native-referenced value: quality_irmsd."
        ),
    },
    "lrmsd": {
        "higher_is_better": False,
        "display_name": "lRMSD (HADDOCK caprieval)",
        "reference_based": True,
        "description": (
            "HADDOCK's own ligand RMSD (Å) from caprieval, against HADDOCK's "
            "reference, not your native.  Lower = better.  "
            "Native-referenced value: quality_lrmsd."
        ),
    },
    "fnat": {
        "higher_is_better": True,
        "display_name": "Fnat (HADDOCK caprieval)",
        "bounds": (0.0, 1.0),
        "reference_based": True,
        "description": (
            "HADDOCK's own fraction of reference contacts from caprieval, "
            "against HADDOCK's reference, not your native.  Higher = better.  "
            "Native-referenced value: quality_fnat."
        ),
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
        "reference_based": True,
        "description": "Rosetta interface RMSD (Å). Lower = better.",
    },
    "rms": {
        "higher_is_better": False,
        "display_name": "RMS",
        "reference_based": True,
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
    # ── Reference-free consensus ranking ────────────────────────────────
    "consrank_score": {
        "higher_is_better": True,
        "display_name": "CONSRANK Score",
        "description": (
            "Contact-map consensus score from Iter-CONSRANK. Higher means "
            "the pose agrees more closely with the ranked pose pool; it is "
            "not a native-structure accuracy score."
        ),
    },
    # ── DockQ quality metrics (from ppinsight.quality) ─────────────────
    # These use a ``quality_`` prefix to distinguish them from the
    # engine-native metrics above (e.g. HADDOCK's own ``dockq`` score).
    # They are produced by ``add_quality_to_scores()`` which compares
    # predicted complexes against a known native structure via DockQ.
    "quality_dockq": {
        "reference_based": True,
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
        "reference_based": True,
        "higher_is_better": True,
        "display_name": "Fnat",
        "bounds": (0.0, 1.0),
        "description": (
            "Fraction of native contacts (Fnat) recovered in the docked "
            "model.  Higher = better."
        ),
    },
    "quality_irmsd": {
        "reference_based": True,
        "higher_is_better": False,
        "display_name": "iRMSD",
        "description": (
            "Interface RMSD (Å) between docked model and native structure.  "
            "Lower = better."
        ),
    },
    "quality_lrmsd": {
        "reference_based": True,
        "higher_is_better": False,
        "display_name": "lRMSD",
        "description": (
            "Ligand RMSD (Å) between docked model and native structure.  "
            "Lower = better."
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


def _require_reference_free(metric: str, purpose: str) -> None:
    """Refuse binder / non-binder analyses on reference-based metrics."""
    if is_reference_based(metric):
        raise ValueError(
            f"{purpose} cannot use '{metric}': it compares each pose with a "
            "reference native structure, and a pair that does not bind has no "
            "native structure, so the value is meaningless for half of the "
            "comparison.  Use a score every pair has, such as an engine's "
            "own score."
        )


def is_reference_based(metric: str) -> bool:
    """True when *metric* compares a pose with a reference structure.

    These values (DockQ, Fnat, RMSDs) only exist for a pair that has a
    reference complex, so they cannot score a pair that does not bind.
    """
    meta = METRIC_METADATA.get(metric.lower())
    return bool(meta and meta.get("reference_based"))


HADDOCK_CAPRIEVAL_METRICS = ("dockq", "fnat", "irmsd", "lrmsd")


def haddock_caprieval_note(metrics) -> str | None:
    """Explain HADDOCK's own caprieval metrics when any appear in *metrics*."""
    present = sorted({
        str(m).lower() for m in metrics
        if str(m).lower() in HADDOCK_CAPRIEVAL_METRICS
    })
    if not present:
        return None
    return (
        f"NOTE: {', '.join(present)} {'is' if len(present) == 1 else 'are'} "
        "HADDOCK's own caprieval value(s), measured against the reference "
        "HADDOCK was given.  PPInsight runs HADDOCK without one, so caprieval "
        "compares each model to HADDOCK's own best-scoring model (which "
        "therefore shows DockQ = 1.0).  They are not measured against your "
        "native complex.  For DockQ against a native structure, run "
        "'ppinsight quality <scores.tsv> <native.pdb> -o <merged.tsv>' and "
        "use the quality_* metrics (e.g. --metric quality_dockq)."
    )


def _apply_metric_bounds(ax, metric: str, axis: str = "y") -> None:
    """Pin a bounded metric's axis to its full range (DockQ: 0 to 1)."""
    bounds = get_metric_bounds(metric)
    if bounds is None:
        return
    if axis == "y":
        ax.set_ylim(bounds)
    else:
        ax.set_xlim(bounds)


def _engines_title(subject: str, models, pair_tuples) -> str:
    """'<subject> by engine', naming the engine when there is only one."""
    models = list(models)
    if len(models) == 1:
        title = f"{subject} — {_engine_label(models[0])}"
    else:
        title = f"{subject} by engine"
    if len(pair_tuples) == 1:
        title += f" ({_pair_label(*pair_tuples[0])})"
    return title


def _rare_good_tier(value_arrays, metric: str, limit: int):
    """CAPRI tier whose poses are rare enough to always draw, if any.

    For DockQ-scale metrics, returns the lowest ``(threshold, name)`` tier
    that has at most a fifth of *limit* poses in every engine with more than
    *limit* poses — acceptable first, then medium, then high.  ``None`` when
    no engine needs sampling, the metric is not on the DockQ scale, or even
    high-quality poses are plentiful (the sample then shows them anyway).
    """
    import numpy as np

    if metric.lower() not in DOCKQ_METRICS:
        return None
    large = [np.asarray(v, dtype=float) for v in value_arrays if len(v) > limit]
    if not large:
        return None
    for threshold, name, _ in DOCKQ_CAPRI_THRESHOLDS:
        if all((v >= threshold).sum() <= limit // 5 for v in large):
            return threshold, name
    return None


def _draw_selection(values, metric: str, limit: int, tier=None):
    """Choose which poses to draw when an engine has more than *limit*.

    Returns ``(background, highlight, sampled)``.  *sampled* is the size of
    a fixed-seed uniform random sample, so the drawn density matches the
    real one.  *highlight* holds rare good poses that a sample would likely
    miss: every pose in *tier* or better for DockQ-scale metrics (see
    :func:`_rare_good_tier`), or the :data:`ALWAYS_DRAWN_BEST` best poses for
    other metrics.  They are drawn only in the highlight style, so
    *background* is the sample minus *highlight*.  At or below *limit*,
    every pose is in *background*.
    """
    import numpy as np

    values = np.asarray(values, dtype=float)
    empty = np.array([], dtype=int)
    if len(values) <= limit:
        return np.arange(len(values)), empty, len(values)
    sample = np.sort(
        np.random.default_rng(0).choice(len(values), size=limit, replace=False)
    )
    if metric.lower() in DOCKQ_METRICS:
        if tier is None:
            return sample, empty, limit
        tail = np.flatnonzero(values >= tier[0])
    else:
        order = np.argsort(values)
        if get_metric_direction(metric):
            order = order[::-1]
        tail = np.sort(order[:ALWAYS_DRAWN_BEST])
    return np.setdiff1d(sample, tail), tail, limit


def _highlight_legend_handle(metric: str, tier=None) -> Line2D:
    """Legend entry for the always-drawn rare good poses of large runs."""
    if tier is not None:
        text = f"Every pose ≥ CAPRI {tier[1]} ({tier[0]:.2f}), always drawn"
    else:
        text = f"Each engine's {ALWAYS_DRAWN_BEST} best poses, always drawn"
    return Line2D([0], [0], marker="o", linestyle="None", markersize=6,
                  markerfacecolor=".45", markeredgecolor="black",
                  markeredgewidth=0.8, label=text)


def _draw_deterministic_stripplot(
    ax,
    data: pd.DataFrame,
    order: list[str],
    metric: str = "",
    color: str = ".15",
) -> tuple[dict[str, int], Line2D | None]:
    """Draw a deterministic strip overlay on top of violin bodies.

    Points sit above the violin fill so they stay visible; opacity scales
    with the largest per-engine pose count so small runs show every pose
    and large ensembles do not turn into a solid block.  Above
    :data:`MAX_DRAWN_POSES` per engine a uniform random sample of dots is
    drawn, with rare good poses drawn in full on top (see
    :func:`_draw_selection`); the violin itself still uses every pose.

    Returns the number of sampled dots per engine and, when rare good poses
    were overlaid, the legend entry describing them (else ``None``).
    """
    import numpy as np

    sampled_parts = []
    highlight_parts = []
    shown: dict[str, int] = {}
    tier = _rare_good_tier(
        [g.to_numpy(float) for _, g in data.groupby("model")["score_value"]],
        metric, MAX_DRAWN_POSES,
    )
    for model, group in data.groupby("model"):
        background, highlight, n_sampled = _draw_selection(
            group["score_value"].to_numpy(), metric, MAX_DRAWN_POSES, tier,
        )
        sampled_parts.append(group.iloc[background])
        highlight_parts.append(group.iloc[highlight])
        shown[model] = n_sampled
    sampled = pd.concat(sampled_parts)
    highlights = pd.concat(highlight_parts)

    rng_state = np.random.get_state()
    np.random.seed(0)
    try:
        sns.stripplot(
            data=sampled,
            x="model",
            y="score_value",
            size=2.2,
            alpha=_point_alpha(max(shown.values())),
            jitter=0.12,
            ax=ax,
            color=color,
            order=order,
            zorder=3,
            rasterized=max(shown.values()) > 2000,
        )
        if not highlights.empty:
            # Rare good poses from a large run: drawn in full on top.
            sns.stripplot(
                data=highlights, x="model", y="score_value", size=3.5,
                alpha=0.9, jitter=0.12, ax=ax, color=color, order=order,
                zorder=4,
            )
    finally:
        np.random.set_state(rng_state)
    handle = None if highlights.empty else _highlight_legend_handle(metric, tier)
    return shown, handle


def _set_engine_ticks(
    ax, data: pd.DataFrame, order: list[str], shown: dict | None = None,
) -> None:
    """Label categorical engine ticks with display names and pose counts."""
    counts = data.groupby("model")["score_value"].count()
    shown = shown or {}
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([
        _engine_tick_label(m, int(counts.get(m, 0)), shown.get(m))
        for m in order
    ])


def _text_scale(fig: plt.Figure) -> float:
    """Text enlargement for wide multi-panel figures.

    A 3-panel figure is shown at the same on-screen width as a 1-panel one,
    so its text would look about half the size without this.
    """
    return float(min(1.5, max(1.0, (fig.get_figwidth() / 8.0) ** 0.5)))


def _finish_figure(
    fig: plt.Figure,
    axes,
    *,
    title: str,
    xlabel: str | None = None,
    legend_handles: list | None = None,
    legend_ncol: int = 3,
    legend_title: str | None = None,
    figure_title: bool = False,
) -> None:
    """Title the figure and stack the x-label and legend below the panels.

    Multi-panel figures get their text scaled up (see :func:`_text_scale`)
    and a figure-level title (also used when *figure_title* is True, for
    single panels whose own title holds other text).  Positions below the
    panels are measured from the drawn tick labels, so the x-label and
    legend never overlap them.
    """
    axes = list(axes)
    scale = _text_scale(fig) if len(axes) > 1 else 1.0
    multi = len(axes) > 1 or figure_title
    if scale != 1.0:
        for ax in axes:
            ax.tick_params(axis="both", labelsize=TICK_LABEL_SIZE * scale)
            for text in (ax.title, ax.xaxis.label, ax.yaxis.label, *ax.texts):
                text.set_fontsize(text.get_fontsize() * scale)
            legend = ax.get_legend()
            if legend is not None:
                for text in (*legend.get_texts(), legend.get_title()):
                    text.set_fontsize(text.get_fontsize() * scale)

    fig_h = fig.get_figheight()
    gap = 0.1 / fig_h  # 0.1 inch, as a figure fraction
    if multi:
        title_h = (PLOT_TITLE_SIZE * scale / 72.0 * 1.8) / fig_h
        fig.tight_layout(pad=1.0, rect=(0, 0, 1, 1 - title_h))
        fig.suptitle(title, fontsize=PLOT_TITLE_SIZE * scale,
                     y=1 - title_h / 2, va="center")
    else:
        axes[0].set_title(title, fontsize=PLOT_TITLE_SIZE, pad=TITLE_PAD)
        if xlabel is not None:
            axes[0].set_xlabel(xlabel)
        fig.tight_layout(pad=1.0)

    renderer = fig.canvas.get_renderer()
    height_px = fig.bbox.height
    bottom = min(ax.get_tightbbox(renderer).y0 for ax in axes) / height_px
    if multi and xlabel is not None:
        label = fig.text(0.5, bottom - gap, xlabel, ha="center", va="top",
                         fontsize=AXIS_LABEL_SIZE * scale)
        bottom = label.get_window_extent(renderer).y0 / height_px
    if legend_handles:
        fig.legend(
            handles=legend_handles,
            loc="upper center",
            bbox_to_anchor=(0.5, bottom - gap),
            ncol=min(legend_ncol, len(legend_handles)),
            frameon=True,
            fancybox=True,
            framealpha=0.9,
            fontsize=BELOW_LEGEND_FONT_SIZE * scale,
            title=legend_title,
            title_fontsize=BELOW_LEGEND_TITLE_SIZE * scale,
            handlelength=2.2,
        )


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
    """Build one CAPRI classification row per pose, from DockQ.

    Each pose is placed in a tier by its DockQ (see
    :data:`DOCKQ_CAPRI_THRESHOLDS`).  The native-referenced
    ``quality_dockq`` is used when the pose has one; HADDOCK's own caprieval
    ``dockq`` is used only when it does not.  Uses explicit pose metadata
    when present and falls back to row-order reconstruction for older
    unified score files that lack ``pose_id``.
    """
    df = scores_df.copy()
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
    df["score_value"] = pd.to_numeric(df["score_value"], errors="coerce")

    if "score_type" not in df.columns:
        return pd.DataFrame(columns=["model", "capri_quality"])

    df = df[df["score_type"].str.lower().isin(DOCKQ_METRICS)].dropna(
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

    import math

    def _nan_to_none(val):
        if val is not None and isinstance(val, float) and math.isnan(val):
            return None
        return val

    def _row_quality(row):
        # Native-referenced quality_dockq wins; HADDOCK's own caprieval
        # dockq is used only when no native-referenced value exists for this
        # pose.  (Checking for the column alone would let HADDOCK's column
        # hide quality_dockq values for every other engine.)
        native = _nan_to_none(row.get("quality_dockq"))
        dockq = native if native is not None else _nan_to_none(row.get("dockq"))
        return capri_quality(dockq=dockq)

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
        Unified scores (long format) with DockQ rows (``quality_dockq``, or
        HADDOCK's own ``dockq``).

    Returns
    -------
    DataFrame
        Columns: ``model``, ``high``, ``medium``, ``acceptable``,
        ``incorrect``, ``total``, ``pct_acceptable_plus``.
    """
    pose_table = _capri_pose_table(scores_df)
    if pose_table.empty or "capri_quality" not in pose_table.columns:
        raise ValueError(_NO_DOCKQ_FOR_TIERS)

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


_NO_DOCKQ_FOR_TIERS = (
    "CAPRI tiers are assigned from DockQ, but the scores have no DockQ rows "
    "(score_type 'quality_dockq', or HADDOCK's own 'dockq')."
)

CAPRI_TIER_LABELS = {
    "incorrect": "Incorrect (DockQ < 0.23)",
    "acceptable": "Acceptable (0.23 ≤ DockQ < 0.49)",
    "medium": "Medium (0.49 ≤ DockQ < 0.80)",
    "high": "High (DockQ ≥ 0.80)",
}


def _capri_tier_source(classified: pd.DataFrame) -> str:
    """Legend title naming where the DockQ values behind the tiers came from."""
    def _present(name: str) -> bool:
        return name in classified.columns and classified[name].notna().any()

    if _present("quality_dockq") and _present("dockq"):
        source = "DockQ vs the native where available, else HADDOCK's own"
    elif _present("dockq"):
        source = "HADDOCK's own caprieval DockQ, not vs the native"
    else:
        source = "DockQ measured against the native structure"
    return f"CAPRI tier ({source})"


def quality_bar_chart(
    scores_df: pd.DataFrame,
    plot_title: str | None = None,
    output: str | None = None,
) -> plt.Figure:
    """100 % stacked bar chart: percentage of poses in each CAPRI tier.

    Each engine gets one bar that fills from 0 % to 100 %.  The bar is
    divided into four colour-coded tiers (incorrect, acceptable, medium,
    high), assigned from each pose's DockQ with the thresholds in the
    legend.  The pose count (``n = X``) is shown under each engine name so
    the reader knows how many individual poses the percentages are based on.

    The plot covers exactly one protein pair: poses from different targets
    are never pooled into one percentage.

    Inspired by the cumulative bar charts in CAPRI assessments (Lensink
    MF *et al.*, *Proteins* 84(S1):323-348, 2016).

    Parameters
    ----------
    scores_df : DataFrame
        Unified scores for one protein pair with DockQ rows
        (``quality_dockq``, or HADDOCK's own ``dockq``).  Other metrics are
        ignored.
    plot_title : str or None
    output : str or None

    Returns
    -------
    matplotlib.figure.Figure

    Raises
    ------
    ValueError
        If *scores_df* has no DockQ rows or covers more than one protein pair.
    """
    import numpy as np

    lowered = scores_df.copy()
    lowered.columns = [c.strip().lower() for c in lowered.columns]
    _, _, pair_tuples = _detect_pairs(lowered)
    if len(pair_tuples) > 1:
        raise ValueError(
            "quality_bar shows one target at a time, but the scores contain "
            f"{len(pair_tuples)} protein pairs."
        )

    # ── Classify per pose (from DockQ) ─────────────────────────────────
    classified = _capri_pose_table(scores_df)
    if classified.empty or "capri_quality" not in classified.columns:
        raise ValueError(_NO_DOCKQ_FOR_TIERS)

    models = sorted(classified["model"].unique())
    n = len(models)

    tier_colours = CAPRI_TIER_COLORS
    tier_order = ["incorrect", "acceptable", "medium", "high"]

    fig_w = max(6, n * 2.2)
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

    # Pose count sits under each engine name so it never crowds the bars.
    ax.set_xticks(x)
    ax.set_xticklabels(
        [_engine_tick_label(m, total_table[m]) for m in models],
        rotation=30 if n > 4 else 0,
        ha="right" if n > 4 else "center",
    )
    ax.set_ylim(0, 100)
    ax.set_ylabel("% of poses")
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0f}%"))
    sns.despine(ax=ax)

    # The legend states each tier's DockQ range and where DockQ came from.
    for handle, tier in zip(handles, tier_order, strict=False):
        handle.set_label(CAPRI_TIER_LABELS[tier])
    if len(models) == 1:
        title = f"Percentage of {_engine_label(models[0])} poses in each CAPRI tier"
    else:
        title = "Percentage of poses in each CAPRI tier, by engine"
    if pair_tuples:
        title += f" ({_pair_label(*pair_tuples[0])})"
    _finish_figure(
        fig, [ax],
        title=plot_title or title,
        xlabel="Docking Engine",
        legend_handles=handles,
        legend_ncol=2,
        legend_title=_capri_tier_source(classified),
    )

    if output:
        fig.savefig(output, bbox_inches="tight", dpi=150)
        print(f"Plot saved to {output}")
    else:
        plt.show()
    return fig


# ---------------------------------------------------------------------------
# Strip plot
# ---------------------------------------------------------------------------

def strip_plot(
    scores_df: pd.DataFrame,
    metric: str,
    plot_title: str | None = None,
    output: str | None = None,
) -> plt.Figure:
    """Strip plot: one dot per pose, with a box summary per engine.

    Nothing is smoothed, so isolated near-native poses, gaps, and the true
    data range are all visible.  The box shows the median and interquartile
    range; whiskers reach 1.5 × IQR.  Pose counts sit under each engine name.
    CAPRI tier thresholds are drawn when *metric* is on the DockQ scale, and
    bounded metrics always span their full range (DockQ: 0 to 1).

    Very large runs (more than :data:`MAX_DRAWN_POSES` poses for an engine)
    draw a fixed-seed uniform random sample of that size, so the dot
    density stays true, and overlay every rare good pose in an outlined
    style (see :func:`_draw_selection`).  The box and the pose count still
    use every pose, and the tick label says how many dots were sampled.

    When the data contains **≥ 2 protein pairs**, each pair gets its own
    side-by-side panel so poses from different targets are never pooled.

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
    if "score_type" in df.columns:
        df = df[df["score_type"].str.lower() == metric.lower()]
    df["score_value"] = pd.to_numeric(df["score_value"], errors="coerce")
    df = df.dropna(subset=["score_value"])
    if df.empty:
        raise ValueError(f"No rows match metric='{metric}'")

    models = sorted(df["model"].unique())
    metric_name = get_metric_display_name(metric)
    model_palette = _model_palette(models)
    pa, pb, pair_tuples = _detect_pairs(df)

    if len(pair_tuples) >= 2:
        panels = [
            (_pair_label(a, b), df[(df[pa] == a) & (df[pb] == b)])
            for a, b in pair_tuples
        ]
    else:
        panels = [(None, df)]

    # A lone panel needs room for the two-column legend underneath.
    panel_w = max(4.5 if len(panels) > 1 else 6.5, 1.9 * len(models) + 1.0)
    fig, axes = plt.subplots(
        1, len(panels), figsize=(panel_w * len(panels), 5.5),
        sharey=True, squeeze=False,
    )
    axes = list(axes[0])
    rng = np.random.default_rng(0)
    capri_drawn = False
    highlighted = False
    tier = _rare_good_tier(
        [g.to_numpy(float) for _, data in panels
         for _, g in data.groupby("model")["score_value"]],
        metric, MAX_DRAWN_POSES,
    )

    for idx, (ax, (panel_title, data)) in enumerate(zip(axes, panels, strict=False)):
        panel_models = [m for m in models if m in set(data["model"])]
        positions = np.arange(len(panel_models))
        box_data = []
        tick_labels = []
        for pos, model in zip(positions, panel_models, strict=False):
            vals = data.loc[data["model"] == model, "score_value"].to_numpy(float)
            background, highlight, sampled = _draw_selection(
                vals, metric, MAX_DRAWN_POSES, tier,
            )
            drawn = vals[background]
            ax.scatter(
                pos + rng.uniform(-0.18, 0.18, len(drawn)), drawn,
                s=_point_size(len(drawn)),
                color=model_palette[model],
                alpha=_point_alpha(len(drawn)),
                edgecolors="none",
                zorder=2,
                clip_on=False,
                rasterized=len(drawn) > 2000,
            )
            if len(highlight):
                highlighted = True
                ax.scatter(
                    pos + rng.uniform(-0.18, 0.18, len(highlight)), vals[highlight],
                    s=18, color=model_palette[model], alpha=0.9,
                    edgecolors="black", linewidths=0.5, zorder=4, clip_on=False,
                )
            box_data.append(vals)
            tick_labels.append(_engine_tick_label(model, len(vals), sampled))
        ax.boxplot(
            box_data,
            positions=positions,
            widths=0.5,
            showfliers=False,
            patch_artist=True,
            boxprops={"facecolor": "none", "edgecolor": "#333333", "linewidth": 1.2},
            medianprops={"color": "#111111", "linewidth": 2.0},
            whiskerprops={"color": "#333333", "linewidth": 1.1},
            capprops={"color": "#333333", "linewidth": 1.1},
            zorder=3,
        )
        capri_drawn = _draw_capri_lines(ax, metric) or capri_drawn
        ax.set_xticks(positions)
        ax.set_xticklabels(tick_labels)
        ax.set_xlim(-0.6, len(panel_models) - 0.4)
        _apply_metric_bounds(ax, metric)
        if panel_title:
            ax.set_title(panel_title, fontsize=FACET_TITLE_SIZE)
        ax.set_xlabel("")
        ax.set_ylabel(metric_name if idx == 0 else "")
        _add_faint_horizontal_gridlines(ax)
        sns.despine(ax=ax)

    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="None", markersize=6,
               markerfacecolor=".45", markeredgecolor="none",
               label="One dot = one pose"),
        Line2D([0], [0], color="#111111", lw=2.0,
               label="Box: median and IQR (whiskers 1.5 × IQR)"),
    ]
    if highlighted:
        legend_handles.append(_highlight_legend_handle(metric, tier))
    if capri_drawn:
        legend_handles += _capri_legend_handles()

    _finish_figure(
        fig, axes,
        title=plot_title or _engines_title(metric_name, models, pair_tuples),
        xlabel="Docking Engine",
        legend_handles=legend_handles,
        legend_ncol=2,
    )

    if output:
        fig.savefig(output, bbox_inches="tight", dpi=180)
        print(f"Plot saved to {output}")
    else:
        plt.show()
    return fig


# ---------------------------------------------------------------------------
# Funnel plot (native-reference rank audit only)
# ---------------------------------------------------------------------------

FUNNEL_DOCKQ_METRIC = "quality_dockq"


def _funnel_inputs(scores_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split a unified scores table into engine scores and native DockQ.

    The funnel plot only ever pairs an engine's own score with
    native-referenced DockQ (``quality_dockq``) for the same ``pose_id``.
    Raises ``ValueError`` with an actionable message when that is not
    possible.
    """
    df = scores_df.copy()
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
    for required in ("model", "score_type", "score_value", "pose_id"):
        if required not in df.columns:
            raise ValueError(
                f"Funnel plot needs a '{required}' column so each engine "
                "score can be matched to the DockQ value of the same pose."
            )
    df["score_value"] = pd.to_numeric(df["score_value"], errors="coerce")
    score_types = df["score_type"].astype(str).str.lower()
    is_dockq = score_types == FUNNEL_DOCKQ_METRIC
    if not is_dockq.any():
        raise ValueError(
            "Funnel plot needs native-referenced DockQ rows (score_type "
            f"'{FUNNEL_DOCKQ_METRIC}') in the same file as the engine scores."
        )

    _, _, pair_tuples = _detect_pairs(df)
    if len(pair_tuples) > 1:
        raise ValueError(
            "Funnel plot compares the poses of one target, but the scores "
            f"contain {len(pair_tuples)} protein pairs."
        )

    quality = df[is_dockq].rename(columns={"score_value": "DockQ"})
    quality = quality[["model", "pose_id", "DockQ"]]
    engine_scores = df[~score_types.str.startswith("quality_")]
    return engine_scores, quality


def funnel_plot(
    scores_df: pd.DataFrame,
    plot_title: str | None = None,
    output: str | None = None,
) -> plt.Figure:
    """Funnel plot for the native-reference rank audit — nothing else.

    One panel per engine.  Each point is one pose: x is its DockQ against
    the native structure (``quality_dockq``, written by ``ppinsight
    quality``) and y is the **engine's own** primary ranking score for that
    same pose (HADDOCK ``score``, LightDock ``luciferin_score``, Rosetta
    ``i_sc``; other engines use their first native metric).  Scores and
    DockQ are joined on ``pose_id`` only, so pose identity is preserved.

    The y-axis is oriented so the engine's preferred direction is always
    down.  A scoring function that recognises near-native poses produces a
    funnel: the best-scored poses (bottom) sit at high DockQ (right).  Each
    panel reports how many poses were evaluated, the best DockQ, and that
    pose's rank by the engine's own score.

    This plot deliberately takes no metric argument.  It is not a general
    metric-versus-metric scatter: it always pairs an engine's own score
    with native-referenced DockQ for the same poses, because that is the
    only question it answers ("does the engine's score put its near-native
    poses at the top?").

    Parameters
    ----------
    scores_df : DataFrame
        Unified scores for **one** protein pair containing both engine-native
        rows and ``quality_dockq`` rows that share ``pose_id`` values (the
        output of ``ppinsight quality <scores.tsv> <native.pdb> -o …``).
    plot_title : str or None
    output : str or None

    Returns
    -------
    matplotlib.figure.Figure

    Raises
    ------
    ValueError
        If DockQ rows are missing, the scores span several pairs, or no
        pose IDs match between engine scores and DockQ rows.
    """
    import warnings as _warnings

    from ppinsight.quality import rank_engine_poses_by_dockq

    engine_scores, quality = _funnel_inputs(scores_df)
    ranked = rank_engine_poses_by_dockq(engine_scores, quality)
    if ranked.empty:
        raise ValueError(
            "Funnel plot found no pose_id present in both the engine scores "
            f"and the '{FUNNEL_DOCKQ_METRIC}' rows, so no pose could be matched."
        )

    models = sorted(ranked["model"].unique())
    unmatched = sorted(set(engine_scores["model"].unique()) - set(models))
    if unmatched:
        _warnings.warn(
            "No DockQ values matched the poses of: "
            + ", ".join(_engine_label(m) for m in unmatched)
            + ".  These engines are left out of the funnel plot.",
            stacklevel=2,
        )

    palette = _model_palette(models)
    fig, axes = plt.subplots(
        1, len(models), figsize=(4.8 * len(models), 5.2),
        sharex=True, squeeze=False,
    )
    axes = axes[0]

    for ax, model in zip(axes, models, strict=False):
        sub = ranked[ranked["model"] == model]
        score_type = str(sub["native_score_type"].iloc[0])
        _draw_capri_lines(ax, FUNNEL_DOCKQ_METRIC, orientation="v")
        ax.scatter(
            sub["DockQ"], sub["native_score"],
            s=_point_size(len(sub)),
            color=palette[model],
            alpha=max(0.35, _point_alpha(len(sub))),
            edgecolors="white", linewidths=0.4,
            zorder=2, clip_on=False, rasterized=len(sub) > 2000,
        )
        top = sub.loc[sub["engine_rank"].idxmin()]
        best = sub.loc[sub["DockQ"].idxmax()]
        ax.scatter(
            [top["DockQ"]], [top["native_score"]], s=110,
            facecolors="none", edgecolors="black", linewidths=1.6,
            zorder=4, clip_on=False,
        )
        ax.scatter(
            [best["DockQ"]], [best["native_score"]], s=150, marker="*",
            color="black", edgecolors="white", linewidths=0.6,
            zorder=5, clip_on=False,
        )

        if get_metric_direction(score_type):
            ax.invert_yaxis()
        ax.set_xlim(0.0, 1.0)
        # Plain words: an arrow glyph would rotate with the vertical label.
        ax.set_ylabel(f"{get_metric_display_name(score_type)} (best at bottom)")
        poses_ranked = int(sub["poses_ranked"].iloc[0])
        ax.set_title(
            f"{_engine_label(model)}\n"
            f"n = {len(sub)} of {poses_ranked} ranked poses\n"
            f"best DockQ {best['DockQ']:.2f} at engine rank "
            f"{int(best['engine_rank'])}",
            fontsize=FACET_TITLE_SIZE - 1,
        )
        ax.grid(True, color=".92", linewidth=0.6)
        ax.set_axisbelow(True)
        sns.despine(ax=ax)

    _, _, pair_tuples = _detect_pairs(engine_scores)
    title = plot_title or "Engine's own score vs DockQ" + (
        f" ({_pair_label(*pair_tuples[0])})" if pair_tuples else ""
    )
    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="None", markersize=9,
               markerfacecolor="none", markeredgecolor="black",
               markeredgewidth=1.6, label="Engine's top-ranked pose"),
        Line2D([0], [0], marker="*", linestyle="None", markersize=12,
               color="black", label="Best-DockQ pose"),
        *_capri_legend_handles(),
    ]
    _finish_figure(
        fig, axes,
        title=title,
        xlabel="DockQ vs native structure",
        legend_handles=legend_handles,
        legend_ncol=5 if len(models) > 1 else 3,
        figure_title=True,
    )

    if output:
        fig.savefig(output, bbox_inches="tight", dpi=180)
        print(f"Plot saved to {output}")
    else:
        plt.show()
    return fig


# ---------------------------------------------------------------------------
# Top-N plot (native-reference rank audit only)
# ---------------------------------------------------------------------------

TOPN_CHECKPOINTS = (1, 5, 10, 20, 50, 100, 200, 500, 1000, 5000, 10000, 50000)


def topn_table(scores_df: pd.DataFrame) -> pd.DataFrame:
    """Best DockQ within each engine's own top-N poses at standard N values.

    Uses the same inputs and pose-ID matching as :func:`funnel_plot`; N runs
    over :data:`TOPN_CHECKPOINTS` up to the deepest rank with a DockQ value.

    Returns
    -------
    DataFrame
        Columns ``model``, ``N``, ``best_dockq``, ``poses_with_dockq``.
    """
    import numpy as np

    from ppinsight.quality import rank_engine_poses_by_dockq

    engine_scores, quality = _funnel_inputs(scores_df)
    ranked = rank_engine_poses_by_dockq(engine_scores, quality)
    if ranked.empty:
        raise ValueError(
            "Top-N plot found no pose_id present in both the engine scores "
            f"and the '{FUNNEL_DOCKQ_METRIC}' rows, so no pose could be matched."
        )
    max_rank = int(ranked["engine_rank"].max())
    checkpoints = [n for n in TOPN_CHECKPOINTS if n <= max_rank]
    if not checkpoints or checkpoints[-1] != max_rank:
        checkpoints.append(max_rank)

    rows = []
    for model, sub in ranked.groupby("model"):
        ranks = sub["engine_rank"].to_numpy()
        dockq = sub["DockQ"].to_numpy()
        for n in checkpoints:
            within = dockq[ranks <= n]
            rows.append({
                "model": model,
                "N": n,
                "best_dockq": round(float(within.max()), 4) if len(within) else np.nan,
                "poses_with_dockq": int(len(within)),
            })
    return pd.DataFrame(rows)


def topn_plot(
    scores_df: pd.DataFrame,
    plot_title: str | None = None,
    output: str | None = None,
) -> plt.Figure:
    """Best DockQ within each engine's own top-N poses — rank audit only.

    For every N, each engine's curve shows the best DockQ (against the
    native structure, ``quality_dockq``) among its N best poses by its
    **own** score.  Comparing engines at the same N removes the bias that
    pose-averaged metrics have when engines output very different numbers
    of poses, and it is the curve behind CAPRI-style "top-1 / top-10
    success".  A curve that reaches a CAPRI line at small N means the
    engine ranks a pose of that quality near the top.

    Inputs and pose-ID matching are the same as :func:`funnel_plot`: one
    protein pair, engine scores plus ``quality_dockq`` rows sharing
    ``pose_id``.  Like the funnel plot, it takes no metric argument.

    Parameters
    ----------
    scores_df : DataFrame
    plot_title : str or None
    output : str or None

    Returns
    -------
    matplotlib.figure.Figure
    """
    import numpy as np
    from matplotlib.ticker import FuncFormatter, NullLocator

    from ppinsight.quality import rank_engine_poses_by_dockq

    engine_scores, quality = _funnel_inputs(scores_df)
    ranked = rank_engine_poses_by_dockq(engine_scores, quality)
    if ranked.empty:
        raise ValueError(
            "Top-N plot found no pose_id present in both the engine scores "
            f"and the '{FUNNEL_DOCKQ_METRIC}' rows, so no pose could be matched."
        )

    models = sorted(ranked["model"].unique())
    palette = _model_palette(models)
    max_rank = int(ranked["engine_rank"].max())

    fig, ax = plt.subplots(figsize=(8, 5.2))
    _draw_capri_lines(ax, FUNNEL_DOCKQ_METRIC, orientation="h", zorder=1)
    handles = []
    for model in models:
        sub = ranked[ranked["model"] == model].sort_values("engine_rank")
        ranks = sub["engine_rank"].to_numpy()
        best = np.maximum.accumulate(sub["DockQ"].to_numpy())
        x = np.append(ranks, ranks[-1])
        y = np.append(best, best[-1])
        line, = ax.step(
            x, y, where="post", color=palette[model], lw=2.4, zorder=3,
            label=(
                f"{_engine_label(model)} (DockQ for {len(sub):,} of "
                f"{int(sub['poses_ranked'].iloc[0]):,} ranked poses)"
            ),
        )
        handles.append(line)

    ax.set_xscale("log")
    ax.set_xlim(1, max(2, max_rank))
    ticks = [n for n in TOPN_CHECKPOINTS if n <= max(2, max_rank)]
    ax.set_xticks(ticks)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("Best DockQ within top N")
    ax.grid(True, color=".92", linewidth=0.6)
    ax.set_axisbelow(True)
    sns.despine(ax=ax)

    _, _, pair_tuples = _detect_pairs(engine_scores)
    title = plot_title or "Best DockQ in each engine's top N poses" + (
        f" ({_pair_label(*pair_tuples[0])})" if pair_tuples else ""
    )
    _finish_figure(
        fig, [ax],
        title=title,
        xlabel="N = number of top-ranked poses, by each engine's own score",
        legend_handles=[*handles, *_capri_legend_handles()],
        legend_ncol=1 if len(models) > 1 else 2,
    )

    if output:
        fig.savefig(output, bbox_inches="tight", dpi=180)
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
    are drawn as vertical dashed lines when *metric* is on the DockQ scale
    (``dockq`` or ``quality_dockq``).

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
        ax.step(vals, cdf, where="post", color=color, lw=2,
                label=f"{_engine_label(model)} (n = {len(vals)})")

    # CAPRI threshold lines for HADDOCK dockq and native-referenced DockQ
    if metric.lower() in {"dockq", "quality_dockq"}:
        for thr, label, color in DOCKQ_CAPRI_THRESHOLDS:
            ax.axvline(thr, color=color, lw=1.3, ls="--", alpha=0.75,
                       label=f"CAPRI {label} ({thr:.2f})")

    if metric_bounds is not None:
        ax.set_xlim(metric_bounds)
    ax.set_xlabel(metric_name)
    ax.set_ylabel("Cumulative fraction of poses")
    ax.set_ylim(0, 1.05)
    _, _, pair_tuples = _detect_pairs(sub)
    if len(pair_tuples) == 1:
        where = f" ({_pair_label(*pair_tuples[0])})"
    elif len(pair_tuples) > 1:
        where = f" ({len(pair_tuples)} protein pairs pooled)"
    else:
        where = ""
    ax.set_title(plot_title or f"Empirical CDF – {metric_name}{where}", pad=TITLE_PAD)
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
    print(_with_engine_labels(tbl).to_string(index=False))
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
    print(_with_engine_labels(tbl).to_string(index=False))
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
    _require_reference_free(metric, "Binder / non-binder classification")
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
    model_labels = [_engine_label(m) for m in models]
    models_str = (
        ", ".join(model_labels) if n_models <= 5
        else ", ".join(model_labels[:5]) + ", …"
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
    haddock_present = [m for m in metrics if m.lower() in HADDOCK_CAPRIEVAL_METRICS]
    if haddock_present:
        print(
            f"             {', '.join(haddock_present)} = HADDOCK's own caprieval "
            "values (vs HADDOCK's reference: its own best model unless one\n"
            "             was supplied), not vs your native.\n"
            "             Native-referenced versions: quality_* from "
            "'ppinsight quality'."
        )
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
        "strip", "✅ yes",
        "Works with any setup — every pose as a point, with a box summary.",
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
            "'ppinsight collect' with --label-file <pairs_file>.",
        ))

    # scatter — needs ≥2 models sharing a metric on the same pairs
    if n_models < 2:
        rows.append((
            "scatter", "❌ no",
            f"Needs ≥ 2 models (you have {n_models}). "
            "Use --models MODEL_A MODEL_B.",
        ))
    elif shared_metrics and len(shared_pairs) >= MIN_SCATTER_PAIRS:
        rows.append((
            "scatter",
            "✅ yes",
            f"Shared metrics: {', '.join(sorted(shared_metrics)[:5])}.  "
            f"Shared pairs: {len(shared_pairs)}.  Use --models to pick two."
        ))
    elif shared_metrics and shared_pairs:
        rows.append((
            "scatter",
            "⚠️  limited",
            f"Only {len(shared_pairs)} shared pair(s); correlations are "
            f"unstable below {MIN_SCATTER_PAIRS}.  The plot still draws, "
            "with a warning.",
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

    # quality_bar — needs DockQ (tiers come from DockQ) and one pair
    has_quality_metrics = False
    if "score_type" in scores_df.columns:
        st_lower = set(scores_df["score_type"].dropna().str.lower().unique())
        has_quality_metrics = bool(st_lower & DOCKQ_METRICS)
    if not has_quality_metrics:
        rows.append((
            "quality_bar", "❌ no",
            "Needs DockQ for CAPRI tiers: run 'ppinsight quality' to add "
            "quality_dockq (or use HADDOCK's own dockq).",
        ))
    elif n_pairs > 1:
        rows.append((
            "quality_bar", "⚠️  limited",
            "One target per plot — plot with --pair proteinA:proteinB.",
        ))
    else:
        rows.append((
            "quality_bar", "✅ yes",
            "% of poses in each CAPRI tier (from DockQ). No --metric needed.",
        ))

    rows.append((
        "cdf", "✅ yes",
        "Works with any numeric metric; DockQ adds CAPRI threshold lines.",
    ))

    # funnel — rank audit: engine's own score vs native DockQ, same poses
    has_native_dockq = (
        "score_type" in scores_df.columns
        and scores_df["score_type"].astype(str).str.lower()
        .eq(FUNNEL_DOCKQ_METRIC).any()
    )
    if not has_native_dockq:
        rows.append((
            "funnel", "❌ no",
            "Needs native-referenced DockQ rows. Run 'ppinsight quality "
            "<scores.tsv> <native.pdb> -o <merged.tsv>' first.",
        ))
    elif "pose_id" not in scores_df.columns:
        rows.append((
            "funnel", "❌ no",
            "Needs a pose_id column to match engine scores to DockQ.",
        ))
    elif n_pairs > 1:
        rows.append((
            "funnel", "⚠️  limited",
            "One target per plot — plot with --pair proteinA:proteinB.",
        ))
    else:
        rows.append((
            "funnel", "✅ yes",
            "Rank audit: each engine's own score vs DockQ for the same "
            "poses. No --metric needed.",
        ))
    # topn shares funnel's input; repeat its verdict with its own wording.
    funnel_status, funnel_detail = rows[-1][1], rows[-1][2]
    rows.append((
        "topn", funnel_status,
        "Rank audit: best DockQ within each engine's top N poses. "
        "No --metric needed." if funnel_status == "✅ yes" else funnel_detail,
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
            "(collect with --label-file)."
        )
        extras.append(
            "--split-label    ❌ needs labels "
            "(collect with --label-file)."
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
            "--capri-quality  ❌ needs DockQ (quality_dockq or dockq)."
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

    if "reference native structure" in lower:
        print(
            "Hint: pick a score that exists without a native structure, e.g. "
            "--metric score (HADDOCK), luciferin_score (LightDock), or i_sc "
            "(Rosetta).  Each of these belongs to one engine, so the result "
            "covers that engine only.  Use --list-metrics to see what your "
            "file contains.",
            file=sys.stderr,
        )
    elif "label" in lower:
        # Missing interaction labels — need to annotate at collection time
        print(
            "Hint: your scores file has no 'label' column.  Re-run "
            "'ppinsight collect' with --label-file <pairs_file> to annotate "
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
            "Re-run 'ppinsight collect' with --pair or --label-file.",
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


def _print_tier_error(exc: Exception) -> None:
    """Print a CAPRI-tier error with the step that fixes it."""
    print(f"ERROR: {exc}", file=sys.stderr)
    msg = str(exc).lower()
    if "protein pairs" in msg:
        hint = (
            "pick one target with --pair proteinA:proteinB "
            "(use --list-pairs to see the names)."
        )
    else:
        hint = (
            "run 'ppinsight quality <scores.tsv> <native.pdb> -o "
            "<merged.tsv>' to add DockQ against your native structure, then "
            "use the merged file."
        )
    print(f"Hint: {hint}", file=sys.stderr)


def _print_haddock_note(metrics) -> None:
    """Explain HADDOCK's own caprieval metrics on stderr when they appear."""
    note = haddock_caprieval_note(metrics)
    if note:
        print(note, file=sys.stderr)


def _run_rank_audit_cli(args, analysis_df: pd.DataFrame, pair_filter) -> None:
    """Print the rank-audit table, then draw the funnel or top-N plot."""
    kind = args.plot_type
    if args.metric:
        print(
            f"ERROR: --plot-type {kind} does not take --metric.  It always "
            "uses each engine's own ranking score and DockQ measured "
            "against the native structure, for the same poses.",
            file=sys.stderr,
        )
        print("Hint: drop --metric and run the command again.", file=sys.stderr)
        sys.exit(2)

    from ppinsight.quality import summarize_engine_ranking

    try:
        engine_scores, quality = _funnel_inputs(analysis_df)
        summary = summarize_engine_ranking(engine_scores, quality)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        msg = str(exc).lower()
        if "protein pairs" in msg:
            hint = (
                "pick one target with --pair proteinA:proteinB "
                "(use --list-pairs to see the names)."
            )
        elif "pose_id" in msg and "column" in msg:
            hint = (
                "re-run 'ppinsight collect' so the scores file records a "
                "pose_id for every pose."
            )
        else:
            hint = (
                "run 'ppinsight quality <scores.tsv> <native.pdb> -o "
                "<merged.tsv>' on the collected scores file and plot the "
                "merged file.  Evaluating a separate folder of copied or "
                "renamed pose files gives DockQ rows whose pose_id no longer "
                "matches the engine scores."
            )
        print(f"Hint: {hint}", file=sys.stderr)
        sys.exit(1)

    if summary.empty:
        print(
            "ERROR: no pose_id is shared between the engine scores and the "
            f"'{FUNNEL_DOCKQ_METRIC}' rows.",
            file=sys.stderr,
        )
        print(
            "Hint: run 'ppinsight quality <scores.tsv> <native.pdb> -o "
            "<merged.tsv>' on the collected scores file so DockQ rows keep "
            "the engine's pose IDs.",
            file=sys.stderr,
        )
        sys.exit(1)

    print()
    print("── Rank audit: each engine's own score vs DockQ (same poses) ──")
    print(_with_engine_labels(summary).to_string(index=False))
    print()
    if kind == "topn":
        table = topn_table(analysis_df)
        wide = table.pivot(index="N", columns="model", values="best_dockq")
        wide.columns = [_engine_label(m) for m in wide.columns]
        print("── Best DockQ within each engine's top N poses ──")
        print(wide.to_string())
        print()

    plot_out = args.output
    if plot_out is None and not args.show:
        plot_out = _default_plot_output_path(
            kind, "engine_score_vs_dockq", pair_filter,
        )
    if plot_out is not None:
        os.makedirs(os.path.dirname(plot_out) or ".", exist_ok=True)
    draw = topn_plot if kind == "topn" else funnel_plot
    with _warnings_to_stderr():
        draw(analysis_df, plot_title=None, output=plot_out)


@contextmanager
def _warnings_to_stderr():
    """Show plot warnings as plain ``WARNING:`` lines instead of tracebacks."""
    import warnings as _warnings

    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter("always")
        yield
    for warning in caught:
        print(f"WARNING: {warning.message}", file=sys.stderr)


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
            "'ppinsight collect --label-file'; --pairs is an alias).  "
            "Lets you visually assess "
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
        choices=["violin", "strip", "roc", "scatter",
                 "quality_bar", "cdf", "funnel", "topn"],
        default="violin",
        help=(
            "Plot style.  'violin' (default): smoothed score distributions "
            "with individual data points — best for large pose ensembles.  "
            "'strip': every pose as a point with a box summary — best when "
            "you want to see each pose (up to a few hundred per engine).  "
            "'roc': ROC curves with AUC — needs interaction / "
            "non-interaction labels and a score every pair has (an "
            "engine's own score, not DockQ, which needs a native "
            "structure).  "
            "'scatter': model-agreement plot — same metric, same pairs, "
            "one model per axis (use with --models).  "
            "'quality_bar': percentage of poses in each CAPRI tier per "
            "engine, with tiers from DockQ — one protein pair at a time.  "
            "'cdf': empirical CDF per engine — CAPRI thresholds shown "
            "for DockQ.  "
            "'funnel': rank audit for one target — each engine's own "
            "score against DockQ vs the native structure, for the same "
            "poses.  Needs the merged file from 'ppinsight quality "
            "<scores.tsv> <native.pdb> -o <merged.tsv>'.  It only ever "
            "plots engine score vs DockQ, so it takes no --metric.  "
            "'topn': rank audit for one target — the best DockQ within "
            "each engine's top N poses by its own score, for every N, so "
            "engines that output different numbers of poses are compared "
            "fairly.  Same input as 'funnel'; takes no --metric."
        ),
    )
    parser.add_argument(
        "--models",
        nargs=2,
        metavar=("MODEL_A", "MODEL_B"),
        default=None,
        help=(
            "Two model names for --plot-type scatter (e.g. --models "
            "haddock rosetta; display names such as HADDOCK3 also work).  "
            "Each protein pair becomes one point: x = model A's mean score "
            "for that pair, y = model B's.  Points along the diagonal mean "
            "the engines agree.  Pearson r and Spearman ρ are annotated.  "
            "With fewer than 10 shared pairs the plot still draws, but a "
            "warning is printed and shown on the figure because the "
            "correlations are unstable."
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
            "Print a CAPRI quality tier summary: count of poses classified "
            "as high / medium / acceptable / incorrect for each engine.  "
            "Tiers come from each pose's DockQ: high ≥ 0.80, medium ≥ 0.49, "
            "acceptable ≥ 0.23 (Basu & Wallner 2016, "
            "DOI: 10.1371/journal.pone.0161879).  DockQ measured against "
            "your native (quality_dockq) is used when present; HADDOCK's own "
            "dockq only when it is not.  Does not require --metric."
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
                metrics = available_metrics(scores_df)
                for m in metrics:
                    print(m)
                _print_haddock_note(metrics)
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
                        "Hint: re-run 'ppinsight collect' with --pair or --label-file.",
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

            if args.metric:
                _print_haddock_note([args.metric])
            elif args.capri_quality or args.plot_type == "quality_bar":
                _print_haddock_note(available_metrics(analysis_df))

            # --capri-quality does NOT need --metric
            if args.capri_quality:
                try:
                    summary = capri_summary_table(analysis_df)
                except ValueError as exc:
                    _print_tier_error(exc)
                    sys.exit(1)
                print(_with_engine_labels(summary).to_string(index=False))
                return

            # funnel and topn never take --metric: they always use each
            # engine's own score and native-referenced DockQ, same poses.
            if args.plot_type in ("funnel", "topn") and not args.table:
                _run_rank_audit_cli(args, analysis_df, pair_filter)
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
                    _print_tier_error(exc)
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

            # Binder / non-binder analyses need a score every pair has;
            # refuse reference-based metrics before printing anything else.
            if args.classify or (args.plot_type == "roc" and not args.table):
                try:
                    _require_reference_free(
                        args.metric,
                        "Binder / non-binder classification" if args.classify
                        else "A ROC curve",
                    )
                except ValueError as exc:
                    _cli_error_with_hint(str(exc), args.metric, scores_df)

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
                print(_with_engine_labels(summary).to_string(index=False))
                return

            # --table: tabular summary only (already printed above)
            if args.table:
                return

            if args.rank:
                try:
                    rnk = ranking_table(analysis_df, metric=args.metric)
                except ValueError as exc:
                    _cli_error_with_hint(str(exc), args.metric, scores_df)
                print(_with_engine_labels(rnk).to_string(index=False))
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
                    with _warnings_to_stderr():
                        fig = model_agreement_scatter(
                            analysis_df,
                            metric=args.metric,
                            model_x=args.models[0],
                            model_y=args.models[1],
                            output=plot_out,
                        )
                    _save_or_show(fig, plot_out)
                    return

                if pt == "strip":
                    plot_out = _resolve_plot_output(pt)
                    fig = strip_plot(
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
