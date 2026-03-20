"""
visualizer – compare docking scores across PPI prediction models.

Supports three input modes (from simplest to most flexible):

1. **Unified scores file** (``scores.tsv`` / ``scores.csv``)
   Columns: ``proteinA, proteinB, model, score_type, score_value[, output_path, timestamp]``
   → :func:`load_scores` → :func:`compare_scores`

2. **Per-model score files** (one CSV/TSV per model, each with a numeric
   column the user wants to compare)
   → :func:`to_plot` → :func:`compare_scores`

3. **Native tool outputs** read directly
   - HADDOCK ``capri_ss.tsv``  → columns like *score*, *dockq*, *irmsd*, *fnat*, …
   - Rosetta ``docking_scores.csv`` → columns *run*, *score*
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

from matplotlib import pyplot as plt
import pandas as pd


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
    for name, df in zip(models, frames):
        if score_type not in df.columns:
            raise LookupError(
                f"The score type '{score_type}' does not exist for the "
                f"'{name}' model. Available columns: {list(df.columns)}"
            )
        col = pd.to_numeric(df[score_type], errors="coerce").dropna()
        means.append(col.mean())
        stds.append(col.std() if len(col) > 1 else 0.0)

    fig, ax = plt.subplots()
    ax.bar(models, means, yerr=stds, capsize=5, edgecolor="black")
    ax.set_xlabel("Interaction model")
    ax.set_ylabel(f"{score_type}")
    ax.set_title(plot_title or f"{score_type} by model")

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

    # Build per-model frames for compare_scores
    grouped = df.groupby("model")
    frames: list[pd.DataFrame] = []
    model_names: list[str] = []
    for name, group in grouped:
        # Create a tiny DataFrame with a column named after the metric
        frame = pd.DataFrame({metric: pd.to_numeric(group["score_value"], errors="coerce")})
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
            if c not in ("model", "proteina", "proteinb", "output_path", "timestamp")
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
    import numpy as np

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
    labels_present = sorted(df["label"].dropna().unique())
    models = sorted(df["model"].unique())

    # Compute means and stds per (model, label)
    grouped = df.groupby(["model", "label"])["score_value"]
    stats = grouped.agg(["mean", "std", "count"]).reset_index()

    n_labels = len(labels_present)
    n_models = len(models)
    x = np.arange(n_models)
    width = 0.35
    offsets = np.linspace(-width * (n_labels - 1) / 2, width * (n_labels - 1) / 2, n_labels)

    # Colour map for labels
    label_colors = {
        "interaction": "#2196F3",       # blue
        "non-interaction": "#FF5722",   # red-orange
        "unknown": "#9E9E9E",           # grey
    }

    fig, ax = plt.subplots(figsize=(max(8, n_models * 1.5), 6))

    for j, lbl in enumerate(labels_present):
        means = []
        stds = []
        counts = []
        for m in models:
            sub = stats[(stats["model"] == m) & (stats["label"] == lbl)]
            if len(sub):
                means.append(sub["mean"].values[0])
                stds.append(sub["std"].values[0] if sub["count"].values[0] > 1 else 0)
                counts.append(int(sub["count"].values[0]))
            else:
                means.append(0)
                stds.append(0)
                counts.append(0)

        color = label_colors.get(lbl, f"C{j}")
        bars = ax.bar(
            x + offsets[j], means, width, yerr=stds,
            capsize=4, label=f"{lbl} (n={sum(counts)})",
            color=color, edgecolor="black", alpha=0.85,
        )

    ax.set_xlabel("Docking Model")
    ax.set_ylabel(metric)
    ax.set_title(plot_title or f"{metric} — interaction vs non-interaction")
    ax.set_xticks(x)
    ax.set_xticklabels(models)
    ax.legend()
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

def violin_plot(
    scores_df: pd.DataFrame,
    metric: str,
    split_by_label: bool = False,
    plot_title: str | None = None,
    output: str | None = None,
) -> plt.Figure:
    """Violin plot showing the full score distribution per model.

    Unlike a bar chart (which only shows the mean), a violin plot reveals
    the **shape** of each model's score distribution — you can see whether
    scores are tightly clustered or spread out, and whether there are
    multiple peaks.

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
    import numpy as np

    df = scores_df.copy()
    if "score_type" in df.columns:
        df = df[df["score_type"].str.lower() == metric.lower()]
    df["score_value"] = pd.to_numeric(df["score_value"], errors="coerce")
    df = df.dropna(subset=["score_value"])

    models = sorted(df["model"].unique())
    has_label = "label" in df.columns and split_by_label

    fig, ax = plt.subplots(figsize=(max(8, len(models) * 2), 6))

    if has_label:
        labels_present = sorted(df["label"].dropna().unique())
        label_colors = {"interaction": "#4CAF50", "non-interaction": "#F44336"}
        positions = np.arange(len(models))

        for i, lbl in enumerate(labels_present):
            data = [
                df[(df["model"] == m) & (df["label"] == lbl)]["score_value"].values
                for m in models
            ]
            # Only plot if there's data
            data = [d if len(d) > 0 else np.array([0.0]) for d in data]
            offset = -0.2 + i * 0.4
            parts = ax.violinplot(
                data, positions=positions + offset,
                showmeans=True, showmedians=True, widths=0.35,
            )
            color = label_colors.get(lbl, f"C{i}")
            for pc in parts["bodies"]:
                pc.set_facecolor(color)
                pc.set_alpha(0.7)
            # Dummy patch for legend
            ax.bar([], [], color=color, label=lbl, alpha=0.7)

        ax.set_xticks(positions)
        ax.set_xticklabels(models)
        ax.legend()
    else:
        data = [
            df[df["model"] == m]["score_value"].values for m in models
        ]
        parts = ax.violinplot(data, showmeans=True, showmedians=True)
        ax.set_xticks(range(1, len(models) + 1))
        ax.set_xticklabels(models)

    ax.set_xlabel("Docking Model")
    ax.set_ylabel(metric)
    ax.set_title(plot_title or f"{metric} distribution by model")
    fig.tight_layout()

    if output:
        fig.savefig(output, bbox_inches="tight", dpi=150)
        print(f"Plot saved to {output}")
    else:
        plt.show()
    return fig


def score_heatmap(
    scores_df: pd.DataFrame,
    metric: str,
    agg: str = "mean",
    plot_title: str | None = None,
    output: str | None = None,
) -> plt.Figure:
    """Heatmap of aggregated scores: one row per protein pair, one column per model.

    This is useful when you have many pairs and want to see at a glance
    which pairs score well across which engines.

    Parameters
    ----------
    scores_df : DataFrame
        Must have ``model``, ``score_type``, ``score_value``, ``proteinA``,
        ``proteinB``.
    metric : str
        Which ``score_type`` to show.
    agg : str
        Aggregation function (``"mean"``, ``"median"``, ``"max"``,
        ``"min"``).  Applied per (pair, model) group.
    plot_title : str | None
    output : str | None

    Returns
    -------
    matplotlib.figure.Figure
    """
    import numpy as np

    df = scores_df.copy()
    if "score_type" in df.columns:
        df = df[df["score_type"].str.lower() == metric.lower()]
    df["score_value"] = pd.to_numeric(df["score_value"], errors="coerce")
    df = df.dropna(subset=["score_value"])

    if "proteinA" not in df.columns or "proteinB" not in df.columns:
        raise ValueError("Need proteinA and proteinB columns for heatmap")

    df["pair"] = df["proteinA"] + " vs " + df["proteinB"]
    pivot = df.pivot_table(
        index="pair", columns="model", values="score_value", aggfunc=agg,
    )

    fig, ax = plt.subplots(figsize=(max(6, len(pivot.columns) * 2),
                                     max(4, len(pivot) * 0.4)))
    im = ax.imshow(pivot.values, aspect="auto", cmap="RdYlGn")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    fig.colorbar(im, ax=ax, label=f"{metric} ({agg})")
    ax.set_title(plot_title or f"{metric} heatmap ({agg})")
    fig.tight_layout()

    if output:
        fig.savefig(output, bbox_inches="tight", dpi=150)
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
        ["model", "proteinA", "proteinB", "label"], as_index=False,
    )["score_value"].mean()
    agg["actual"] = (agg["label"].str.lower() == "interaction").astype(int)

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot([0, 1], [0, 1], "k--", alpha=0.3, label="Random (AUC=0.50)")

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
        ax.plot(fpr_list, tpr_list, label=f"{model} (AUC={auc:.3f})", linewidth=2)

    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title(f"ROC Curve — {metric}")
    ax.legend(loc="lower right")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    fig.tight_layout()

    if output:
        fig.savefig(output, bbox_inches="tight", dpi=150)
        print(f"Plot saved to {output}")
    else:
        plt.show()
    return fig


def rank_comparison_scatter(
    scores_df: pd.DataFrame,
    metric: str,
    model_x: str,
    model_y: str,
    output: str | None = None,
) -> plt.Figure:
    """Scatter plot comparing how two models rank the same protein pairs.

    Each point is one pair.  Points are coloured by label (if available).
    If both models agree on which pairs score high, points will cluster
    along the diagonal.

    Parameters
    ----------
    scores_df : DataFrame
    metric : str
    model_x, model_y : str
        Names of the two models to compare.
    output : str | None
    """
    df = scores_df.copy()
    if "score_type" in df.columns:
        df = df[df["score_type"].str.lower() == metric.lower()]
    df["score_value"] = pd.to_numeric(df["score_value"], errors="coerce")
    df = df.dropna(subset=["score_value"])

    agg = df.groupby(
        ["model", "proteinA", "proteinB"] + (["label"] if "label" in df.columns else []),
        as_index=False,
    )["score_value"].mean()

    dx = agg[agg["model"] == model_x].set_index(["proteinA", "proteinB"])
    dy = agg[agg["model"] == model_y].set_index(["proteinA", "proteinB"])
    common = dx.index.intersection(dy.index)
    if len(common) == 0:
        raise ValueError(f"No common pairs between {model_x} and {model_y}")

    fig, ax = plt.subplots(figsize=(7, 7))
    x_vals = dx.loc[common, "score_value"].values
    y_vals = dy.loc[common, "score_value"].values

    if "label" in dx.columns:
        labels = dx.loc[common, "label"].values
        color_map = {"interaction": "#4CAF50", "non-interaction": "#F44336"}
        colors = [color_map.get(str(l).lower(), "grey") for l in labels]
        for lbl, c in color_map.items():
            ax.scatter([], [], c=c, label=lbl, s=50)
        ax.scatter(x_vals, y_vals, c=colors, alpha=0.7, s=40, edgecolors="black", linewidths=0.5)
        ax.legend()
    else:
        ax.scatter(x_vals, y_vals, alpha=0.7, s=40, edgecolors="black", linewidths=0.5)

    ax.set_xlabel(f"{model_x} — {metric}")
    ax.set_ylabel(f"{model_y} — {metric}")
    ax.set_title(f"Pair-level agreement: {model_x} vs {model_y}")
    fig.tight_layout()

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
        "description": "LightDock luciferin/scoring (DFIRE by default). Higher = better fit.",
    },
    # ── HADDOCK ────────────────────────────────────────────────────────
    "score": {
        "higher_is_better": False,
        "description": "HADDOCK overall score (weighted energy). Lower = better.",
    },
    "dockq": {
        "higher_is_better": True,
        "description": "DockQ quality (0–1). Higher = closer to native complex.",
    },
    "irmsd": {
        "higher_is_better": False,
        "description": "Interface RMSD (Å). Lower = better.",
    },
    "lrmsd": {
        "higher_is_better": False,
        "description": "Ligand RMSD (Å). Lower = better.",
    },
    "fnat": {
        "higher_is_better": True,
        "description": "Fraction of native contacts recovered. Higher = better.",
    },
    # ── Rosetta ────────────────────────────────────────────────────────
    "interface_score": {
        "higher_is_better": False,
        "description": "Rosetta interface energy (REU) from PyRosetta wrapper. Lower = better.",
    },
    "i_sc": {
        "higher_is_better": False,
        "description": (
            "Rosetta interface score (I_sc / dG_separated) — the primary "
            "RosettaDock quality metric.  Lower = stronger binding."
        ),
    },
    "total_score": {
        "higher_is_better": False,
        "description": "Rosetta total energy (REU). Lower = better.",
    },
    "irms": {
        "higher_is_better": False,
        "description": "Rosetta interface RMSD (Å). Lower = better.",
    },
    "rms": {
        "higher_is_better": False,
        "description": "Rosetta ligand RMSD (Å). Lower = better.",
    },
    "dg_separated": {
        "higher_is_better": False,
        "description": (
            "Rosetta dG_separated — binding energy after rigid-body "
            "separation.  Equivalent to I_sc.  Lower = stronger binding."
        ),
    },
    "cluster_size": {
        "higher_is_better": True,
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
        "description": (
            "DockQ score (0–1) from comparison with native structure.  "
            "Higher = better.  Thresholds: acceptable ≥ 0.23, "
            "medium ≥ 0.49, high ≥ 0.80."
        ),
    },
    "quality_fnat": {
        "higher_is_better": True,
        "description": (
            "Fraction of native contacts (Fnat) recovered in the docked "
            "model.  Higher = better."
        ),
    },
    "quality_irmsd": {
        "higher_is_better": False,
        "description": (
            "Interface RMSD (Å) between docked model and native structure.  "
            "Lower = better."
        ),
    },
    "quality_lrmsd": {
        "higher_is_better": False,
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
    agg = df.groupby(["model", "proteinA", "proteinB", "label"], as_index=False)["score_value"].mean()

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
            "Hint: the heatmap plot requires proteinA and proteinB "
            "columns.  Re-run 'ppinsight collect' with --pair or --pairs.",
            file=sys.stderr,
        )
    sys.exit(1)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    """CLI entry-point for score comparison plots.

    Examples::

        compare_scores scores.tsv --metric dockq
        compare_scores model_a.csv model_b.csv --metric score --names HADDOCK Rosetta
    """
    parser = argparse.ArgumentParser(
        prog="compare_scores",
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
        required=True,
        help=(
            "Score type to analyse (e.g. 'dockq', 'score', 'irmsd', "
            "'luciferin_score', 'i_sc').  Must match a value in the "
            "'score_type' column of the scores file.  Use --list-metrics "
            "to see what is available."
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
        "--title", "-t",
        default=None,
        help="Custom plot title (default: auto-generated from metric name).",
    )
    parser.add_argument(
        "--output", "-o",
        default=None,
        help=(
            "Save the plot to a file (png/svg/pdf) instead of showing "
            "interactively.  Useful for batch report generation."
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
        "--split-label",
        action="store_true",
        help=(
            "Split bars/violins by interaction label (interaction vs "
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
        "--normalize-global",
        action="store_true",
        help=(
            "With --normalize: normalise across all models jointly (instead "
            "of per-model, which is the default).  Use global normalisation "
            "when you want an absolute ranking across engines; use per-model "
            "(default) when each engine's internal distribution matters more."
        ),
    )
    parser.add_argument(
        "--plot-type",
        choices=["bar", "violin", "heatmap", "roc", "scatter"],
        default="bar",
        help=(
            "Plot style.  'bar' (default): mean ± std per model — good for "
            "a quick overview.  'violin': full score distributions — reveals "
            "bimodality or outliers.  'heatmap': pair × model grid — useful "
            "when you have many pairs.  'roc': ROC curves with AUC — the "
            "gold standard for binary classification (needs labels).  "
            "'scatter': head-to-head rank comparison of two models (use "
            "with --models)."
        ),
    )
    parser.add_argument(
        "--models",
        nargs=2,
        default=None,
        help=(
            "Two model names for --plot-type scatter (e.g. --models "
            "lightdock haddock).  Each pair is plotted as a point with "
            "one model's score on each axis, revealing agreement or "
            "disagreement between engines."
        ),
    )
    parser.add_argument(
        "--agg",
        choices=["best", "topN_mean", "median", "mean"],
        default=None,
        help=(
            "Aggregate multiple poses per pair before plotting.  Same "
            "strategies as 'ppinsight collect --agg' — direction-aware, so "
            "'best' correctly picks the lowest HADDOCK score or the highest "
            "LightDock luciferin score.  Use 'best' or 'topN_mean' when "
            "comparing top predictions; omit for full-distribution plots."
        ),
    )

    args = parser.parse_args(argv)

    # ── unified file mode (single file) ──────────────────────────
    if len(args.files) == 1:
        try:
            scores_df = load_scores(args.files[0])
        except ValueError:
            # Not a unified file — fall through to per-model mode with one file
            scores_df = None

        if scores_df is not None:
            # Apply aggregation if requested
            if args.agg:
                from ppinsight.collect_scores import aggregate_scores
                scores_df = aggregate_scores(scores_df, strategy=args.agg)

            # Apply normalization if requested.  Direction-aware flipping
            # is always enabled — the pipeline knows which metrics are
            # lower-is-better from METRIC_METADATA.
            if args.normalize:
                scores_df = normalize_scores(
                    scores_df,
                    method=args.normalize,
                    per_model=not args.normalize_global,
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

            if args.classify:
                # Always auto-detect direction from METRIC_METADATA.
                # The pipeline knows which scoring functions are
                # lower-is-better vs higher-is-better.
                try:
                    summary = classification_summary(
                        scores_df,
                        metric=args.metric,
                        threshold=args.threshold,
                    )
                except ValueError as exc:
                    _cli_error_with_hint(str(exc), args.metric, scores_df)
                print(summary.to_string(index=False))
                return

            pair = tuple(args.pair.split(":", 1)) if args.pair else None

            # Dispatch by plot type
            pt = args.plot_type

            try:
                if pt == "violin":
                    violin_plot(
                        scores_df,
                        metric=args.metric,
                        split_by_label=args.split_label,
                        plot_title=args.title,
                        output=args.output,
                    )
                    return

                if pt == "heatmap":
                    score_heatmap(
                        scores_df,
                        metric=args.metric,
                        plot_title=args.title,
                        output=args.output,
                    )
                    return

                if pt == "roc":
                    roc_curve_plot(
                        scores_df,
                        metric=args.metric,
                        output=args.output,
                    )
                    return

                if pt == "scatter":
                    if not args.models or len(args.models) != 2:
                        print(
                            "ERROR: --plot-type scatter requires "
                            "--models MODEL_A MODEL_B",
                            file=sys.stderr,
                        )
                        print(
                            "Hint: pass two model names from the 'model' "
                            "column of your scores file.",
                            file=sys.stderr,
                        )
                        sys.exit(2)
                    rank_comparison_scatter(
                        scores_df,
                        metric=args.metric,
                        model_x=args.models[0],
                        model_y=args.models[1],
                        output=args.output,
                    )
                    return

                # Default: bar chart
                if args.split_label:
                    compare_scores_by_label(
                        scores_df,
                        metric=args.metric,
                        plot_title=args.title,
                        output=args.output,
                    )
                    return

                compare_scores_unified(
                    scores_df,
                    metric=args.metric,
                    pair=pair,
                    plot_title=args.title,
                    output=args.output,
                )
            except ValueError as exc:
                _cli_error_with_hint(str(exc), args.metric, scores_df)
            return

    # ── per-model file mode (multiple files) ─────────────────────
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

    compare_scores(frames, names, args.metric, plot_title=args.title, output=args.output)


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
