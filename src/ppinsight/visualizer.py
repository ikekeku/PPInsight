"""
Generates a pyplot bar chart comparing chosen models for a given score type
"""

import os
from typing import List

from matplotlib import pyplot as plt
import pandas as pd


# NOTE: Use `typing.List[...]` instead of PEP 585-style `list[...]` to keep
# compatibility with Python 3.8 (the project supports >=3.8). Using `list[...]`
# here caused an import-time TypeError on Python 3.8 where `list` is not
# subscriptable. This is a minimal, behavior-preserving change.
def to_plot(models: List[os.PathLike]) -> List[pd.DataFrame]:
    """Create a list of dataframes for the models you want to plot."""
    # Use the more idiomatic form for isinstance checks
    if not isinstance(models, list):
        raise TypeError("Format your input as a list of file paths")

    frames: List[pd.DataFrame] = []
    for m in models:
        df = pd.read_csv(m, header=0)
        frames.append(df)
    return frames


def compare_scores(frames: List[pd.DataFrame], models: List[str], score_type: str, plot_title: str):
    """Generate bar chart showing the value of score_type for each model in frames."""
    scores: List[float] = []
    for df in frames:
        if score_type not in df.columns:
            raise LookupError(f"The score type {score_type} does not exist for the {df} model")
        # Keep original behavior (index 0) for a single-value per-frame summary
        scores.append(df['score_a'][0])
    plt.bar(models, scores)
    plt.title(plot_title)
    plt.xlabel('Interaction model')
    plt.ylabel(f'{score_type} value')
    chart = plt.show()
    return chart
