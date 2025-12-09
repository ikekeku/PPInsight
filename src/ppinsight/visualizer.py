"""
Generates a pyplot bar chart comparing chosen models for a given score type
"""

import os

from matplotlib import pyplot as plt
import pandas as pd

def to_plot(models: list[os.PathLike]) -> list[pd.DataFrame]:
    """Create a list of dataframes for the models you want to plot"""
    if isinstance(models,list) is False:
        raise TypeError('Format your input as a list of file paths')
    frames = []
    for m in models:
        df = pd.read_csv(m,header=0)
        frames.append(df)
    return frames

def compare_scores(frames: list[pd.DataFrame], models: list[str], score_type: str, plot_title: str):
    """Generate bar chart showing the value of score_type for each model in frames"""
    scores = []
    for df in frames:
        if score_type not in df.columns:
            raise LookupError(f'The score type {score_type} does not exist for the {df} model')
        scores.append(df['score_a'][0])
    plt.bar(models, scores)
    plt.title(plot_title)
    plt.xlabel('Interaction model')
    plt.ylabel(f'{score_type} value')
    chart = plt.show()
    return chart
