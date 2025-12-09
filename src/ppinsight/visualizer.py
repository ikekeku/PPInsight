"""
Generates a pyplot bar chart comparing chosen models for a given score type
"""

from matplotlib import pyplot as plt
import os
import pandas as pd
from pathlib import Path

def to_plot(models: list[os.PathLike]) -> list[pd.DataFrame]:
    """Create a list of dataframes for the models you want to plot"""
    if type(models) is not list:
        raise TypeError('Format your input as a list of file paths')
    else:
        pass
    frames = []
    for m in models:
        df = pd.read_csv(m)
        frames.append(df)
    return frames

def compare_scores(frames: list[pd.DataFrame], score_type: str, plot_title: str):
    """Generate bar chart showing the value of score_type for each model in frames"""
    scores = []
    for file in frames:
        if score_type not in file.columns:
            raise LookupError(f'The score type {score_type} does not exist for the {file} model')
        else:
            pass
        ind = file.columns.get_iloc(score_type)
        cell = file.iloc[0,ind]
        scores.append(cell)
    plt.bar(frames, scores)
    plt.title(plot_title)
    plt.xlabel('Interaction model')
    plt.ylabel(f'{score_type} value')
    chart = plt.show()
    return chart