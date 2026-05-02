"""
Score analysis module for docking results.

This module handles:
- Sorting and ranking docking results
- Calculating average scores
- Statistical analysis
- Decoy clustering (optional, requires scipy + PyRosetta)
"""

import csv
import logging
import statistics
import warnings
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

PRIMARY_SCORE_COLUMN = 'i_sc'
SECONDARY_SCORE_COLUMN = 'total_score'


def _get_primary_score(result):
    """Return the primary Rosetta ranking score for a decoy result."""
    if PRIMARY_SCORE_COLUMN in result:
        return float(result[PRIMARY_SCORE_COLUMN])
    return float(result['score'])


def _get_total_score(result):
    """Return the Rosetta total energy for a decoy result."""
    if SECONDARY_SCORE_COLUMN in result:
        return float(result[SECONDARY_SCORE_COLUMN])
    return float(result['score'])


def get_top_scores(results, top_n=20):
    """
    Get the top N scoring structures.

    Args:
        results: List of docking results (from run_docking)
        top_n: Number of top scores to return (default: 20)

    Returns:
        List of top N results, sorted by I_sc when available
    """
    # Sort by the primary Rosetta metric (lower is better).
    sorted_results = sorted(results, key=_get_primary_score)

    # Return top N
    return sorted_results[:min(top_n, len(sorted_results))]


def calculate_statistics(scores):
    """
    Calculate statistical measures for a set of scores.

    Args:
        scores: List of score values

    Returns:
        Dictionary with statistical measures:
            - mean: Average score
            - median: Median score
            - stdev: Standard deviation
            - min: Minimum score
            - max: Maximum score
    """
    if not scores:
        return None

    stats = {
        'mean': statistics.mean(scores),
        'median': statistics.median(scores),
        'min': min(scores),
        'max': max(scores)
    }

    if len(scores) > 1:
        stats['stdev'] = statistics.stdev(scores)
    else:
        stats['stdev'] = 0.0

    return stats


def analyze_scores(results, top_n=20, verbose=False):
    """
    Analyze docking results and calculate final score.

    This function:
    1. Sorts results by I_sc when available
    2. Extracts top N scores
    3. Calculates average of top N
    4. Provides statistics

    Args:
        results: List of docking results (from run_docking)
        top_n: Number of top scores to average (default: 20)
        verbose: If True, print detailed analysis (default: False)

    Returns:
        Dictionary with analysis results:
            - final_score: Average of top N scores
            - top_n: Number of scores used
            - best_score: Single best score
            - worst_score: Worst score in top N
            - statistics: Stats for all scores
            - top_results: Top N result dictionaries

    Example:
        >>> results = run_docking(complex_pose, n_runs=100)
        >>> analysis = analyze_scores(results, top_n=20, verbose=True)
        >>> print(f"Final docking score: {analysis['final_score']:.2f}")
    """
    if not results:
        raise ValueError("No results to analyze")

    primary_metric = (
        PRIMARY_SCORE_COLUMN
        if any(PRIMARY_SCORE_COLUMN in result for result in results)
        else 'score'
    )
    primary_label = 'I_sc' if primary_metric == PRIMARY_SCORE_COLUMN else 'Score'

    # Get top N results
    top_results = get_top_scores(results, top_n)
    actual_top_n = len(top_results)

    # Extract scores
    all_scores = [_get_primary_score(result) for result in results]
    top_scores = [_get_primary_score(result) for result in top_results]

    # Calculate final score (average of top N)
    final_score = sum(top_scores) / actual_top_n

    # Calculate statistics
    all_stats = calculate_statistics(all_scores)
    top_stats = calculate_statistics(top_scores)

    # Prepare analysis results
    analysis = {
        'final_score': final_score,
        'primary_metric': primary_metric,
        'top_n': actual_top_n,
        'total_runs': len(results),
        'best_score': top_scores[0],
        'worst_top_score': top_scores[-1],
        'all_statistics': all_stats,
        'top_statistics': top_stats,
        'top_results': top_results
    }

    # Print detailed analysis if requested
    if verbose:
        print("=" * 70)
        print("DOCKING SCORE ANALYSIS")
        print("=" * 70)
        print(f"\nTotal structures analyzed: {len(results)}")
        print(f"Top {actual_top_n} structures used for averaging")
        print()
        print(f"Best {primary_label}:                {analysis['best_score']:8.2f}")
        print(
            f"Worst {primary_label} in top {actual_top_n}:"
            f"     {analysis['worst_top_score']:8.2f}"
        )
        print(
            f"Average {primary_label} of top {actual_top_n}:"
            f"  {analysis['final_score']:8.2f}"
        )
        print()
        print(f"Statistics (all {primary_label} values):")
        print(f"  Mean:       {all_stats['mean']:8.2f}")
        print(f"  Median:     {all_stats['median']:8.2f}")
        print(f"  Std Dev:    {all_stats['stdev']:8.2f}")
        print(f"  Range:      {all_stats['min']:8.2f} to {all_stats['max']:8.2f}")
        print()
        print("=" * 70)
        print(f"FINAL DOCKING SCORE: {final_score:.2f}")
        print("=" * 70)

    return analysis


def print_top_scores(results, top_n=10):
    """
    Print a table of top N scores.

    Args:
        results: List of docking results
        top_n: Number of top scores to print (default: 10)
    """
    top_results = get_top_scores(results, top_n)
    primary_metric = (
        PRIMARY_SCORE_COLUMN
        if any(PRIMARY_SCORE_COLUMN in result for result in top_results)
        else 'score'
    )
    primary_label = 'I_sc' if primary_metric == PRIMARY_SCORE_COLUMN else 'Score'
    show_total_score = any(
        SECONDARY_SCORE_COLUMN in result for result in top_results
    )

    print(f"\nTop {len(top_results)} Docking Scores:")
    print("-" * 56 if show_total_score else "-" * 40)
    if show_total_score:
        print(f"{'Rank':<6} {'Run':<6} {primary_label:>12} {'Total Score':>14}")
        print("-" * 56)
    else:
        print(f"{'Rank':<6} {'Run':<6} {primary_label:>12}")
        print("-" * 40)

    for i, result in enumerate(top_results, 1):
        run_num = result['run']
        score = _get_primary_score(result)
        if show_total_score:
            total_score = _get_total_score(result)
            print(f"{i:<6} {run_num:<6} {score:>12.2f} {total_score:>14.2f}")
        else:
            print(f"{i:<6} {run_num:<6} {score:>12.2f}")

    print("-" * 56 if show_total_score else "-" * 40)


def export_scores_to_csv(results, output_path):
    """
    Export all Rosetta per-decoy scores to a CSV file.

    Args:
        results: List of docking results
        output_path: Path for output CSV file
    """

    with open(output_path, 'w', encoding="utf8", newline='') as f:
        writer = csv.writer(f)

        # Header
        writer.writerow(['run', 'total_score', 'i_sc'])

        # Data
        for result in results:
            writer.writerow([
                result['run'],
                _get_total_score(result),
                _get_primary_score(result),
            ])


# ---------------------------------------------------------------------------
# Decoy clustering (optional)
# ---------------------------------------------------------------------------

def cluster_and_rank(
    scores_csv: str | Path,
    pdb_dir: str | Path,
    *,
    score_col: str = "i_sc",
    top_n: int = 200,
    rmsd_cutoff: float = 4.0,
    output_csv: str | Path | None = None,
) -> pd.DataFrame | None:
    """Cluster Rosetta decoys by Cα-RMSD and (optionally) write results.

    This is a thin convenience wrapper around
    :func:`ppinsight.rosetta.cluster.cluster_decoys` that:

    * reads a CSV/TSV score file
    * runs clustering
    * saves the annotated DataFrame to *output_csv* (default:
      ``clustered_scores.csv`` next to *scores_csv*)
    * returns the clustered DataFrame

    If ``scipy`` or ``pyrosetta`` are not installed the function issues a
    warning and returns *None* instead of raising.

    Parameters
    ----------
    scores_csv : path-like
        Score file (CSV produced by the pipeline with explicit Rosetta
        metrics, or a ``.sc`` file previously converted to CSV).
    pdb_dir : path-like
        Directory containing the decoy PDB files.
    score_col, top_n, rmsd_cutoff
        Forwarded to :func:`cluster_decoys`.
    output_csv : path-like, optional
        Where to write the clustered output.  Defaults to
        ``<scores_csv_dir>/clustered_scores.csv``.

    Returns
    -------
    DataFrame or None
    """
    try:
        from ppinsight.rosetta.cluster import cluster_decoys
    except ImportError as exc:
        warnings.warn(
            f"Cannot cluster decoys — missing dependency: {exc}. "
            "Install scipy and pyrosetta for clustering support.",
            stacklevel=2,
        )
        return None

    scores_csv = Path(scores_csv)
    pdb_dir = Path(pdb_dir)

    # Read scores
    sep = "\t" if scores_csv.suffix == ".tsv" else ","
    df = pd.read_csv(scores_csv, sep=sep)
    df.columns = [c.strip().lower() for c in df.columns]

    try:
        result = cluster_decoys(
            df, pdb_dir,
            score_col=score_col,
            top_n=top_n,
            rmsd_cutoff=rmsd_cutoff,
        )
    except Exception as exc:
        warnings.warn(
            f"Clustering failed: {exc}",
            stacklevel=2,
        )
        return None

    # Write output — defaults to "clustered_scores.csv" next to the
    # original scores file, which is what collect_scores looks for first.
    if output_csv is None:
        output_csv = scores_csv.parent / "clustered_scores.csv"
    result.to_csv(output_csv, index=False)
    logger.info("Clustered scores written to %s", output_csv)

    return result
