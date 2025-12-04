"""
Score analysis module for docking results.

This module handles:
- Sorting and ranking docking results
- Calculating average scores
- Statistical analysis
"""

from typing import List, Dict, Any
import statistics


def get_top_scores(results, top_n=20):
    """
    Get the top N scoring structures.
    
    Args:
        results: List of docking results (from run_docking)
        top_n: Number of top scores to return (default: 20)
    
    Returns:
        List of top N results, sorted by score (best first)
    """
    # Sort by score (lower is better)
    sorted_results = sorted(results, key=lambda x: x['score'])
    
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
    1. Sorts results by score
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
    
    # Get top N results
    top_results = get_top_scores(results, top_n)
    actual_top_n = len(top_results)
    
    # Extract scores
    all_scores = [r['score'] for r in results]
    top_scores = [r['score'] for r in top_results]
    
    # Calculate final score (average of top N)
    final_score = sum(top_scores) / actual_top_n
    
    # Calculate statistics
    all_stats = calculate_statistics(all_scores)
    top_stats = calculate_statistics(top_scores)
    
    # Prepare analysis results
    analysis = {
        'final_score': final_score,
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
        print(f"Best score:               {analysis['best_score']:8.2f}")
        print(f"Worst score in top {actual_top_n}:     {analysis['worst_top_score']:8.2f}")
        print(f"Average of top {actual_top_n}:        {analysis['final_score']:8.2f}")
        print()
        print("Statistics (all scores):")
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
    
    print(f"\nTop {len(top_results)} Docking Scores:")
    print("-" * 40)
    print(f"{'Rank':<6} {'Run':<6} {'Score':>12}")
    print("-" * 40)
    
    for i, result in enumerate(top_results, 1):
        run_num = result['run']
        score = result['score']
        print(f"{i:<6} {run_num:<6} {score:>12.2f}")
    
    print("-" * 40)


def export_scores_to_csv(results, output_path):
    """
    Export all scores to a CSV file.
    
    Args:
        results: List of docking results
        output_path: Path for output CSV file
    """
    import csv
    
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Header
        writer.writerow(['run', 'score'])
        
        # Data
        for result in results:
            writer.writerow([result['run'], result['score']])
