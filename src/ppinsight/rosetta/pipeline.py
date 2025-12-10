"""
Main pipeline class for protein-protein docking.

This module provides a high-level interface for the complete docking workflow.
"""

# pylint: disable=too-many-instance-attributes

from pathlib import Path
from .prepare_structure import prepare_structures
from .dock import run_docking, save_docked_structure
from .analyze import analyze_scores, export_scores_to_csv


class DockingPipeline:
    """
    High-level interface for protein-protein docking.

    This class encapsulates the complete docking workflow:
    1. Structure preparation and relaxation
    2. Multiple docking runs
    3. Score analysis and ranking

    Attributes:
        protein1_pdb: Path to first protein PDB file
        protein2_pdb: Path to second protein PDB file
        n_runs: Number of docking runs to perform
        top_n: Number of top scores to average for final result

    Example:
        >>> path1 = "../../PPInsight/examples/ppinsight_data/input_files/2UUY_lig.pdb"
        >>> path2 = "../../PPInsight/examples/ppinsight_data/input_files/2UUY_rec.pdb"
        >>> pipeline = DockingPipeline(path1, path2, n_runs=50)
        >>> result = pipeline.run()
        >>> print(f"Final docking score: {result['final_score']:.2f}")

        >>> # Save top structures
        >>> pipeline.save_top_structures("output_dir/", top_n=5)
    """

    # pylint: disable=too-many-arguments,too-many-positional-arguments
    def __init__(self, protein1_pdb, protein2_pdb, n_runs=10, top_n=20,
                 relax=True, jump_distance=15.0, verbose=True):
        """
        Initialize the docking pipeline.

        Args:
            protein1_pdb: Path to first protein PDB file
            protein2_pdb: Path to second protein PDB file
            n_runs: Number of docking runs (default: 10)
            top_n: Number of top scores to average (default: 20)
            relax: If True, relax structures before docking (default: True)
            jump_distance: Initial separation distance in Å (default: 15.0)
            verbose: If True, print progress messages (default: True)
        """
        self.protein1_pdb = Path(protein1_pdb)
        self.protein2_pdb = Path(protein2_pdb)
        if n_runs < 1:
            raise ValueError("n_runs must be at least 1")
        self.n_runs = n_runs

        if top_n < 1:
            raise ValueError("top_n must be at least 1")

        self.top_n = top_n
        self.relax = relax
        self.jump_distance = jump_distance
        self.verbose = verbose

        # Results storage
        self.complex_pose = None
        self.docking_results = None
        self.analysis = None

    def prepare(self):
        """
        Prepare protein structures for docking.

        This step:
        - Loads both PDB files
        - Relaxes structures (if enabled)
        - Combines them into a complex

        Returns:
            Combined Pose object
        """
        if self.verbose:
            print("=" * 70)
            print("STEP 1/3: STRUCTURE PREPARATION")
            print("=" * 70)

        self.complex_pose = prepare_structures(
            self.protein1_pdb,
            self.protein2_pdb,
            relax=self.relax,
            jump_distance=self.jump_distance,
            verbose=self.verbose
        )

        if self.verbose:
            print("✓ Structure preparation complete")

        return self.complex_pose

    def dock(self):
        """
        Run docking simulations.

        This step runs multiple independent docking simulations
        starting from the prepared complex.

        Returns:
            List of docking results
        """
        if self.complex_pose is None:
            raise RuntimeError("Must call prepare() before dock()")

        if self.verbose:
            print("\n" + "=" * 70)
            print("STEP 2/3: DOCKING")
            print("=" * 70)

        self.docking_results = run_docking(
            self.complex_pose,
            n_runs=self.n_runs,
            save_all=True,
            verbose=self.verbose
        )

        if self.verbose:
            print(f"✓ Docking complete ({self.n_runs} runs)")

        return self.docking_results

    def analyze(self):
        """
        Analyze docking results.

        This step:
        - Ranks structures by score
        - Calculates average of top N scores
        - Provides statistical analysis

        Returns:
            Analysis dictionary with final score and statistics
        """
        if self.docking_results is None:
            raise RuntimeError("Must call dock() before analyze()")

        if self.verbose:
            print("\n" + "=" * 70)
            print("STEP 3/3: SCORE ANALYSIS")
            print("=" * 70)

        self.analysis = analyze_scores(
            self.docking_results,
            top_n=self.top_n,
            verbose=self.verbose
        )

        return self.analysis

    def run(self):
        """
        Run the complete docking pipeline.

        Returns:
            Analysis dictionary with final score and statistics
        """
        self.prepare()
        self.dock()
        result = self.analyze()

        if self.verbose:
            print("\n✓ Pipeline complete!")

        return result

    def save_top_structures(self, output_dir, top_n=5):
        """
        Save the top N docked structures to PDB files.

        Args:
            output_dir: Directory to save structures
            top_n: Number of top structures to save (default: 5)
        """
        if self.analysis is None:
            raise RuntimeError("Must run pipeline before saving structures")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        top_results = self.analysis['top_results'][:top_n]

        for i, result in enumerate(top_results, 1):
            output_path = output_dir / f"docked_top_{i}_score_{result['score']:.2f}.pdb"
            save_docked_structure(result['pose'], output_path)

            if self.verbose:
                print(f"Saved: {output_path}")

    def save_scores(self, output_path):
        """
        Save all scores to a CSV file.
        """
        if self.docking_results is None:
            raise RuntimeError("Must run docking before saving scores")

        export_scores_to_csv(self.docking_results, output_path)

        if self.verbose:
            print(f"Scores saved to: {output_path}")

    def get_best_structure(self):
        """
        Get the best scoring docked structure.

        Returns:
            Dictionary with 'pose' and 'score' for best structure
        """
        if self.analysis is None:
            raise RuntimeError("Must run pipeline before getting best structure")

        return self.analysis['top_results'][0]

    def print_summary(self):
        """
        Print a summary of docking results.
        """
        if self.analysis is None:
            raise RuntimeError("Must run pipeline before printing summary")

        print("\n" + "=" * 70)
        print("DOCKING PIPELINE SUMMARY")
        print("=" * 70)
        print(f"Protein 1: {self.protein1_pdb.name}")
        print(f"Protein 2: {self.protein2_pdb.name}")
        print(f"Number of docking runs: {self.n_runs}")
        print(f"Top N for averaging: {self.top_n}")
        print()
        print(f"Best score: {self.analysis['best_score']:.2f}")
        print(f"Final score (avg of top {self.top_n}): {self.analysis['final_score']:.2f}")
        print("=" * 70)
