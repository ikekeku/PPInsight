#!/usr/bin/env python3
"""
Checkpointed docking pipeline - saves progress and can resume.

This version saves intermediate results so you don't lose progress if
the run crashes or you need to stop it.

Uses PDB files for Pose objects and JSON for results.
"""

import json
from pathlib import Path
import pyrosetta

from ppinsight._rosetta_prepare_structure import prepare_structures, load_structure
from ppinsight._rosetta_dock import (
    run_single_docking, 
    setup_simple_docking, 
    setup_full_docking_protocol
)
from ppinsight._analyze import analyze_scores


class CheckpointedDockingPipeline:
    """
    Docking pipeline with checkpoint support.
    
    Automatically saves progress after each step and can resume from
    where it left off if interrupted.
    
    Example:
        >>> pipeline = CheckpointedDockingPipeline("p1.pdb", "p2.pdb", n_runs=50)
        >>> result = pipeline.run()
        # If this crashes at run 30...
        
        >>> # Just run again - it will resume from run 31!
        >>> pipeline = CheckpointedDockingPipeline("p1.pdb", "p2.pdb", n_runs=50)
        >>> result = pipeline.run()
    """
    
    def __init__(self, protein1_pdb, protein2_pdb, n_runs=10, top_n=20,
                 relax=True, jump_distance=15.0, verbose=True,
                 checkpoint_dir="checkpoints", use_full_protocol=False):
        """
        Initialize checkpointed docking pipeline.
        
        Args:
            protein1_pdb: Path to first protein PDB
            protein2_pdb: Path to second protein PDB
            n_runs: Number of docking runs
            top_n: Number of top scores to average
            relax: If True, relax structures before docking
            jump_distance: Initial separation distance (Å)
            verbose: If True, print progress
            checkpoint_dir: Directory to save checkpoints
            use_full_protocol: If True, use full docking protocol (slower but better)
        """
        self.protein1_pdb = Path(protein1_pdb)
        self.protein2_pdb = Path(protein2_pdb)
        self.n_runs = n_runs
        self.top_n = top_n
        self.relax = relax
        self.jump_distance = jump_distance
        self.verbose = verbose
        self.use_full_protocol = use_full_protocol
        
        self.checkpoint_dir = Path(checkpoint_dir)
        self.checkpoint_dir.mkdir(exist_ok=True)
        
        self.complex_pose = None
        self.docking_results = None
        self.analysis = None
    
    def _get_checkpoint_path(self, name, ext=".pdb"):
        """Get path for a checkpoint file."""
        return self.checkpoint_dir / f"{name}{ext}"
    
    def _save_pose_checkpoint(self, name, pose):
        """Save a Pose object as PDB."""
        checkpoint_path = self._get_checkpoint_path(name, ".pdb")
        try:
            pose.dump_pdb(str(checkpoint_path))
            if self.verbose:
                print(f"  Checkpoint saved: {checkpoint_path}")
            return True
        except Exception as e:
            if self.verbose:
                print(f"  Warning: Could not save checkpoint {checkpoint_path}: {e}")
            return False
    
    def _load_pose_checkpoint(self, name):
        """Load a Pose object from PDB."""
        checkpoint_path = self._get_checkpoint_path(name, ".pdb")
        if checkpoint_path.exists():
            try:
                return load_structure(checkpoint_path)
            except Exception as e:
                if self.verbose:
                    print(f"  Warning: Could not load checkpoint {checkpoint_path}: {e}")
                return None
        return None
    
    def _save_json_checkpoint(self, name, data):
        """Save data as JSON."""
        checkpoint_path = self._get_checkpoint_path(name, ".json")
        with open(checkpoint_path, 'w') as f:
            json.dump(data, f, indent=2)
        if self.verbose:
            print(f"  Checkpoint saved: {checkpoint_path}")
    
    def _load_json_checkpoint(self, name):
        """Load data from JSON."""
        checkpoint_path = self._get_checkpoint_path(name, ".json")
        if checkpoint_path.exists():
            with open(checkpoint_path, 'r') as f:
                return json.load(f)
        return None
    
    def prepare(self, force=False):
        """
        Prepare structures with checkpointing.
        
        Args:
            force: If True, ignore checkpoint and prepare fresh
        
        Returns:
            Combined Pose object
        """
        if self.verbose:
            print("=" * 70)
            print("STEP 1/3: STRUCTURE PREPARATION")
            print("=" * 70)
        
        if not force:
            self.complex_pose = self._load_pose_checkpoint("prepared_complex")
        
        if self.complex_pose is not None:
            if self.verbose:
                print("✓ Loaded prepared structure from checkpoint")
        else:
            if self.verbose:
                print("Preparing structures...")
            
            self.complex_pose = prepare_structures(
                self.protein1_pdb,
                self.protein2_pdb,
                relax=self.relax,
                jump_distance=self.jump_distance,
                verbose=self.verbose
            )
            
            self._save_pose_checkpoint("prepared_complex", self.complex_pose)
            
            if self.verbose:
                print("✓ Structure preparation complete")
        
        return self.complex_pose
    
    def dock(self, force=False):
        """
        Run docking with checkpointing after each run.
        
        Args:
            force: If True, ignore checkpoint and start fresh
        
        Returns:
            List of docking results
        """
        if self.complex_pose is None:
            raise RuntimeError("Must call prepare() before dock()")
        
        if self.verbose:
            print("\n" + "=" * 70)
            print("STEP 2/3: DOCKING")
            print("=" * 70)
        
        # Create docking output directory
        docking_dir = self.checkpoint_dir / "docking_structures"
        docking_dir.mkdir(exist_ok=True)
        
        # Load existing results metadata
        if not force:
            results_meta = self._load_json_checkpoint("docking_results")
        else:
            results_meta = None
        
        if results_meta is None:
            results_meta = []
            completed_runs = 0
        else:
            completed_runs = len(results_meta)
            if self.verbose:
                print(f"✓ Loaded {completed_runs} previous runs from checkpoint")
        
        # Reconstruct docking_results with poses
        self.docking_results = []
        for meta in results_meta:
            pose_file = docking_dir / f"docked_run_{meta['run']:04d}.pdb"
            if pose_file.exists():
                try:
                    pose = load_structure(pose_file)
                    self.docking_results.append({
                        'run': meta['run'],
                        'score': meta['score'],
                        'pose': pose
                    })
                except:
                    # Skip corrupted files
                    if self.verbose:
                        print(f"  Warning: Could not load {pose_file}")
        
        # Run remaining docking simulations
        if completed_runs < self.n_runs:
            if self.verbose:
                print(f"Running docking simulations {completed_runs + 1} to {self.n_runs}...")
            
            # Setup protocol once for efficiency
            if self.use_full_protocol:
                if self.verbose:
                    print("Using full docking protocol (with packing)")
                docking_protocol = setup_full_docking_protocol()
            else:
                if self.verbose:
                    print("Using simple docking protocol")
                docking_protocol = setup_simple_docking()
            
            scorefxn = pyrosetta.create_score_function("ref2015")
            # Disable disulfide scoring
            from pyrosetta import rosetta
            scorefxn.set_weight(rosetta.core.scoring.ScoreType.dslf_fa13, 0.0)
            
            for i in range(completed_runs, self.n_runs):
                if self.verbose:
                    print(f"  Run {i + 1}/{self.n_runs}...")
                
                try:
                    # Run single docking
                    docked_pose, score = run_single_docking(
                        self.complex_pose,
                        docking_protocol,
                        scorefxn,
                        randomize=True
                    )
                    
                    # Only save if we got a valid score (not penalty)
                    if score < 9000.0:  # Valid score
                        result = {
                            'run': i + 1,
                            'score': score,
                            'pose': docked_pose
                        }
                        
                        self.docking_results.append(result)
                        
                        # Try to save docked structure
                        pose_file = docking_dir / f"docked_run_{i + 1:04d}.pdb"
                        save_success = False
                        try:
                            docked_pose.dump_pdb(str(pose_file))
                            save_success = True
                        except Exception as e:
                            if self.verbose:
                                print(f"    Warning: Could not save PDB: {e}")
                            # Continue without saving PDB, just keep the score
                        
                        # Save metadata
                        meta_entry = {
                            'run': i + 1,
                            'score': score,
                        }
                        if save_success:
                            meta_entry['pose_file'] = str(pose_file)
                        
                        results_meta.append(meta_entry)
                        self._save_json_checkpoint("docking_results", results_meta)
                        
                    else:
                        if self.verbose:
                            print(f"    Run {i + 1} failed (penalty score), skipping...")
                    
                    if self.verbose and (i + 1) % 5 == 0:
                        print(f"  Progress: {len(self.docking_results)} successful runs out of {i + 1} attempts")
                    
                except Exception as e:
                    if self.verbose:
                        print(f"  ERROR in run {i + 1}: {e}")
                        print(f"  Skipping this run and continuing...")
                    # Don't raise - just continue to next run
        
        if self.verbose:
            print(f"✓ Docking complete ({len(self.docking_results)} successful runs)")
        
        if len(self.docking_results) == 0:
            raise RuntimeError("No successful docking runs! Check your input structures.")
        
        return self.docking_results
    
    def analyze(self):
        """
        Analyze docking results.
        
        Returns:
            Analysis dictionary
        """
        if self.docking_results is None or len(self.docking_results) == 0:
            raise RuntimeError("No docking results to analyze")
        
        if self.verbose:
            print("\n" + "=" * 70)
            print("STEP 3/3: SCORE ANALYSIS")
            print("=" * 70)
        
        self.analysis = analyze_scores(
            self.docking_results,
            top_n=min(self.top_n, len(self.docking_results)),
            verbose=self.verbose
        )
        
        # Save analysis results
        analysis_summary = {
            'final_score': self.analysis['final_score'],
            'best_score': self.analysis['best_score'],
            'worst_top_score': self.analysis['worst_top_score'],
            'total_runs': self.analysis['total_runs'],
            'top_n': self.analysis['top_n']
        }
        self._save_json_checkpoint("analysis_results", analysis_summary)
        
        return self.analysis
    
    def run(self, force_prepare=False, force_dock=False):
        """
        Run complete pipeline with checkpointing.
        
        Args:
            force_prepare: If True, ignore prepare checkpoint
            force_dock: If True, ignore docking checkpoint
        
        Returns:
            Analysis dictionary
        """
        self.prepare(force=force_prepare)
        self.dock(force=force_dock)
        result = self.analyze()
        
        if self.verbose:
            print("\n✓ Pipeline complete!")
        
        return result
    
    def clear_checkpoints(self):
        """Clear all checkpoints to start fresh."""
        import shutil
        if self.checkpoint_dir.exists():
            shutil.rmtree(self.checkpoint_dir)
            self.checkpoint_dir.mkdir()
            if self.verbose:
                print("✓ Checkpoints cleared")
    
    def get_progress(self):
        """
        Get current progress information.
        
        Returns:
            Dictionary with progress info
        """
        results_meta = self._load_json_checkpoint("docking_results")
        if results_meta:
            completed = len(results_meta)
        else:
            completed = 0
        
        return {
            'completed_runs': completed,
            'total_runs': self.n_runs,
            'progress_percent': (completed / self.n_runs) * 100,
            'remaining_runs': self.n_runs - completed
        }


# Example usage
if __name__ == "__main__":
    # Create pipeline with checkpointing
    pipeline = CheckpointedDockingPipeline(
        protein1_pdb="protein1.pdb",
        protein2_pdb="protein2.pdb",
        n_runs=20,
        verbose=True,
        checkpoint_dir="docking_checkpoints"
    )
    
    # Check progress
    progress = pipeline.get_progress()
    print(f"Progress: {progress['completed_runs']}/{progress['total_runs']} runs")
    
    # Run pipeline - will save checkpoints
    result = pipeline.run()
    
    print(f"\nFinal docking score: {result['final_score']:.2f}")
    
    # If you want to start completely fresh:
    # pipeline.clear_checkpoints()
    # result = pipeline.run()


# from ppinsight.cp_pipeline import CheckpointedDockingPipeline

# path1 = "/home/mayagh/projects/CSE583_projects/PPInsight/examples/ppinsight_data/input_files/2UUY_lig.pdb"
# path2 = "/home/mayagh/projects/CSE583_projects/PPInsight/examples/ppinsight_data/input_files/2UUY_rec.pdb"

# pipeline = CheckpointedDockingPipeline(path1, path2, n_runs=50)
# pipeline.clear_checkpoints()  