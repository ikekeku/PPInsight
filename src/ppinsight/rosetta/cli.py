"""
Docstring for ppinsight.rosetta.cli
This file accepts command-line arguments to run the DockingPipeline from ppinsight.rosetta.pipeline
and outputs the final docking score and saves the top structures.

It takes as input two protein PDB files, number of docking runs, output directory, 
and number of top structures to save.
"""

import argparse
from ppinsight.rosetta.pipeline import DockingPipeline

def main():
    parser = argparse.ArgumentParser(
        description="Run Rosetta protein–protein docking with PPInsight"
    )
    parser.add_argument("protein1", help="Path to first protein PDB")
    parser.add_argument("protein2", help="Path to second protein PDB")
    parser.add_argument("--n-runs", type=int, default=50)
    parser.add_argument("--outdir", default="output")
    parser.add_argument("--top-n", type=int, default=5)

    args = parser.parse_args()

    pipeline = DockingPipeline(
        args.protein1,
        args.protein2,
        n_runs=args.n_runs,
    )

    result = pipeline.run()
    pipeline.save_top_structures(args.outdir, top_n=args.top_n)

    print(f"Final docking score: {result['final_score']:.2f}")

if __name__ == "__main__":
    main()
