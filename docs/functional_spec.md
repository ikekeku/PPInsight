# Functional Specification

## Background
This project builds a tool to benchmark protein–protein interaction (PPI) prediction models. Experimental PPI determination is slow and costly. Many computational models exist, but outputs are hard to compare. The system retrieves protein sequences/structures, runs multiple PPI models, and visualizes their scores side-by-side.

## User profile
- **Primary users:** Graduate students and researchers in biology, bioengineering, or computer science.  
- **Knowledge:** Basic molecular biology (what a protein is, why interactions matter) and  basic computing (can run Python scripts, use command line, browse the web). No deep ML or structural-biology expertise required.

## Data sources
- **Protein sequences/structures:** FASTA / PDB / mmCIF files from NCBI, UniProt, or AlphaFold.  
- **Model outputs:** text / CSV / JSON files containing interaction scores (e.g., DockQ, RMSD, E-value).

## Use Case 1 - Compare models on one protein pair
- **Objective:** Determine which PPI model gives the best score for a given pair (e.g., proteinA-proteinB).
- **Interaction:** User enters protein names → system fetches FASTA/PDB → user selects 2–3 models → system runs models → writes scores to table → shows plot of scores per model

## Use Case 2 - Curate protein data
- **Objective:** Generate FASTA / PDB files for a batch of proteins.  
- **Interaction:** User uploads CSV of pairs → system fetchs FASTA/PDB file
