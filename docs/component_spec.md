# Component Specification
This file gives concise component-level contracts for the main pieces of
PPInsight. Each component lists responsibilities, inputs, outputs and
supplementary components relevant to users/developers.

## Software components

### 1. Database Search
- **What it does:** Fetches protein sequences and available 3D structures using UniProt accession IDs and organizes them for downstream modeling.  
- **Inputs:** List of UniProt accession IDs (strings), optional output file paths for FASTA and CSV, and optional output directory for PDB files.  
- **Outputs:** Combined FASTA file, structured CSV metadata file with fields: ID, Name, Description, Sequence Length, Sequence, and downloaded PDB structure files when available.
- **Use of other components:** Uses the UniProt REST API for sequence and structure cross-references and Biopython’s `PDBList` for downloading PDB files.

### 2. Model Evaluator (for example, `ppinsight haddock`, `ppinsight lightdock`, `ppinsight rosetta`)
- **What it does:** Stages inputs, writes tool-specific configs, and runs the chosen model on a pair of proteins and calculates scores that quantify the protein interaction.
- **Inputs:** Pair identifiers plus corresponding PDB structures (either fetched accession-named files or user-provided local files), path hints (absolute or repo basenames), and runtime options (cores, mode, container, runname, etc.)
- **Outputs:**
    - Run directory with `data/` containing staged inputs
    - Tool config file (e.g., `*.cfg` for HADDOCK, `setup.json` for LightDock)
    - Execution summary with opts passed (executed command or container invocation)
    - Associated output files models ran mapping proteinA/proteinB (e.g., scoring and docked models [`.pdb`s])
  - **Use of other components:** Requires inputs from the database search and pairs-file conventions. May require use of computing cluster (Hyak) to run computationally expensive models.

#### 2.5. Scorer / Aggregator
- **What it does:** Collects model outputs and writes normalized `scores.tsv` used by the
  visualizer.
- **Inputs:** One or more model output folders.
- **Outputs:** `scores.tsv` with columns: proteinA, proteinB, model, score_type,
  score_value, output_path, timestamp
-**Notes**:
  - Implementations should normalize metric names (dockq, rmsd, evalue)
  - Add provenance fields (command-line, container image) for reproducibility

### 3. Visualization Manager
- **What it does:** Generates plots comparing model scores.  
- **Inputs:** Delimited text file (e.g., `scores.tsv`) and user filters (protein pair, metric, models).
- **Outputs:** List of pandas DataFrames containing score values, one DataFrame for each model. Pyplot visualizations, exportable as PNG, SVG, or PDF.

## Interactions to accomplish Use Case 1 (refer to functional_spec.md)
1. User fetches proteins by accession ID, producing FASTA/CSV and accession-named PDB files.
2. User specifies a protein pair directly or through a pairs file.
3. User runs one or more docking engines (single-engine command or `ppinsight batch`).
4. Docking commands write run directories and engine artifacts; `batch` also writes `batch_results.csv`.
5. User runs `ppinsight collect` on run directories to build unified `scores.tsv`.
6. The visualization manager reads `scores.tsv` and renders summaries/plots for comparison and export.
7. User inspects outputs and exports results.

## Preliminary plan (priority)
- Task 1: fetch proteins by accession ID and save files (i.e., implement database search)
- Task 2: run PPI predictors on file pairs and store scores (i.e., model selection and data analysis)
- Task 3: plot interaction scores (i.e., visualization manager)
- Elective Task #1: Add batch mode for multiple pairs.  

