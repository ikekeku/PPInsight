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

### 2. Model Evaluator (e.g., pdb_to_haddock, pdb_to_lightdock, rosetta_docking)
- **What it does:** Stages inputs, writes tool-specific configs, and runs the chosen model on a pair of proteins and calculates scores that quantify the protein interaction.
- **Inputs:** The request output by the model selector (`.fasta` and `.pdb` files output by the database search) and respective paths (absolute or repo basenames), runtime options (cores, mode, container, runname, etc.)
- **Outputs:**
    - Run directory with `data/` containing staged inputs
    - Tool config file (e.g., `*.cfg` for HADDOCK, `setup.json` for LightDock)
    - Execution summary with opts passed (executed command or container invocation)
    - Associated output files models ran mapping proteinA/proteinB (e.g., scoring and docked models [`.pdb`s])
- **Use of other components:** Requires inputs from the database search and model selector components. May require use of computing cluster (Hyak) to run computationally expensive models.

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
1. User specifies a protein pair.  
2. A database search is conducted to fetch and store FASTA and PDB files.
4. User determines models to run
5. For each selected model:
   - Input data is converted into the model’s required format.
   - Necessary configuration files are generated.
   - Jobs are dispatched to Hyak if parallel computation is required.
8. Results from the model runs are appended to `scores.tsv`.  
9. The visualization manager reads `scores.tsv` and renders a chart (e.g., x = model, y = score).  
10. User inspects plot and exports results.

## Preliminary plan (priority)
- Task 1: fetch proteins by search term and save files (i.e., implement database search)
- Task 2: run PPI predictors on file pairs and store scores (i.e., model selection and data analysis)
- Task 3: plot interaction scores (i.e., visualization manager)
- Elective Task #1: Add batch mode for multiple pairs.  

