# Component Specification

## Software components

### 1. Database Search
- **What it does:** Fetches protein sequences and available 3D structures using UniProt accession IDs and organizes them for downstream modeling.  
- **Inputs:** List of UniProt accession IDs (strings), optional output file paths for FASTA and CSV, and optional output directory for PDB files.  
- **Outputs:** Combined FASTA file, structured CSV metadata file with fields: ID, Name, Description, Sequence Length, Sequence, and downloaded PDB structure files when available.
- **Use of other components:** Uses the UniProt REST API for sequence and structure cross-references and Biopython’s `PDBList` for downloading PDB files.

### 2. Model Selector
- **What it does:** Picks protein-protein interaction (PPI) models to run.  
- **Inputs:** Model of choice, chosen through a dropdown menu.
- **Outputs:** A job request to run proteins of interest through selected models.

### 3. Model Evaluator
- **What it does:** Runs the chosen model on a pair of proteins and calculates scores that quantify the protein interaction.
- **Inputs:** The request output by the model selector, `.fasta` and `.pdb` files output by the database search.
- **Outputs:** For each model run, produces a delimited text file with column headers for each protein and for the score value for each score type.
- **Use of other components:** Requires inputs from the database search and model selector components. May require use of computing cluster (Hyak) to run computationally expensive models.

### 4. Visualization Manager
- **What it does:** Generates interactive plots comparing model scores.  
- **Inputs:** Delimited text file (e.g., `scores.tsv`) and user filters (protein pair, metric, models).  
- **Outputs:** Vega-Lite / Altair visualizations, exportable as PNG, SVG, or PDF. Option to export raw spreadsheet values instead (e.g., download `scores.tsv`).

## Interactions to accomplish Use Case 1 (refer to functional_spec.md)
1. User specifies a protein pair.  
2. A database search is conducted to fetch and store FASTA and PDB files.
4. User selects models to run
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

