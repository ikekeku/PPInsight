# Component Specification

## Software components

### 1. Database Search
- **What it does:** Retrieves and organizes protein data.  
- **Inputs:** Protein name (or IDs), species, optional parameters (e.g., number of search results found).  
- **Outputs:** `*.fasta` and/or `*.pdb` files, metadata table,  organized folders.

### 2. Model Selector
- **What it does:** Pick protein-protein interaction (PPI) models to run.  
- **Inputs:** Model of choice, chosen through a dropdown menu.  
- **Outputs:** A job request to run proteins of interest through selected models.

### 3. Visualization Manager
- **What it does:** Generates interactive plots comparing model scores.  
- **Inputs:** Delimited text file (e.g., `scores.tsv`) and user filters (protein pair, metric, models).  
- **Outputs:** Vega-Lite / Altair visualizations, exportable as PNG, SVG, or PDF. Option to export raw spreadsheet values instead (e.g., download `scores.tsv`).

## Interactions to accomplish Use Case 1 (refer to functional_spec.md)
1. User specifies a protein pair.  
2. A database search is conducted to fetch and store FASTA/PDB files.  
3. User selects models to run and results are appended to `scores.tsv`.  
4. The visualization manager reads `scores.tsv` and renders a chart (e.g., x = model, y = score).  
5. User inspects plot and exports results.

## Preliminary plan (priority)
- Task 1: fetch proteins by search term and save files (i.e., implement database search)
- Task 2: run PPI predictors on file pairs and store scores (i.e., model selection and data analysis)
- Task 3: plot interaction scores (i.e., visualization manager)
- Elective Task #1: Add batch mode for multiple pairs.  

