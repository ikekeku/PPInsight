# Use Cases 
This document describes concrete, actionable use cases supported by the
repository. Each use case lists actors, preconditions, steps, outputs, and
success criteria.

### USE CASE A / *Task 1: Fetch proteins by search term and save files*
- User: Researcher
- System: PPInsight fetcher (`protein_fetch.py`)
- Preconditions: network access, list of UniProt accessions or search terms
- Steps:
  1. User provides accession list or search term.
  2. `ppinsight.protein_fetch.get_uniprot_data(...)` is invoked.
  3. The fetcher queries UniProt, writes a FASTA (optional), metadata CSV,
     and downloads PDB files into `pdb_files/` when available.
- Outputs: FASTA, CSV metadata, PDB files
- Success criteria: CSV has one row per accession; PDB files exist when
  available; function returns structured metadata and PDB mapping.


### USE CASE B / *Task 2: Run PPI predictors on file pairs and store scores*
- Actors: Researcher
- System: PPInsight pipeline scripts (e.g., `pdb_to_haddock` and `pdb_to_lightdock`)
- Preconditions: receptor and ligand PDBs are available (from Use Case 1 or
  provided locally)
- Steps (single model example):
  1. User determines model to use (HADDOCK, LightDock, or Rosetta).
  2. Run the corresponding helper:
     - HADDOCK: `pdb_to_haddock.py receptor ligand [--run]`
     - LightDock: `pdb_to_lightdock.py receptor ligand [--generate]`
     - Rosetta: use `ppinsight.rosetta_docking.DockingPipeline`
  3. The runner stages inputs under `.../output_files/<ProteinAvsProteinB>/<runname>/data/`.
  4. Configs are created (`.cfg` or `setup.json`) and are optionally executed
     locally or inside a container.
- Outputs: staged run folder, config file, optional model outputs (e.g., swarm files and
  docked PDBs), and a command summary
- Success criteria: run directory created with `data/` and a config file.
#### *Task 2.5: Scorer / Aggregator*
- System: Aggregator script
- Preconditions: model runs completed, output files available
- Steps:
  1. Parse model outputs to extract numeric metrics (DockQ, RMSD, Fnat, E-values).
  2. Normalize metric names and write rows to `scores.tsv` with provenance.
  3. Add entries for unsuccessful runs with an error tag.
- Outputs: `scores.tsv` or `scores.csv`
- Success criteria: each row contains required columns and a parsable numeric
  score_value when applicable.

### USE CASE C / *Task 3: Plot interaction scores*
- User: Researcher
- System: Visualizer (`visualizer.py`)
- Preconditions: `scores.tsv` is present
- Steps:
  1. Use `ppinsight.visualizer` or to render plots.
  2. Filter by protein pair and metric. Choose models to display.
  3. Export plots as PNG/SVG.
- Outputs: images and interactive HTML visualizations
- Success criteria: plots render without errors and reflect `scores.tsv` data.


### USE CASE D: *Batch-dock all pairs from an interaction table*
- User: Researcher
- System: PPInsight batch pipeline (`parse_pairs.py`, `batch_dock.py`,
  `collect_scores.py`, `visualizer.py`)
- Preconditions: an annotation table listing protein pairs with
  interaction / non-interaction labels; network access (or pre-fetched PDB
  files); at least one docking engine installed (LightDock, HADDOCK, or
  Rosetta)
- Steps:
  1. Parse the annotation table into a flat pairs file:
     `ppinsight parse table.tsv -o pairs.csv --stats`
  2. Run docking for every pair across one or more engines:
     `ppinsight batch pairs.csv --engines lightdock haddock --pdb-dir pdb_files/`
     (use `--dry-run` first to verify resolution and `--limit N` for testing)
  3. Collect scores from all output directories into a unified scores file,
     annotated with interaction labels:
     `ppinsight collect output_dirs... -o scores.tsv --pairs pairs.csv --summary`
  4. Compare engines with classification metrics and ROC curves:
     `ppinsight compare scores.tsv -m luciferin_score --classify`
     `ppinsight compare scores.tsv -m score --plot-type roc`
  5. Optionally normalise and aggregate for cross-engine comparison:
     `ppinsight compare scores.tsv -m score --normalize minmax --agg best`
- Outputs: `pairs.csv`, `batch_results.csv`, `scores.tsv`, confusion-matrix
  tables, ROC plots, bar/violin/heatmap visualisations
- Success criteria: every pair × engine combination produces a row in
  `batch_results.csv` with status `success` or a logged error;
  `scores.tsv` contains labelled scores for all successful runs; ROC AUC
  is computable for engines with both interaction and non-interaction pairs.


<!-- ##### EVERYONE RESTRUCTURE THE USE CASES CLEARLY UNDER THE ACCORDING TASK!!! 
SETUP 
- Accessing a remote server (HYAK?)  

- Interpreting a protein name and pulling FASTA files from the internet:
The input of this function will be the protein name and the species of interest. The output of this function will be the according FASTA file for the protein of interest. The function itself will be able to input the two designations into the search bar of an according database in order to pull the file output.
- Retrieve PDB structure from PDB ID or protein name
- If multiple FASTA files are returned (e.g., isoforms), display all candidates and let the user pick the one(s) of interest
- Option to upload a local FASTA file
-  Save the FASTA files of one or multiple proteins into one folder
- Taking the designated FASTA files and putting them into models to run
- Having chosen the proteins and saved the FASTA files, compare two models of the interaction between the selected proteins:
The inputs of this component are the models that the user wants to compare and the FASTA files of the selected proteins. This requires a menu of models from which the user is able to select at least (or exactly?) two. It also requires either a place where the user uploads the FASTA files, or where the program pulls the FASTA files from an earlier component. The output will be a CSV (or similar) file that contains the scores which will then serve as inputs for a visualization/plotting step.
- Output of the model result and keep as a value in your local server 
- Plot the scores that compare models
- Filter plots by score type or protein pair
- Export plot option: allow PNG or SVG export of graphs
- The system can compile all results into a text summary report, listing models and metrics, as well as key findings (input by user)
   -->