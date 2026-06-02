# Use Cases 
This document describes concrete, actionable use cases supported by the
repository. Each use case lists actors, preconditions, steps, outputs, and
success criteria.

### USE CASE A / *Task 1: Fetch proteins by accession ID and save files*
- User: Researcher
- System: PPInsight fetch command (`ppinsight fetch`)
- Preconditions: network access, list of UniProt accession IDs
- Steps:
  1. User provides one or more UniProt accession IDs.
  2. `ppinsight fetch` queries UniProt and resolves available PDB cross-references.
  3. The command writes a FASTA (optional), metadata CSV, and accession-named
     PDB files into the chosen `--pdb-dir` when structures are available.
- Outputs: FASTA, CSV metadata, PDB files
- Success criteria: CSV has one row per accession; PDB files exist when
  available; function returns structured metadata and PDB mapping.


### USE CASE B / *Task 2: Run PPI predictors on file pairs and store scores*
- Actors: Researcher
- System: PPInsight docking commands (`ppinsight haddock`, `ppinsight lightdock`, `ppinsight rosetta`)
- Preconditions: receptor and ligand PDBs are available (from Use Case 1 or
  provided locally)
- Steps (single model example):
  1. User determines model to use (HADDOCK, LightDock, or Rosetta).
  2. Run the corresponding helper:
     - HADDOCK: `ppinsight haddock receptor ligand --run`
     - LightDock: `ppinsight lightdock receptor ligand`
     - Rosetta: `ppinsight rosetta receptor ligand`
  3. The runner stages inputs under `data/output/<engine>_runs/<ProteinA>_vs_<ProteinB>/`.
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
- System: Visualization command (`ppinsight compare`)
- Preconditions: `scores.tsv` is present
- Steps:
  1. Use `ppinsight compare` to render summaries and plots.
  2. Filter by protein pair and metric. Choose models to display.
  3. Export plots as PNG/SVG.
- Outputs: terminal summaries and PNG/SVG/PDF visualizations
- Success criteria: plots render without errors and reflect `scores.tsv` data.


### USE CASE D: *Batch-dock all pairs from an interaction table*
- User: Researcher
- System: PPInsight batch pipeline (`ppinsight parse`, `ppinsight batch`,
  `ppinsight collect`, `ppinsight compare`)
- Preconditions: an annotation table listing protein pairs with
  interaction / non-interaction labels; network access (or pre-fetched PDB
  files); at least one docking engine installed (LightDock, HADDOCK, or
  Rosetta)
- Steps:
  1. Parse the annotation table into a flat pairs file:
     `ppinsight parse table.tsv -o pairs.csv --stats`
  2. Run docking for every pair across one or more engines:
      `ppinsight batch pairs.csv --engines lightdock haddock --pdb-dir data/input/`
     (use `--dry-run` first to verify resolution and `--limit N` for testing)
  3. Collect scores from the engine run directories into a unified scores file.
     Use `--pair ProteinA:ProteinB` for a single known pair.
      Use `--label-file pairs.csv` (legacy alias: `--pairs`) when you want
      labels annotated from the parsed pairs file.
     Examples:
     `ppinsight collect output_dir... --pair ProteinA:ProteinB -o scores.tsv`
      `ppinsight collect output_dir... --label-file pairs.csv -o scores.tsv`
  4. Compare engines with classification metrics and ROC curves:
     `ppinsight compare scores.tsv -m luciferin_score --classify`
     `ppinsight compare scores.tsv -m score --plot-type roc`
  5. Optionally normalise and aggregate for cross-engine comparison:
     `ppinsight compare scores.tsv -m score --normalize minmax --agg best`
- Outputs: `pairs.csv`, `batch_results.csv`, `scores.tsv`, confusion-matrix
    tables, ROC plots, and compare-CLI visualisations
- Success criteria: every pair × engine combination produces a row in
  `batch_results.csv` with status `success` or a logged error;
  `scores.tsv` contains labelled scores for all successful runs; ROC AUC
  is computable for engines with both interaction and non-interaction pairs.
