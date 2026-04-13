# `data/` — PPInsight working directory

This is the **default I/O root** for real pipeline runs.  When you don't
pass explicit `--input-dir`, `--pdb-dir`, or `--output-root` flags, the
pipeline reads from `data/input/` and writes to `data/output/`.

> **Not for examples or tutorials** — those live in `examples/`.
> This directory is where you work when you are actually using the pipeline.

## Directory layout

```
data/
├── input/
│   ├── *.pdb              # Receptor and ligand PDB structures
│   └── pairs/             # Batch-mode CSV / TSV files go here
│       ├── pairs.csv      #   (output of `ppinsight parse`)
│       └── <table>.tsv    #   (your protein interaction table)
│
└── output/
    ├── haddock_runs/      # One subdirectory per HADDOCK docking run
    │   └── <ProteinA>_vs_<ProteinB>/
    ├── lightdock_runs/    # One subdirectory per LightDock run
    │   └── <ProteinA>_vs_<ProteinB>/
    ├── rosetta_runs/      # One subdirectory per Rosetta run
    │   └── <ProteinA>_vs_<ProteinB>/
    ├── scores/            # Unified scores & batch results
    │   ├── scores.tsv
    │   ├── scores.tsv.provenance.json
    │   └── batch_results.csv
    └── plots/             # Saved figures
        └── *.png / *.svg / *.pdf
```

## Where things go

### Inputs — where to put it

| What                        | Default location              | CLI flag to override         |
|-----------------------------|-------------------------------|------------------------------|
| PDB structures              | `data/input/*.pdb`            | `--input-dir` / `--pdb-dir`  |
| Protein interaction table   | `data/input/pairs/*.tsv`      | positional arg to `ppinsight parse` |
| Parsed pairs file           | `data/input/pairs/pairs.csv`  | `-o` on `ppinsight parse`    |

### Outputs — where to find it

| What                        | Default location                         | CLI flag to override         |
|-----------------------------|------------------------------------------|------------------------------|
| Docking outputs             | `data/output/<engine>_runs/`             | `--output-root`              |
| Unified scores file         | `data/output/scores/scores.tsv`          | `-o` on `ppinsight collect`  |
| Run record (provenance)     | `data/output/scores/scores.tsv.provenance.json` | (accompanies scores file) |
| Batch results table         | `data/output/scores/batch_results.csv`   | `-o` on `ppinsight batch`    |
| Plots and figures           | shown interactively; saved when you add `-o <path>` | `-o` on `ppinsight compare`  |

## Pairs file format

The pairs file is the input to `ppinsight batch` (and optionally
`ppinsight collect --pairs`).  You can create one in two ways:

1. **Automatically** — run `ppinsight parse` on a protein interaction table.
2. **By hand** — create a CSV or TSV with the columns below.

### Required columns

| Column       | Type   | Description |
|--------------|--------|-------------|
| `proteinA`   | string | Name or identifier of the first protein (must match a PDB filename stem in `data/input/`, e.g. `2UUY_rec` matches `2UUY_rec.pdb`) |
| `proteinB`   | string | Name or identifier of the second protein (same matching rule) |

### Optional columns

| Column       | Type   | Values / description |
|--------------|--------|----------------------|
| `label`      | string | `interaction` or `non-interaction`.  Used by `--classify` and ROC analysis to separate true positives from true negatives.  If omitted, defaults to empty. |
| `family`     | string | Protein family grouping (e.g. `Ephrin`, `FGFR`).  Carried through to outputs for filtering. |
| `references` | string | Citation numbers from the source table (e.g. `5,8`).  Informational only. |

Column names are **case-insensitive** and leading/trailing whitespace is
stripped, so `ProteinA`, `proteinA`, and `proteinA ` all work.

### Minimal example (`pairs.csv`)

```csv
proteinA,proteinB,label
2UUY_rec,2UUY_lig,interaction
EPHA1,ERBB2,interaction
EPHA1,FGFR3,non-interaction
```

### Full example (`pairs.tsv`)

```tsv
proteinA	proteinB	label	family	references
EPHA1	ERBB2	interaction	Ephrin	6
EPHA1	FGFR2	interaction	Ephrin	7
EPHA1	FGFR3	non-interaction	Ephrin	
EPHA2	MET	interaction	Ephrin	10,12
```

### Notes

- **File extension determines delimiter**: `.tsv` → tab-separated,
  `.csv` → comma-separated.  The pipeline detects this automatically.
- **Protein names must match PDB filenames** in `data/input/` (or
  wherever `--pdb-dir` points).  For example, a row with
  `proteinA=2UUY_rec` expects a file named `2UUY_rec.pdb`.
- If you skip the `label` column, classification (`--classify`) and
  ROC curves will not work, but docking and score collection will
  proceed normally.
- `ppinsight parse` produces this exact format from an
  protein interaction table.  If your data is already
  flat (one row per pair), just create the CSV/TSV directly — no need
  to run `parse`.

## Typical workflow

```bash
# 1. Place your PDB files in data/input/
cp my_structures/*.pdb data/input/

# 2. Place your protein interaction table in data/input/pairs/, then parse it
ppinsight parse data/input/pairs/my_proteins.tsv \
    -o data/input/pairs/pairs.csv --stats

# 3. Batch-dock (reads PDBs from data/input/, writes to data/output/)
ppinsight batch data/input/pairs/pairs.csv \
    --engines lightdock haddock \
    --pdb-dir data/input/ \
    --dry-run            # verify first, remove --dry-run for real run

# 4. Collect scores (writes to data/output/scores/ by default)
ppinsight collect data/output/lightdock_runs data/output/haddock_runs \
    --pairs data/input/pairs/pairs.csv

# 5. Compare / classify / plot
ppinsight compare data/output/scores/scores.tsv -m dockq --classify
ppinsight compare data/output/scores/scores.tsv -m score \
    --plot-type roc -o data/output/plots/roc.png
```
