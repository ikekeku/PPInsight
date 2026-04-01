# PPInsight

A Python toolkit for benchmarking computational protein–protein interaction
(PPI) prediction models.  PPInsight automates protein data retrieval, runs
multiple docking engines, collects scores into a unified format, and
produces comparative visualisations.

![Workflow schematic](docs/images/CSE583_PPInsight-Schematic.png)

## Contents

- [Install & quick start](#install--quick-start)
- [CLI reference](#cli-reference)
- [Project structure](#project-structure)
- [Team members](#team-members)
- [License](#license)
- [Resources](#resources)

---

## Install & quick start

### One-command setup

```bash
bash setup.sh          # creates conda env, installs PyRosetta, pip install -e .
conda activate ppinsight
```

> `setup.sh` downloads PyRosetta (~1.5 GB) automatically.  To skip that
> step (e.g. if you only need HADDOCK / LightDock), run
> `bash setup.sh --no-rosetta`.

<details>
<summary>Manual step-by-step setup</summary>

```bash
conda env create -f environment.yml
conda activate ppinsight
python -c "import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()"
pip install -e .
```
</details>

### Run tests

```bash
python -m pytest tests/ -v
```

---

## CLI reference

All tools are available both as standalone commands **and** as subcommands
of the `ppinsight` umbrella CLI:

```bash
ppinsight <command> [options]
# or equivalently:
<command> [options]
```

### 1. Fetch protein data

```bash
ppinsight fetch P69905 P68871 --fasta sequences.fasta --csv metadata.csv
# standalone:
protein_fetch P69905 P68871 --fasta sequences.fasta --csv metadata.csv
```

### 2. Run docking pipelines

```bash
# HADDOCK (stage config + optional execution)
ppinsight haddock 2UUY_rec 2UUY_lig --runname example_2UUY
ppinsight haddock 2UUY_rec 2UUY_lig --run              # also execute

# LightDock
ppinsight lightdock 2UUY_rec 2UUY_lig --steps 100 --generate

# Rosetta (requires PyRosetta)
ppinsight rosetta 2UUY_rec 2UUY_lig --n-runs 5
```

All docking CLIs accept `--input-dir <dir>` to restrict file search to a
specific directory instead of the whole repository.

### 3. Collect scores

Aggregate docking outputs into a single unified scores file:

```bash
ppinsight collect examples/haddock3/run1-test examples/lightdock/simulation \
    --pair e2aP:hpr -o scores.tsv
# standalone:
collect_scores examples/haddock3/run1-test -o scores.tsv
```

### 4. Visualise & compare

```bash
# Violin plot (default) — always prints a tabular summary first
ppinsight compare scores.tsv --metric dockq

# Filter by protein pair and save to PNG
ppinsight compare scores.tsv --metric dockq --pair 2UUY_rec:2UUY_lig -o plot.png

# Tabular summary only (no plot)
ppinsight compare scores.tsv --metric dockq --table

# CAPRI quality classification (high/medium/acceptable/incorrect)
ppinsight compare scores.tsv --capri-quality
ppinsight compare scores.tsv --plot-type quality_bar -o quality.png

# Compare per-model files side by side
ppinsight compare haddock_scores.tsv rosetta_scores.csv \
    --metric score --names HADDOCK Rosetta

# Discover what's in a scores file
ppinsight compare scores.tsv --list-metrics
ppinsight compare scores.tsv --list-pairs
```

### 5. Parse interaction tables & batch-dock

```bash
# Parse an annotation table into a flat pairs file
ppinsight parse data/input/pairs/RTK_Interactome.tsv -o data/input/pairs/pairs.csv --stats

# Batch-dock all pairs (dry-run first, then for real)
ppinsight batch data/input/pairs/pairs.csv --engines lightdock haddock --dry-run
ppinsight batch data/input/pairs/pairs.csv --engines lightdock --pdb-dir data/input/

# Limit to first N pairs for a quick test
ppinsight batch data/input/pairs/pairs.csv --engines lightdock --limit 5
```

### 6. Evaluate docking quality (DockQ)

```bash
# Score a docked model against a native structure
ppinsight quality model.pdb native.pdb

# Score every model in a directory
ppinsight quality --model-dir docked_models/ --native native.pdb -o quality.tsv
```

> Requires the optional `quality` extra: `pip install ppinsight[quality]`

---

## Project structure

```
src/ppinsight/
  cli.py               # Umbrella ppinsight CLI
  protein_fetch.py      # UniProt / PDB data retrieval
  pdb_to_lightdock.py   # LightDock docking wrapper
  pdb_to_haddock.py     # HADDOCK3 staging & execution
  pdb_to_rosetta.py     # PyRosetta docking wrapper
  collect_scores.py     # Score aggregator (unified TSV/CSV)
  visualizer.py         # Plotting & compare_scores CLI
  batch_dock.py         # Batch docking for all pairs
  parse_pairs.py        # Annotation table → pairs file
  quality.py            # DockQ quality evaluation
  utils.py              # Shared utilities (path resolution)
  rosetta/              # PyRosetta pipeline internals
data/                   # User I/O: input PDBs & pairs → output docking runs
tests/                  # pytest test suite
examples/               # Sample data, tutorials, and test fixtures
docs/                   # Specs, use cases, and slide deck
```

### User input & output directories

The **`data/`** directory at the repository root is the default working
directory for real pipeline runs.  Drop your files here; the pipeline
reads from `data/input/` and writes to `data/output/`.

> `examples/` is for sample data, tutorials, and test fixtures only.

```
data/
├── input/
│   ├── *.pdb              # Drop receptor / ligand PDB structures here
│   └── pairs/             # Drop batch CSV / TSV files here
│       ├── pairs.csv      #   (output of `ppinsight parse`)
│       └── <table>.tsv    #   (your annotation / interaction table)
│
└── output/
    ├── haddock_runs/      # HADDOCK docking outputs
    ├── lightdock_runs/    # LightDock docking outputs
    └── rosetta_runs/      # Rosetta docking outputs
```

| What you have             | Where to put it              | CLI flag to override         |
|---------------------------|------------------------------|------------------------------|
| PDB structures            | `data/input/`                | `--input-dir` / `--pdb-dir`  |
| Interaction table (TSV)   | `data/input/pairs/`          | positional arg               |
| Parsed pairs file         | `data/input/pairs/`          | `-o` on `ppinsight parse`    |
| Docking outputs           | `data/output/`               | `--output-root`              |
| Unified scores file       | working directory            | `-o` on `ppinsight collect`  |
| Plots / figures           | working directory            | `-o` on `ppinsight compare`  |

> **Tip:** `data/input/pairs/` is the drop folder for batch-mode CSV/TSV
> files.  Place your annotation table there, then run
> `ppinsight parse` → `ppinsight batch` → `ppinsight collect`.
> See [`data/README.md`](data/README.md) for a full walkthrough.

### Notes

- For HPC runs, pass an absolute `--output-root` or `--output-dir` to
  place outputs on a shared filesystem (e.g. `/gscratch/...`).
- External tools (HADDOCK, LightDock, Rosetta/PyRosetta) are required for
  full pipeline execution.  The HADDOCK helper supports container
  execution (`--container docker|apptainer`).
- Tests monkeypatch external calls; run `python -m pytest` to execute the
  full suite.

---

## Team members

| Name | Contributions |
|------|--------------|
| **Ike Keku** | Co-developed PPI predictor pipelines (`pdb_to_haddock`, `pdb_to_lightdock`), oversaw general workflow |
| **Rita Kamenetskiy** | Co-developed protein data fetching (`protein_fetch`), packaging, CI tests |
| **Fiona McLary** | Developed visualiser module |
| **Walter Avila** | Co-developed `protein_fetch` and corresponding tests |
| **Maya Gatt Harari** | Developed Rosetta integration (`ppinsight/rosetta`) |

---

## License

See [LICENSE](LICENSE).

## Resources

- [Biopython Entrez guide](https://biopython.org/wiki/Entrez)
- [RCSB API docs](https://data.rcsb.org/)
- [DockQ](https://github.com/bjornwallner/DockQ)
- [HADDOCK docs](https://wenmr.science.uu.nl/haddock2.4/)
- [LightDock](https://github.com/lightdock/lightdock)
- [Rosetta / PyRosetta](https://rosettacommons.org)

### CAPRI quality classification thresholds

The `--capri-quality` flag and `quality_bar` plot classify docked
predictions into four tiers (high / medium / acceptable / incorrect)
using published community thresholds:

| Source | DOI | Used for |
|--------|-----|----------|
| Lensink MF, Velankar S, Wodak SJ (2016). *Proteins* 84(S1):323-348 | [10.1002/prot.25007](https://doi.org/10.1002/prot.25007) | fnat + Lrms + irms thresholds (Table 3) |
| Basu S, Wallner B (2016). *PLoS ONE* 11(8):e0161879 | [10.1371/journal.pone.0161879](https://doi.org/10.1371/journal.pone.0161879) | DockQ fallback thresholds |
