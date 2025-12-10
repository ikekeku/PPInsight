# PPInsight
This project seeks to streamline the evaluation of computational protein-protein interaction (PPI) prediction models. Because experimental PPI determination is time-consuming and costly, computational approaches hold promise for large-scale interactome studies. Existing models frequently lack standardised evaluation measures and instruments for comparison. To close this gap, we propose creating a Python-based benchmarking and visualization tool that automates the extraction of reference interaction data, runs various predictive models, and produces comparative performance graphs. This approach will allow researchers to easily evaluate model correctness and consistency.
Tasks include developing a Python tool to automate protein sequence retrieval, convert data into structural formats, and visualize interaction scores of PPI prediction models.

![Schematic detailling workflow for PPInsight project, grouped by task](docs/images/CSE583_PPInsight-Schematic.png)

# Team Members and Contributions
- Ike Keku
- Rita Kamenetskiy: Co-developed protein data fetching module ('protein_fetch.py'), organized and functionalized the packaging, and functionalized the continuous integration tests
- Fiona McLary: Developed visualizer module 
- Walter Avila: Co-developed protein data fetching module (`protein_fetch.py`) and corresponding tests
- Maya Gatt Harari

## Task 1 — fetch proteins by search term and save files
Implemented in `src/ppinsight/protein_fetch.py` via:

`get_uniprot_data(accession_ids, fasta_file=None, csv_file=None, pdb_dir="fetched_data")`

This function retrieves protein sequences and structures from UniProt using a user-inputted list of accession IDs. It:

* Fetches FASTA sequences and saves them to a file
* Extracts structured metadata (ID, name, description, sequence length, sequence) and saves to CSV
* Queries UniProt JSON records for linked PDB IDs and automatically downloads the first available structure per protein into `pdb_files/`

**Returns:**

* `structured_data`: list of dictionaries with fields
  `ID`, `Name`, `Description`, `Sequence Length`, `Sequence`
* `pdb_info`: dictionary mapping
  `accession_id → PDB_ID` (or `None` if no structure exists)

**Testing:**

Pytest tests cover successful runs, correct sequence length and structure availability, and invalid accession ID handling.

---

## Task 2 — run PPI predictors on file pairs and store scores  
Once you’ve retrieved the protein files in Task 1, here we swap from “getting the data” to “using models” to predict how strongly two proteins interact and then save those predictions to compare models later. You don’t need to code every model from scratch — instead, make small Python functions (“wrappers”) that take two protein files and return the model’s predicted score. This helps keep your code organized and reusable. This task wraps 2–3 PPI prediction approaches and produces a single `scores.tsv` (or CSV) in the same repo. For each protein pair and model run, write one or more score rows with clear headers (`proteinA, proteinB, model, score_type, score_value, output_path, timestamp`). Start with three complementary predictors: can include one of each or some other combination of a sequence-based model (fast local or embedding-based), a docking server (ClusPro or HADDOCK) for structural docking, and a structure-scoring tool (DockQ) for pose quality.

Selected Models (of interest):
- **Haddock** (https://github.com/haddocking; https://pmc.ncbi.nlm.nih.gov/articles/PMC3966529/)
- **Rosetta** (https://github.com/RosettaCommons/rosetta; https://rosettacommons.org)
- **LightDock** (https://github.com/lightdock/lightdock; https://academic.oup.com/bioinformatics/article/34/1/49/4103399)

Key implementation notes:
- For each protein pair, you’ll run models that predict how well two proteins fit or “dock” together — meaning how their 3D structures might align to form a stable complex. The output is a numerical score that tells you how good that fit is.
- Make lightweight wrappers so `score_pair(a_files, b_files, model_name)` is non-repetitive and logs raw outputs. See ClusPro help and usage: [https://cluspro.org/help.php](https://cluspro.org/help.php) and HADDOCK docs: [https://wenmr.science.uu.nl/haddock2.4/](https://wenmr.science.uu.nl/haddock2.4/).

**CAPRI Metrics** (https://www.sciencedirect.com/science/article/pii/S0022283624001359, https://link.springer.com/article/10.1186/s12859-024-05991-4)
 - **Docking quality score (DockQ)** – a value from 0 to 1 that rates how well two proteins’ shapes align after docking. Higher = better interaction. [https://github.com/bjornwallner/DockQ](https://github.com/bjornwallner/DockQ)  
 - **RMSD (Root Mean Square Deviation)** – measures how far apart the predicted protein complex is from a reference complex. Lower = better accuracy. [https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions](https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions)  
 - **Fraction of native contacts (Fnat score)** - a measure of the accuracy of a predicted protein-protein interaction model. It is calculated by dividing the number of correctly predicted residue contacts in the docked model by the total number of contacts in the original, native complex. (https://en.wikipedia.org/wiki/Native_contact)

**Other values (may be normalized)**, such as:
 - **E-value** – used in sequence-based predictions (like BLAST). It estimates how likely a match happened by chance. Lower = more significant. [https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html](https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html)

You’ll want to save multiple types of scores (e.g., DockQ, RMSD, E-value) for each prediction so that later you can visualize and compare models across these different performance metrics.
Store multiple `score_type` values per run (e.g., `dockq`, `rmsd`, `evalue`) so visualizations can pick which metric to display.

Testing and reproducibility:
- Create simple test files (like three known interacting protein pairs in a pairs.csv file).
- Write “test stubs” — mini functions that pretend to run a model but always return the same score. These let you test your plotting and data-handling code without needing to re-run expensive real models.

---

## Task 3 — plot interaction scores (D3 visualizations)  
This task turns the `scores.tsv` table into interactive visuals. For each metric type (one chart per metric—e.g., DOCKQ, RMSD), build a D3 page where x = model name and y = score value; each protein pair is shown as a series of points/boxes, with tooltips linking to raw outputs. The default view will show a single protein pair (dropdown to switch). Possible immplenetations include allowing grouped views (multiple pairs) or aggregated summaries (boxplot/violin per model).

Key implementation notes:
- Convert `scores.tsv` → `scores.json` via a small script; D3 reads the JSON. Use D3 v7: [https://d3js.org/](https://d3js.org/).  
- UX (user experience): dropdown for protein pair, legend for models, sort controls (by median score), hover tooltip with `score_value`, `score_type`, and `output_path`. Export PNG/SVG for posters.  
- Keep each chart focused: one metric per page (easier comparison and legend clarity).

Testing and examples:
- Supply `examples/plots/*.html` and a small server (e.g., `python -m http.server`) for local review. Include snapshot SVGs for visual regression if desired.

---

### Repo expectations  
Put code and examples in this GitHub repo and update `README.md`, `.ipynb` files (using comments), and `examples/` (pairs.csv, example outputs) as needed. 
### Resources (update as found/needed)
Biopython Entrez guide: [https://biopython.org/wiki/Entrez](https://biopython.org/wiki/Entrez). RCSB API docs: [https://data.rcsb.org/](https://data.rcsb.org/). DockQ: [https://github.com/bjornwallner/DockQ](https://github.com/bjornwallner/DockQ). D3: [https://d3js.org/](https://d3js.org/).

---


