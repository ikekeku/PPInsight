My name is Rita Kamenetskiy and I am a second year bioengineering student. 

# PPInsight
This project seeks to streamline the evaluation of computational protein-protein interaction (PPI) prediction models. Because experimental PPI determination is time-consuming and costly, computational approaches hold promise for large-scale interactome studies. Existing models frequently lack standardised evaluation measures and instruments for comparison. To close this gap, we propose creating a Python-based benchmarking and visualization tool that automates the extraction of reference interaction data, runs various predictive models, and produces comparative performance graphs. This approach will allow researchers to easily evaluate model correctness and consistency.
Tasks include developing a Python tool to automate protein sequence retrieval, convert data into structural formats, and visualize interaction scores of PPI prediction models.

# Team Member
- Ike Keku
- ...
- ...
- ...
- ...

---

## Task 1 — fetch proteins by search term and save files  
This task builds a single function (e.g., `def create_FASTA(protein_name, species, ...)`) that takes a search string for NCBI Entrez (gene/protein name, organism, etc.) and returns matched protein sequence and structural data written to disk in FASTA, CSV/TSV (metadata), and/or PDB/mmCIF when available. The function should prefer reviewed Swiss-Prot entries, let the user choose among ambiguous hits, and write files using a consistent naming scheme (e.g., `GENEID.fasta`, `GENEID.meta.tsv`, `GENEID.pdb`). Handle rate limits and require `Entrez.email` and optional API key.

Key implementation notes:
- Use Biopython Entrez for searching and efetch. See the Biopython docs: [https://biopython.org/docs/1.76/api/Bio.Entrez.html](https://biopython.org/docs/1.76/api/Bio.Entrez.html).  
- For structures use RCSB APIs or mmCIF downloads: [https://data.rcsb.org/](https://data.rcsb.org/). AlphaFold predicted structures are available too: [https://alphafold.ebi.ac.uk/](https://alphafold.ebi.ac.uk/).  
- Outputs: FASTA file, `{id}.meta.tsv` or `{id}.meta.csv` (columns: id, name, organism, source_db, accession, length), and `{id}.pdb` or `{id}.cif` when present.

Testing and edge cases:
- Test ambiguous queries (checking the function can handle vague search terms and returns a list of possible matching proteins instead of failing).
- Test missing structures (verifying that if no PDB or structure is found, the code still saves the sequence (FASTA) and metadata without error).
- Add unit tests that mock Entrez responses (instead of always sending real requests to Biopython’s Entrez API during testing, write tests that simulate Entrez returning specific results or errors. That way the tests run quickly, reliably, and offline, without relying on the live API).

---

## Task 2 — run PPI predictors on file pairs and store scores  
Once you’ve retrieved the protein files in Task 1, here we swap from “getting the data” to “using models” to predict how strongly two proteins interact and then save those predictions to compare models later. You don’t need to code every model from scratch — instead, make small Python functions (“wrappers”) that take two protein files and return the model’s predicted score. This helps keep your code organized and reusable. This task wraps 2–3 PPI prediction approaches and produces a single `scores.tsv` (or CSV) in the same repo. For each protein pair and model run, write one or more score rows with clear headers (`proteinA, proteinB, model, score_type, score_value, output_path, timestamp`). Start with three complementary predictors: can include one of each or some other combination of a sequence-based model (fast local or embedding-based), a docking server (ClusPro or HADDOCK) for structural docking, and a structure-scoring tool (DockQ) for pose quality.

Key implementation notes:
- For each protein pair, you’ll run models that predict how well two proteins fit or “dock” together — meaning how their 3D structures might align to form a stable complex. The output is a numerical score that tells you how good that fit is.
- Make lightweight wrappers so `score_pair(a_files, b_files, model_name)` is non-repetitive and logs raw outputs. See ClusPro help and usage: [https://cluspro.org/help.php](https://cluspro.org/help.php) and HADDOCK docs: [https://wenmr.science.uu.nl/haddock2.4/](https://wenmr.science.uu.nl/haddock2.4/).  
- **Docking quality score (DockQ)** – a value from 0 to 1 that rates how well two proteins’ shapes align after docking. Higher = better interaction. [https://github.com/bjornwallner/DockQ](https://github.com/bjornwallner/DockQ)  
- **RMSD (Root Mean Square Deviation)** – measures how far apart the predicted protein complex is from a reference complex. Lower = better accuracy. [https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions](https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions)  
- **E-value** – used in sequence-based predictions (like BLAST). It estimates how likely a match happened by chance. Lower = more significant. [https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html](https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html)
- You’ll want to save multiple types of scores (e.g., DockQ, RMSD, E-value) for each prediction so that later you can visualize and compare models across these different performance metrics.
- Store multiple `score_type` values per run (e.g., `dockq`, `rmsd`, `evalue`) so visualizations can pick which metric to display.

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
