# Functional Specification

## Background
This project builds a tool to benchmark protein–protein interaction (PPI) prediction models. Experimental PPI determination is slow and costly, deterring researchers from finding the best fit for evaluating protein-protein interactions. Many computational models exist, but outputs are hard to compare. This system retrieves protein sequences/structures, runs multiple PPI models, and visualizes their scores side-by-side. This will allow for those interested in PPI interactions to not only evaluate different PPI models, but also evaluate any model discrepancies. 

## User profile
- **Primary users:** Primary users will include graduate students and researchers in biology, bioengineering, or computer science. Mainly these will be investigators who are involved in the field of protein-protein interactions either as their main field of study or a sub category in which they would not have a lot of expertise specifically with protein-protein interactions but are interested in understanding it from a modeling perspective. 
- **Knowledge:** Basic molecular biology (what a protein is, why interactions matter) and  basic computing (can run Python scripts, use command line, browse the web). No deep ML or structural-biology expertise required. 

## Data sources
**Protein sequences/structures:** FASTA / PDB / mmCIF files from NCBI, UniProt, or AlphaFold.
- **FASTA:** A FASTA protein file is text-based format used to represent protein sequences. The data structure usually starts with a ">" sign followed by a name identification of the protein sequence with an optional description. Example: >crab_rabit ALPHA CRYSTALLIN B CHAIN (ALPHA(B)-CRYSTALLIN). The rest of the data is followed by the acual protein sequence which is represented through single-letter amino acid codes.Example: MDIAIHHPWIRRPFFPFHSPSRLFDQFFGEHLLESDLFPTSTSLSPFYLR PPSFLRAPSWIDTGLSEMRLEKDRFSVNLDVKHFSPEELKVKVLGDVIEV HGKHEERQDEHGFISREFHRKYRIPADVDPLTITSSLSSDGVLTVNGPRK QAPGPERTIPITREEKPAVTAAPKK
- **PDB:** Protein data bank file is in standard text format to store 3D structural data of proteins. 3D structural data is designated through the X,Y,Z position of the atom in 3D space with the designation of the atom type, residue name, chain ID, and occupancy. Below is an example of a small protein molecule, glucagon and its according PDB:
  <img width="641" height="308" alt="image" src="https://github.com/user-attachments/assets/25e138d8-f124-4272-b8a9-37d0a0f35d13" />
  
  Here we see that the first designation of the data is 'ATOM', then the atom serial number, the atom name, a branch indicator if required, residue type, the chain identifier, the residue sequence number, the X, Y, Z coordinate values repectively with the last three fields detailing occupancy, temperature factor (B-factor) and element symbol.
- **mmCIF:** A macromolecular crystallographic information file is a more novel form of storing 3D structural data of macromolecules. It is similar to PDB files but is more flexible and machine-readable. This is a dictionary-based format which utilizes tags and values, making it easier to be read by specific packages. It stores all the original data a PDB file can store but also has more complex and experimental inputs and has no limitation on the number of atoms or chains as its input. 

  
**Model outputs:** 
- **Retrieved Protein Data:** When the software fetches protein data from external databases, it will output: the FASTA and PDB/mmCIF structual files for each protein specified by the user. These data will be returned in a dataframe.

- **Model Execution Outputs:** The PPI models will output standardized artifacts and metrics to support cross-model comparison. The format in which each predictor organizes their output data may differ between models. Roughly, each model will return: DockQ score, RMSD, confidence metrics, protein IDs, predicted complex structures (PDB/mmCIF), and the FASTA file of the aligned protein sequences.

- **Model Execution Outputs:** The tool will generate on-screen and exportable visualizations. The visual will be a bar chart (made using matplotlib) showing each model's score for the input protien pair. Moreover, the plot ouput will include the data table used to produce the graph.

## Use Case 1 - Compare models on one protein pair
- **Objective:** Determine which PPI model gives the best score for a given pair (e.g., proteinA-proteinB).
- **Interaction:** User enters protein names → system fetches FASTA/PDB → user selects 2–3 models → system runs models → writes scores to table → shows plot of scores per model

## Use Case 2 - Compare models in bulk, for many protein pairs
- **Objective:** Determine which PPI model gives the best score for many pairs (e.g., proteinA-proteinB).
- **Interaction:** User uploads CSV of pairs of proteins → system fetches FASTA/PDB → user selects 2–3 models → system runs models → writes scores to table → provides statistic of the different models per pair.
  
## Use Case 3 - Curate protein data
- **Objective:** Generate FASTA / PDB files for a batch of proteins.  
- **Interaction:** User uploads CSV of pairs → system fetchs FASTA/PDB file 
