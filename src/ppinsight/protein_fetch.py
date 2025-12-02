import requests
from Bio import SeqIO
from io import StringIO
from Bio.PDB import PDBList
import csv
import os 
import sys



def get_uniprot_data(accession_ids, fasta_file=None, csv_file=None, pdb_dir="pdb_files"):
    """
    Unified pipeline:
    1. Fetch UniProt FASTA sequences for given accession IDs.
    2. Save FASTA file.
    3. Extract structured data and save to CSV.
    4. Retrieve PDB IDs from UniProt.
    5. Download the first PDB file for each accession.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/"
    sequences_data = []
    structured_data = []

    # Ensure PDB directory exists
    os.makedirs(pdb_dir, exist_ok=True)

    for accession_id in accession_ids:
        # Fetch FASTA
        fasta_url = f"{base_url}{accession_id}.fasta"
        try:
            response = requests.get(fasta_url, timeout=10)
            response.raise_for_status()
            sequences_data.append(response.text)
        except requests.exceptions.RequestException as e:
            raise ValueError(f"Invalid UniProt accession ID: {accession_id}") from e
            raise ValueError(f"Could not fetch FASTA for {accession_id}. Error: {e}")

    full_fasta_string = "".join(sequences_data)

    # Save FASTA file if requested
    if fasta_file and full_fasta_string:
        with open(fasta_file, "w") as f:
            f.write(full_fasta_string)
        print(f"FASTA sequences saved to {fasta_file}")

    # Parse sequences into SeqRecord objects
    records = list(SeqIO.parse(StringIO(full_fasta_string), "fasta")) if full_fasta_string else []

    # Extract structured data
    for entry in records:
        structured_data.append({
            "ID": entry.id,
            "Name": entry.name,
            "Description": entry.description,
            "Sequence Length": len(entry.seq),
            "Sequence": str(entry.seq)
        })

    # Save structured data to CSV if requested
    if csv_file and structured_data:
        with open(csv_file, "w", newline="") as csv_out:
            writer = csv.DictWriter(csv_out, fieldnames=["ID", "Name", "Description", "Sequence Length", "Sequence"])
            writer.writeheader()
            writer.writerows(structured_data)
        print(f"Structured data saved to {csv_file}")

    # Fetch PDB IDs and download first PDB file for each accession
    pdbl = PDBList()
    pdb_info = {}
    for accession_id in accession_ids:
        json_url = f"{base_url}{accession_id}.json"
        try:
            response = requests.get(json_url, timeout=10)
            response.raise_for_status()
            data = response.json()
            pdb_ids = [xref.get("id") for xref in data.get("uniProtKBCrossReferences", []) if xref.get("database") == "PDB"]
            if pdb_ids:
                first_pdb_id = pdb_ids[0]
                pdbl.retrieve_pdb_file(first_pdb_id, pdir=pdb_dir, file_format="pdb")
                pdb_info[accession_id] = first_pdb_id
                print(f"Downloaded PDB file for {accession_id}: {first_pdb_id}")
            else:
                pdb_info[accession_id] = None
                print(f"No PDB IDs found for {accession_id}")
        except requests.exceptions.RequestException as e:
            print(f"Warning: Could not fetch PDB info for {accession_id}. Error: {e}", file=sys.stderr)
            pdb_info[accession_id] = None

    return structured_data, pdb_info