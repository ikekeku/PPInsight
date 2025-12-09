"""
protein_fetch.py
----------------

Tools for retrieving protein sequence and structure data from UniProt and the Protein Data Bank (PDB).

The module provides a single pipeline that:
- Fetches FASTA sequences for UniProt accession IDs
- Optionally saves the raw FASTA output to a file
- Parses sequences into structured records and writes them to CSV
- Retrieves PDB cross-references from UniProt
- Downloads the first available PDB file for each accession

Dependencies:
- requests: HTTP requests to UniProt REST API
- biopython: sequence parsing (SeqIO) and PDB handling (Bio.PDB)
- csv, os, sys: file and system utilities

Example:
    >>> from protein_fetch import get_uniprot_data
    >>> data, pdb_info = get_uniprot_data(
    ...     ["P69905", "P68871"],
    ...     fasta_file="hemoglobin.fasta",
    ...     csv_file="hemoglobin.csv",
    ...     pdb_dir="pdb_files"
    ... )
"""

import csv
import os
import sys
from io import StringIO

import requests
from Bio import SeqIO
from Bio.PDB import PDBList



def get_uniprot_data(accession_ids, fasta_file=None, csv_file=None, pdb_dir="pdb_files"):
    """
    Fetch protein sequence and structure data from UniProt and PDB.

    Args:
        accession_ids (list[str]): List of UniProt accession IDs to process.
        fasta_file (str, optional): Path to save the combined FASTA sequences.
        csv_file (str, optional): Path to save structured sequence data as CSV.
        pdb_dir (str, optional): Directory to store downloaded PDB files. Defaults to "pdb_files".

    Returns:
        tuple:
            structured_data (list[dict]): Parsed sequence records with metadata.
            pdb_info (dict): Mapping of accession IDs to their first PDB ID (or None if unavailable).

    Workflow:
        1. Fetch FASTA sequences from UniProt.
        2. Optionally save FASTA to file.
        3. Parse sequences into structured records.
        4. Optionally save structured data to CSV.
        5. Retrieve PDB cross-references and download the first PDB file for each accession.
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
        except requests.exceptions.RequestException as request_error:
            raise ValueError(f"Invalid UniProt accession ID: {accession_id}") from request_error

    full_fasta_string = "".join(sequences_data)

    # Save FASTA file if requested
    if fasta_file and full_fasta_string:
        with open(fasta_file, "w", encoding ="utf-8") as fasta_out:
            fasta_out.write(full_fasta_string)
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
        with open(csv_file, "w", newline="", encoding="utf-8") as csv_out:
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
        except requests.exceptions.RequestException as request_error:
            print(f"Warning: Could not fetch PDB info for {accession_id}. Error: {request_error}", file=sys.stderr)
            pdb_info[accession_id] = None
            raise ValueError(f"Warning: Could not fetch PDB info for: {accession_id}") from e

    return structured_data, pdb_info
