"""
protein_fetch.py
----------------

Tools for retrieving protein sequence and structure data from UniProt and
the Protein Data Bank (PDB).

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

Example (CLI):
    protein_fetch P69905 P68871
    protein_fetch P69905 P68871 --fasta hemoglobin.fasta --csv hemoglobin.csv

Example (API):
    >>> from ppinsight.protein_fetch import get_uniprot_data
    >>> data, pdb_info = get_uniprot_data(
    ...     ["P69905", "P68871"],
    ...     fasta_file="hemoglobin.fasta",
    ...     csv_file="hemoglobin.csv",
    ...     pdb_dir="pdb_files"
    ... )
"""

import argparse
import csv
import os
import shutil
import sys
from io import StringIO

import requests
from Bio import SeqIO
from Bio.PDB import PDBList


def _alias_pdb_path(accession_id: str, pdb_dir: str) -> str:
    """Return the user-facing PDB filename for an accession.

    Fetch writes accession-based aliases so the same identifiers can be used
    directly in ``proteinA`` / ``proteinB`` pairs files and docking commands.
    """
    return os.path.join(pdb_dir, f"{accession_id}.pdb")


def _copy_to_accession_alias(
    raw_path: str,
    accession_id: str,
    pdb_dir: str,
) -> str | None:
    """Copy a downloaded PDB artifact to an accession-named ``.pdb`` file.

    Returns ``None`` when no readable raw artifact exists (for example when
    upstream retrieval reports an ID but does not materialize a local file).
    """
    if not raw_path:
        return None

    alias_path = _alias_pdb_path(accession_id, pdb_dir)
    raw_abs = os.path.abspath(raw_path)
    alias_abs = os.path.abspath(alias_path)

    if not os.path.exists(raw_abs):
        return None

    if raw_abs != alias_abs:
        shutil.copyfile(raw_abs, alias_abs)
        raw_name = os.path.basename(raw_abs).lower()
        if raw_name.startswith("pdb") and raw_name.endswith(".ent"):
            try:
                os.remove(raw_abs)
            except OSError:
                pass

    return alias_path


def get_uniprot_data(accession_ids, fasta_file=None, csv_file=None,
                     pdb_dir="pdb_files"):
    """
    Fetch protein sequence and structure data from UniProt and PDB.

    Args:
        accession_ids (list[str]): List of UniProt accession IDs to process.
        fasta_file (str, optional): Path to save the combined FASTA sequences.
        csv_file (str, optional): Path to save structured sequence data as CSV.
        pdb_dir (str, optional): Directory to store downloaded PDB files.
        Defaults to "pdb_files".

    Returns:
        tuple:
            structured_data (list[dict]): Parsed sequence
            records with metadata.
            pdb_info (dict): Mapping of accession IDs to their first PDB ID
            (or None if unavailable).

    Workflow:
        1. Fetch FASTA sequences from UniProt.
        2. Optionally save FASTA to file.
        3. Parse sequences into structured records.
        4. Optionally save structured data to CSV.
        5. Retrieve PDB cross-references and download the
        first PDB file for each accession.
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
            raise ValueError(
                f"Invalid UniProt accession ID: {accession_id}"
                ) from request_error

    full_fasta_string = "".join(sequences_data)

    # Save FASTA file if requested
    if fasta_file and full_fasta_string:
        with open(fasta_file, "w", encoding="utf-8") as fasta_out:
            fasta_out.write(full_fasta_string)
        print(f"FASTA sequences saved to {fasta_file}")

    # Parse sequences into SeqRecord objects
    records = list(
        SeqIO.parse(StringIO(full_fasta_string), "fasta")
        ) if full_fasta_string else []

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
            writer = csv.DictWriter(
                csv_out, fieldnames=[
                    "ID", "Name", "Description", "Sequence Length", "Sequence"]
                    )
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
            pdb_ids = [xref.get("id") for xref in data.get(
                "uniProtKBCrossReferences", []) if xref.get("database"
                                                            ) == "PDB"]
            if pdb_ids:
                first_pdb_id = pdb_ids[0]
                raw_pdb_path = pdbl.retrieve_pdb_file(
                    first_pdb_id, pdir=pdb_dir, file_format="pdb")
                alias_path = _copy_to_accession_alias(
                    raw_pdb_path, accession_id, pdb_dir,
                )
                pdb_info[accession_id] = first_pdb_id
                if alias_path:
                    print(
                        f"Downloaded PDB file for {accession_id}: "
                        f"{first_pdb_id} -> {alias_path}"
                    )
                else:
                    print(
                        f"PDB ID found for {accession_id}: {first_pdb_id}, "
                        "but no local file was downloaded"
                    )
            else:
                pdb_info[accession_id] = None
                print(f"No PDB IDs found for {accession_id}")
        except requests.exceptions.RequestException as request_error:
            print(f"Warning: Could not fetch PDB info for {accession_id}."
                  f"Error: {request_error}",
                  file=sys.stderr,
                 )
            pdb_info[accession_id] = None
            raise ValueError(
                f"Warning: Could not fetch PDB info for: {accession_id}"
            ) from request_error

    return structured_data, pdb_info


# ── CLI ──────────────────────────────────────────────────────────────

def main(argv=None):
    """CLI entrypoint for fetching protein data from UniProt.

    Usage matches the style of the other PPInsight scripts::

        protein_fetch P69905 P68871
        protein_fetch P69905 P68871 --fasta hemoglobin.fasta --csv hemoglobin.csv
    """
    parser = argparse.ArgumentParser(
        description="Fetch protein sequence and structure data from UniProt / PDB."
    )
    parser.add_argument(
        "accessions",
        nargs="+",
        help=(
            "One or more UniProt accession IDs (e.g. P69905 P68871).  "
            "The pipeline fetches sequence data and resolves PDB structures "
            "for each accession."
        ),
    )
    parser.add_argument(
        "--fasta",
        default=None,
        help=(
            "Save combined FASTA sequences to this file.  Useful for "
            "sequence alignment or as input to other bioinformatics tools.  "
            "Omit if you only need PDB structures."
        ),
    )
    parser.add_argument(
        "--csv",
        default=None,
        help=(
            "Save structured sequence metadata (ID, length, organism, etc.) "
            "to this CSV file.  Useful for record-keeping in large benchmarks."
        ),
    )
    parser.add_argument(
        "--pdb-dir",
        default="data/input",
        help=(
            "Directory for downloaded PDB files (default: data/input).  "
            "Fetch saves accession-named .pdb files here (for example "
            "P69905.pdb) and cleans up the raw Biopython download name so "
            "the same accession can be used directly in proteinA/proteinB."
        ),
    )

    args = parser.parse_args(argv)

    structured_data, pdb_info = get_uniprot_data(
        accession_ids=args.accessions,
        fasta_file=args.fasta,
        csv_file=args.csv,
        pdb_dir=args.pdb_dir,
    )

    # Print a brief summary
    print(f"\nFetched {len(structured_data)} protein(s):")
    for entry in structured_data:
        pdb_id = pdb_info.get(entry["ID"].split("|")[1]
                              if "|" in entry["ID"]
                              else entry["ID"], "—")
        print(f"  {entry['ID']:30s}  {entry['Sequence Length']:>5d} aa  PDB: {pdb_id}")

    ready = [acc for acc, pdb_id in pdb_info.items() if pdb_id]
    if ready:
        print("\nPair-ready PDB names:")
        for accession_id in ready:
            print(f"  {accession_id} -> {_alias_pdb_path(accession_id, args.pdb_dir)}")
        print(
            "Use these accession stems directly in proteinA/proteinB, for "
            "example: P69905:P68871"
        )


if __name__ == "__main__":
    main()
