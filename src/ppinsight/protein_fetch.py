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
import re
import shutil
import sys
from io import StringIO

import requests
from Bio import SeqIO
from Bio.PDB import PDBList

_UNIPROT_BASE_URL = "https://rest.uniprot.org/uniprotkb/"
_UNIPROT_SEARCH_URL = f"{_UNIPROT_BASE_URL}search"
_ENTRY_NAME_HUMAN_PATTERN = re.compile(r"^[A-Z0-9-]+_HUMAN$")


# UniProt accession formats (6-char and 10-char forms)
_ACCESSION_PATTERN_6 = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$")
_ACCESSION_PATTERN_10 = re.compile(
    r"^[A-NR-Z][0-9](?:[A-Z0-9]{3}[0-9]){2}$"
)


def _looks_like_uniprot_accession(token: str) -> bool:
    """Return True when *token* matches a canonical UniProt accession."""
    normalized = token.strip().upper()
    return bool(
        _ACCESSION_PATTERN_6.match(normalized)
        or _ACCESSION_PATTERN_10.match(normalized)
    )


def _looks_like_human_entry_name(token: str) -> bool:
    """Return True when *token* looks like an exact UniProt *_HUMAN ID."""
    return bool(_ENTRY_NAME_HUMAN_PATTERN.match(token.strip().upper()))


def _query_uniprot_search(
    query: str,
    *,
    size: int,
    fields: str,
) -> list[dict]:
    """Execute a UniProt search query and return the JSON result list."""
    try:
        response = requests.get(
            _UNIPROT_SEARCH_URL,
            params={
                "query": query,
                "format": "json",
                "size": size,
                "fields": fields,
            },
            timeout=10,
        )
        response.raise_for_status()
        payload = response.json()
    except (requests.exceptions.RequestException, ValueError):
        return []

    return list(payload.get("results", []))


def _query_uniprot_primary_accession(query: str) -> str | None:
    """Return the first matching primary accession for a UniProt query."""
    rows = _query_uniprot_search(
        query,
        size=1,
        fields="accession,id,reviewed,organism_name",
    )
    for row in rows:
        accession = row.get("primaryAccession")
        if accession:
            return str(accession).upper()
    return None


def _extract_protein_name(result_row: dict) -> str:
    """Return the recommended/submitted protein name from a UniProt row."""
    desc = result_row.get("proteinDescription") or {}
    recommended = desc.get("recommendedName") or {}
    full_name = recommended.get("fullName") or {}
    if isinstance(full_name, dict):
        value = str(full_name.get("value", "")).strip()
        if value:
            return value

    for submitted in desc.get("submissionNames", []) or []:
        sub_full = submitted.get("fullName") or {}
        value = str(sub_full.get("value", "")).strip()
        if value:
            return value

    return ""


def _extract_gene_names(result_row: dict) -> str:
    """Return comma-separated gene symbols/synonyms from a UniProt row."""
    genes = result_row.get("genes") or []
    tokens: list[str] = []

    for gene in genes:
        gene_name = gene.get("geneName") or {}
        primary = str(gene_name.get("value", "")).strip()
        if primary:
            tokens.append(primary)

        for synonym in gene.get("synonyms", []) or []:
            value = str(synonym.get("value", "")).strip()
            if value:
                tokens.append(value)

    deduped: list[str] = []
    for token in tokens:
        if token not in deduped:
            deduped.append(token)
    return ", ".join(deduped)


def _format_protein_existence(value) -> str:
    """Normalise UniProt protein-existence labels for display."""
    text = str(value or "").strip()
    if not text:
        return ""
    if ":" in text:
        return text.split(":", 1)[1].strip()
    return text


def search_uniprot_matches(term: str, limit: int = 5) -> list[dict[str, str]]:
    """Return top UniProt preview matches for an ambiguous search term."""
    limit = max(1, min(int(limit), 20))
    query = f"{term.strip()} AND organism_id:9606 AND reviewed:true"
    rows = _query_uniprot_search(
        query,
        size=limit,
        fields=(
            "accession,id,protein_name,gene_names,organism_name,length,"
            "protein_existence,annotation_score,reviewed"
        ),
    )

    matches: list[dict[str, str]] = []
    for row in rows:
        matches.append({
            "accession": str(row.get("primaryAccession", "")).strip(),
            "entry_name": str(row.get("uniProtkbId", "")).strip(),
            "protein_name": _extract_protein_name(row),
            "gene_names": _extract_gene_names(row),
            "organism": str(
                (row.get("organism") or {}).get("scientificName", "")
            ).strip(),
            "length": str(row.get("sequence", {}).get("length", "")).strip(),
            "protein_existence": _format_protein_existence(
                row.get("proteinExistence")
            ),
            "annotation_score": str(row.get("annotationScore", "")).strip(),
        })

    return matches


def _candidate_name_queries(identifier: str) -> list[str]:
    """Build robust UniProt search queries for human gene/name inputs."""
    raw = " ".join(identifier.strip().split())
    upper = raw.upper()

    # Accept forms such as "VEGFA human" and "VEGFA_HUMAN".
    token = upper
    if token.endswith(" HUMAN"):
        token = token[:-6].strip()
    if token.endswith("_HUMAN"):
        token = token[:-6].strip()

    gene = re.sub(r"[^A-Z0-9-]", "", token)
    queries: list[str] = []

    if gene:
        queries.extend([
            f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
            f"gene:{gene} AND organism_id:9606 AND reviewed:true",
            f"id:{gene}_HUMAN AND reviewed:true",
            f"gene_exact:{gene} AND reviewed:true",
        ])

    # Keep one permissive fallback for names that are not gene symbols.
    if raw:
        queries.append(f"{raw} AND organism_id:9606 AND reviewed:true")

    # Stable de-duplication.
    deduped: list[str] = []
    for q in queries:
        if q not in deduped:
            deduped.append(q)
    return deduped


def _resolve_to_accession(identifier: str) -> str:
    """Resolve a user-provided identifier into a UniProt accession.

    Accepted forms:
    - canonical accession (e.g. P15692)
    - FASTA-style IDs (e.g. sp|P15692|VEGFA_HUMAN)
    - human gene/name inputs (e.g. VEGFA, VEGFA human)
    """
    raw = identifier.strip()
    if not raw:
        raise ValueError(
            "Invalid UniProt accession ID: ''. Use an accession such as "
            "P15692 or a human gene/name such as VEGFA."
        )

    parts = raw.split("|")
    if len(parts) >= 3 and _looks_like_uniprot_accession(parts[1]):
        return parts[1].upper()

    if _looks_like_uniprot_accession(raw):
        return raw.upper()

    # Deterministic branch for exact-looking entry names such as NRP1_HUMAN.
    entry_name = raw.upper()
    if _looks_like_human_entry_name(entry_name):
        accession = _query_uniprot_primary_accession(
            f"id:{entry_name} AND reviewed:true"
        )
        if accession:
            return accession

    for query in _candidate_name_queries(raw):
        accession = _query_uniprot_primary_accession(query)
        if accession:
            return accession

    raise ValueError(
        f"Invalid UniProt accession ID: {identifier}. Use an accession "
        "(e.g. P15692) or a resolvable human gene/name (e.g. VEGFA or "
        "'VEGFA human')."
    )


def _resolve_identifiers(identifiers: list[str]) -> list[str]:
    """Resolve each provided identifier into a canonical accession."""
    return [_resolve_to_accession(identifier) for identifier in identifiers]


def _sanitize_alias_stem(stem: str) -> str:
    """Return a filesystem-safe stem for generated PDB aliases."""
    text = stem.strip().replace(" ", "_")
    return re.sub(r"[^A-Za-z0-9_.-]", "_", text)


def _alias_pdb_path(alias_stem: str, pdb_dir: str) -> str:
    """Return the user-facing PDB filename for an alias stem."""
    safe = _sanitize_alias_stem(alias_stem)
    return os.path.join(pdb_dir, f"{safe}.pdb")


def _extract_accession_and_entry_name(record_id: str) -> tuple[str, str | None]:
    """Parse UniProt FASTA record identifiers.

    Example:
        ``sp|P15692|VEGFA_HUMAN`` -> ("P15692", "VEGFA_HUMAN")
    """
    parts = record_id.split("|")
    if len(parts) >= 3:
        return parts[1], parts[2]
    return record_id, None


def _alias_stems_for_accession(
    accession_id: str,
    entry_name: str | None,
    pdb_name_mode: str,
) -> list[str]:
    """Return one or more output stems for a fetched PDB file."""
    if pdb_name_mode == "uniprot":
        preferred = entry_name or accession_id
        return [_sanitize_alias_stem(preferred)]
    if pdb_name_mode == "both":
        stems = [_sanitize_alias_stem(accession_id)]
        if entry_name:
            stems.append(_sanitize_alias_stem(entry_name))
        # Keep ordering stable while removing duplicates.
        deduped: list[str] = []
        for stem in stems:
            if stem not in deduped:
                deduped.append(stem)
        return deduped
    return [_sanitize_alias_stem(accession_id)]


def _materialize_pdb_aliases(
    raw_path: str,
    alias_stems: list[str],
    pdb_dir: str,
    *,
    overwrite: bool = False,
) -> tuple[list[str], list[str]]:
    """Copy a downloaded PDB artifact to one or more alias ``.pdb`` files.

    Returns ``(written_paths, skipped_paths)``. Both are empty when no
    readable raw artifact exists (for example when upstream retrieval reports
    an ID but does not materialize a file).
    """
    if not raw_path:
        return [], []

    raw_abs = os.path.abspath(raw_path)

    if not os.path.exists(raw_abs):
        return [], []

    alias_paths: list[str] = []
    skipped_paths: list[str] = []
    for stem in alias_stems:
        alias_path = _alias_pdb_path(stem, pdb_dir)
        alias_abs = os.path.abspath(alias_path)
        if os.path.exists(alias_abs) and not overwrite:
            skipped_paths.append(alias_path)
            continue
        if raw_abs != alias_abs:
            shutil.copyfile(raw_abs, alias_abs)
        alias_paths.append(alias_path)

    raw_name = os.path.basename(raw_abs).lower()
    if raw_name.startswith("pdb") and raw_name.endswith(".ent"):
        try:
            os.remove(raw_abs)
        except OSError:
            pass

    return alias_paths, skipped_paths


def remove_local_pdb_aliases(
    identifiers: list[str],
    *,
    pdb_dir: str,
) -> tuple[list[str], list[str]]:
    """Remove local fetched PDB aliases and return (removed, missing)."""
    removed: list[str] = []
    missing: list[str] = []

    for token in identifiers:
        raw = token.strip()
        if not raw:
            continue

        candidates: list[str] = []
        if os.path.isabs(raw) or os.path.sep in raw:
            candidates.append(raw)
        else:
            stem, ext = os.path.splitext(raw)
            if ext:
                candidates.append(os.path.join(pdb_dir, raw))
            else:
                safe_stem = _sanitize_alias_stem(stem or raw)
                candidates.append(_alias_pdb_path(safe_stem, pdb_dir))
                candidates.append(os.path.join(pdb_dir, f"{safe_stem}.ent"))

        target = next((c for c in candidates if os.path.exists(c)), None)
        if target is None:
            missing.append(raw)
            continue

        os.remove(target)
        removed.append(target)

    return removed, missing


def get_uniprot_data(
    accession_ids,
    fasta_file=None,
    csv_file=None,
    pdb_dir="pdb_files",
    pdb_name_mode="accession",
    force=False,
):
    """
    Fetch protein sequence and structure data from UniProt and PDB.

    Args:
        accession_ids (list[str]): List of UniProt accession IDs to process.
        fasta_file (str, optional): Path to save the combined FASTA sequences.
        csv_file (str, optional): Path to save structured sequence data as CSV.
        pdb_dir (str, optional): Directory to store downloaded PDB files.
        Defaults to "pdb_files".
        pdb_name_mode (str, optional): Naming mode for downloaded PDB files.
        Accepted values: "accession" (default), "uniprot", "both".
        force (bool, optional): If True, overwrite existing local alias
        files. By default existing aliases are preserved.

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
    resolved_accessions = _resolve_identifiers([str(token) for token in accession_ids])

    sequences_data = []
    structured_data = []

    # Ensure PDB directory exists
    os.makedirs(pdb_dir, exist_ok=True)

    for accession_id in resolved_accessions:
        # Fetch FASTA
        fasta_url = f"{_UNIPROT_BASE_URL}{accession_id}.fasta"
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
    accession_to_entry_name: dict[str, str | None] = {}
    for entry in records:
        accession, entry_name = _extract_accession_and_entry_name(entry.id)
        accession_to_entry_name[accession] = entry_name
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
    for accession_id in resolved_accessions:
        json_url = f"{_UNIPROT_BASE_URL}{accession_id}.json"
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
                alias_stems = _alias_stems_for_accession(
                    accession_id,
                    accession_to_entry_name.get(accession_id),
                    pdb_name_mode,
                )
                alias_paths, skipped_paths = _materialize_pdb_aliases(
                    raw_pdb_path,
                    alias_stems,
                    pdb_dir,
                    overwrite=force,
                )
                pdb_info[accession_id] = first_pdb_id
                if alias_paths:
                    print(
                        f"Downloaded PDB file for {accession_id}: "
                        f"{first_pdb_id} -> {', '.join(alias_paths)}"
                    )
                    if skipped_paths:
                        print(
                            "Kept existing local alias file(s) for "
                            f"{accession_id}: {', '.join(skipped_paths)}"
                        )
                elif skipped_paths:
                    print(
                        "Skipped download overwrite for "
                        f"{accession_id}; existing alias file(s): "
                        f"{', '.join(skipped_paths)}"
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
        nargs="*",
        help=(
            "One or more UniProt identifiers.  Accepted forms include "
            "accessions (e.g. P69905), FASTA-style IDs "
            "(sp|P15692|VEGFA_HUMAN), and common human gene/name inputs "
            "such as VEGFA or 'VEGFA human'.  "
            "Name-based inputs are resolved against reviewed human UniProt "
            "records (organism_id=9606).  "
            "The pipeline fetches sequence data and resolves PDB structures "
            "for each resolved accession."
        ),
    )
    parser.add_argument(
        "--search",
        default=None,
        help=(
            "Preview top UniProt reviewed-human matches for one ambiguous "
            "query term and exit (no downloads).  Example: --search "
            "'Neuropilin-1 human'."
        ),
    )
    parser.add_argument(
        "--search-limit",
        type=int,
        default=5,
        help=(
            "Number of rows to show for --search (default: 5, max: 20)."
        ),
    )
    parser.add_argument(
        "--remove",
        nargs="+",
        default=None,
        help=(
            "Delete local fetched PDB aliases from --pdb-dir and exit. "
            "Use accession/entry stems (e.g. P15692 or VEGFA_HUMAN) or "
            "paths."
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
            "Save structured sequence metadata to this CSV file.  Columns: "
            "ID, Name, Description, Sequence Length, Sequence.  Useful for "
            "record-keeping in large benchmarks."
        ),
    )
    parser.add_argument(
        "--pdb-name",
        choices=["accession", "uniprot", "both"],
        default="accession",
        help=(
            "Naming mode for downloaded PDB files (default: accession).  "
            "'accession' writes P15692.pdb, 'uniprot' writes VEGFA_HUMAN.pdb, "
            "'both' writes both aliases."
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
    parser.add_argument(
        "--force",
        action="store_true",
        help=(
            "Overwrite existing local alias files in --pdb-dir.  By "
            "default fetch keeps existing files and skips conflicting "
            "writes."
        ),
    )

    args = parser.parse_args(argv)

    if args.search_limit < 1 or args.search_limit > 20:
        print("ERROR: --search-limit must be between 1 and 20", file=sys.stderr)
        sys.exit(2)

    if args.search and args.remove:
        print("ERROR: --search and --remove cannot be used together", file=sys.stderr)
        sys.exit(2)

    if args.search and args.accessions:
        print(
            "ERROR: --search is preview-only; do not pass accessions with it",
            file=sys.stderr,
        )
        sys.exit(2)

    if args.remove and args.accessions:
        print(
            "ERROR: --remove only manages local files; do not pass accessions",
            file=sys.stderr,
        )
        sys.exit(2)

    if args.search:
        matches = search_uniprot_matches(args.search, limit=args.search_limit)
        if not matches:
            print(
                f"No reviewed human UniProt matches found for: {args.search}"
            )
            return

        print(
            f"Top {len(matches)} reviewed human UniProt match(es) for "
            f"'{args.search}':"
        )
        for row in matches:
            accession = row["accession"] or "-"
            entry_name = row["entry_name"] or "-"
            print(f"{accession} · {entry_name}")

            protein_name = row["protein_name"] or "(name unavailable)"
            gene_names = row["gene_names"] or "-"
            organism = row["organism"] or "-"
            length = row["length"] or "-"
            existence = row["protein_existence"] or "-"
            annotation = row["annotation_score"] or "-"
            print(
                f"{protein_name} · Gene: {gene_names} · {organism} · "
                f"{length} amino acids · {existence} · "
                f"Annotation score: {annotation}/5"
            )
            print()
        return

    if args.remove:
        os.makedirs(args.pdb_dir, exist_ok=True)
        removed, missing = remove_local_pdb_aliases(args.remove, pdb_dir=args.pdb_dir)
        if removed:
            print("Removed local PDB alias file(s):")
            for path in removed:
                print(f"  {path}")
        if missing:
            print("No local file found for:")
            for token in missing:
                print(f"  {token}")
        if not removed and not missing:
            print("No files matched --remove inputs.")
        return

    if not args.accessions:
        print(
            "ERROR: provide at least one accession or use --search/--remove",
            file=sys.stderr,
        )
        sys.exit(2)

    structured_data, pdb_info = get_uniprot_data(
        accession_ids=args.accessions,
        fasta_file=args.fasta,
        csv_file=args.csv,
        pdb_dir=args.pdb_dir,
        pdb_name_mode=args.pdb_name,
        force=args.force,
    )

    # Print a brief summary
    print(f"\nFetched {len(structured_data)} protein(s):")
    for entry in structured_data:
        pdb_id = pdb_info.get(entry["ID"].split("|")[1]
                              if "|" in entry["ID"]
                              else entry["ID"], "—")
        print(f"  {entry['ID']:30s}  {entry['Sequence Length']:>5d} aa  PDB: {pdb_id}")

    accession_to_entry_name: dict[str, str | None] = {}
    for entry in structured_data:
        accession, entry_name = _extract_accession_and_entry_name(entry["ID"])
        accession_to_entry_name[accession] = entry_name

    ready = [
        acc
        for acc, pdb_id in pdb_info.items()
        if pdb_id
        and any(
            os.path.exists(_alias_pdb_path(stem, args.pdb_dir))
            for stem in _alias_stems_for_accession(
                acc,
                accession_to_entry_name.get(acc),
                args.pdb_name,
            )
        )
    ]
    if ready:
        print("\nPair-ready PDB names:")
        for accession_id in ready:
            stems = _alias_stems_for_accession(
                accession_id,
                accession_to_entry_name.get(accession_id),
                args.pdb_name,
            )
            for stem in stems:
                print(f"  {stem} -> {_alias_pdb_path(stem, args.pdb_dir)}")
        print(
            "Use these filename stems directly in proteinA/proteinB, for "
            "example: P69905:P68871"
        )


if __name__ == "__main__":
    main()
