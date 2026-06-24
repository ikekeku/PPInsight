"""Find and optionally download native experimental complexes from RCSB.

This command is pair-aware: it resolves two user-provided identifiers to
UniProt accessions, searches the RCSB PDB assembly index for structures that
contain both partners, summarizes candidate biological assemblies, and can
download either a selected assembly or the top-ranked hit for each pair.
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
import sys
from pathlib import Path

import pandas as pd
import requests

from ppinsight.protein_fetch import _resolve_to_accession, _sanitize_alias_stem
from ppinsight.utils import _project_root, find_column

_RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
_RCSB_CORE_URL = "https://data.rcsb.org/rest/v1/core"
_RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download"
_REQUEST_TIMEOUT = 30

_METHOD_PRIORITY = {
    "X-RAY DIFFRACTION": 0,
    "ELECTRON MICROSCOPY": 1,
    "ELECTRON CRYSTALLOGRAPHY": 2,
    "NEUTRON DIFFRACTION": 3,
    "SOLUTION NMR": 4,
    "SOLID-STATE NMR": 5,
}

_COLUMN_ORDER = [
    "query_proteinA",
    "query_proteinB",
    "label",
    "family",
    "references",
    "accessionA",
    "accessionB",
    "status",
    "error",
    "rank",
    "assembly_identifier",
    "entry_id",
    "assembly_id",
    "experimental_method",
    "resolution_A",
    "release_date",
    "polymer_composition",
    "oligomeric_details",
    "oligomeric_count",
    "protein_entity_count",
    "protein_instance_count",
    "proteinA_entity_ids",
    "proteinA_asym_ids",
    "proteinA_auth_chains",
    "proteinB_entity_ids",
    "proteinB_asym_ids",
    "proteinB_auth_chains",
    "extra_uniprot_ids",
    "assembly_asym_ids",
    "assembly_oper_expression",
    "title",
    "download_url_pdb",
    "download_url_cif",
    "downloaded_path",
]


def _guess_sep(path: str | Path) -> str:
    """Infer CSV vs TSV separator from a filename."""
    return "\t" if str(path).lower().endswith(".tsv") else ","


def _dedupe(values: list[str]) -> list[str]:
    """Return values with stable de-duplication and blank removal."""
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        text = str(value).strip()
        if not text or text in seen:
            continue
        seen.add(text)
        ordered.append(text)
    return ordered


def _stringify(value: object) -> str:
    """Normalize optional values into plain strings."""
    if value is None:
        return ""
    if pd.isna(value):
        return ""
    return str(value).strip()


def _semicolon_join(values: list[str]) -> str:
    """Serialize a list field for tabular output."""
    return ";".join(_dedupe(values))


def _resolution_value(entry_data: dict) -> float | None:
    """Extract the first combined resolution value when present."""
    resolutions = (
        entry_data.get("rcsb_entry_info", {}) or {}
    ).get("resolution_combined") or []
    for item in resolutions:
        try:
            return float(item)
        except (TypeError, ValueError):
            continue
    return None


def _experimental_method(entry_data: dict) -> str:
    """Return a stable semicolon-joined experimental method string."""
    methods: list[str] = []
    for item in entry_data.get("exptl") or []:
        method = _stringify((item or {}).get("method"))
        if method:
            methods.append(method)
    return _semicolon_join(methods)


def _assembly_asym_ids(assembly_data: dict) -> tuple[list[str], list[str]]:
    """Return assembly asym IDs plus their generation operators."""
    asym_ids: list[str] = []
    oper_expressions: list[str] = []
    for item in assembly_data.get("pdbx_struct_assembly_gen") or []:
        item = item or {}
        asym_ids.extend(item.get("asym_id_list") or [])
        operator = _stringify(item.get("oper_expression"))
        if operator:
            oper_expressions.append(operator)
    return _dedupe(asym_ids), _dedupe(oper_expressions)


def _extract_uniprot_ids(entity_data: dict) -> list[str]:
    """Return UniProt accessions attached to a polymer entity."""
    container = entity_data.get("rcsb_polymer_entity_container_identifiers", {}) or {}
    accessions = [str(item).upper() for item in container.get("uniprot_ids") or []]
    for item in container.get("reference_sequence_identifiers") or []:
        item = item or {}
        if _stringify(item.get("database_name")).lower() != "uniprot":
            continue
        accession = _stringify(item.get("database_accession")).upper()
        if accession:
            accessions.append(accession)
    return _dedupe(accessions)


def _method_rank(methods: str) -> int:
    """Return a heuristic priority for experimental methods."""
    if not methods:
        return 99
    return min(
        (
            _METHOD_PRIORITY.get(method.strip().upper(), 99)
            for method in methods.split(";")
        ),
        default=99,
    )


def _int_or_large(value: object) -> int:
    """Parse an integer-like field, falling back to a large sentinel."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return 999_999


def _float_or_large(value: object) -> float:
    """Parse a float-like field, falling back to a large sentinel."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.inf


def _candidate_rank_key(record: dict[str, object]) -> tuple[object, ...]:
    """Heuristic ordering for candidate assemblies."""
    extra_ids = _stringify(record.get("extra_uniprot_ids"))
    extra_count = len([item for item in extra_ids.split(";") if item])
    return (
        extra_count,
        _int_or_large(record.get("protein_entity_count")),
        _int_or_large(record.get("protein_instance_count")),
        _method_rank(_stringify(record.get("experimental_method"))),
        _float_or_large(record.get("resolution_A")),
        _stringify(record.get("entry_id")),
        _stringify(record.get("assembly_id")),
    )


def _search_payload(accession_a: str, accession_b: str, max_hits: int) -> dict:
    """Build an RCSB assembly-level search query for one pair."""
    accessions = _dedupe([accession_a.upper(), accession_b.upper()])
    nodes = [
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": (
                    "rcsb_polymer_entity_container_identifiers."
                    "reference_sequence_identifiers.database_name"
                ),
                "operator": "exact_match",
                "value": "UniProt",
            },
        }
    ]
    for accession in accessions:
        nodes.append(
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": (
                        "rcsb_polymer_entity_container_identifiers."
                        "reference_sequence_identifiers.database_accession"
                    ),
                    "operator": "exact_match",
                    "value": accession,
                },
            }
        )

    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": nodes,
        },
        "return_type": "assembly",
        "request_options": {
            "results_content_type": ["experimental"],
            "paginate": {"start": 0, "rows": max_hits},
        },
    }


def search_native_assemblies(
    accession_a: str,
    accession_b: str,
    *,
    max_hits: int = 10,
) -> list[tuple[str, str]]:
    """Return candidate ``(entry_id, assembly_id)`` tuples from RCSB."""
    response = requests.post(
        _RCSB_SEARCH_URL,
        json=_search_payload(accession_a, accession_b, max_hits),
        timeout=_REQUEST_TIMEOUT,
    )
    response.raise_for_status()
    if getattr(response, "status_code", 200) == 204:
        return []

    try:
        payload = response.json()
    except ValueError as exc:
        raise RuntimeError("RCSB search returned invalid JSON") from exc

    pairs: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for item in payload.get("result_set") or []:
        identifier = _stringify((item or {}).get("identifier"))
        if "-" not in identifier:
            continue
        entry_id, assembly_id = identifier.split("-", 1)
        pair = (entry_id.upper(), assembly_id)
        if pair in seen:
            continue
        seen.add(pair)
        pairs.append(pair)
    return pairs


def _get_json(url: str) -> dict:
    """Fetch JSON from a REST endpoint with consistent error handling."""
    response = requests.get(url, timeout=_REQUEST_TIMEOUT)
    response.raise_for_status()
    payload = response.json()
    if not isinstance(payload, dict):
        raise RuntimeError(f"Unexpected payload type from {url}")
    return payload


def _download_url(entry_id: str, assembly_id: str, file_format: str) -> str:
    """Return the RCSB download URL for a biological assembly."""
    if file_format == "cif":
        return f"{_RCSB_DOWNLOAD_URL}/{entry_id}-assembly{assembly_id}.cif.gz"
    return f"{_RCSB_DOWNLOAD_URL}/{entry_id}.pdb{assembly_id}.gz"


def _load_pairs_file(path: str | Path) -> list[dict[str, str]]:
    """Load a flat pairs CSV/TSV file into query records."""
    df = pd.read_csv(path, sep=_guess_sep(path))
    protein_a_col = find_column(df.columns, "proteinA", "protein_a", "proteina")
    protein_b_col = find_column(df.columns, "proteinB", "protein_b", "proteinb")
    if protein_a_col is None or protein_b_col is None:
        raise ValueError(
            "pairs file must contain proteinA and proteinB columns. "
            f"Got: {list(df.columns)}"
        )

    label_col = find_column(df.columns, "label")
    family_col = find_column(df.columns, "family")
    refs_col = find_column(df.columns, "references", "reference")

    records: list[dict[str, str]] = []
    for _, row in df.iterrows():
        protein_a = _stringify(row.get(protein_a_col))
        protein_b = _stringify(row.get(protein_b_col))
        if not protein_a or not protein_b:
            continue
        records.append(
            {
                "query_proteinA": protein_a,
                "query_proteinB": protein_b,
                "label": _stringify(row.get(label_col)) if label_col else "",
                "family": _stringify(row.get(family_col)) if family_col else "",
                "references": _stringify(row.get(refs_col)) if refs_col else "",
            }
        )
    return records


def _status_row(
    query_record: dict[str, str],
    *,
    accession_a: str = "",
    accession_b: str = "",
    status: str,
    error: str = "",
    entry_id: str = "",
    assembly_id: str = "",
) -> dict[str, object]:
    """Return a standardized status row for no-hit and error cases."""
    row = {
        "query_proteinA": query_record["query_proteinA"],
        "query_proteinB": query_record["query_proteinB"],
        "label": query_record.get("label", ""),
        "family": query_record.get("family", ""),
        "references": query_record.get("references", ""),
        "accessionA": accession_a,
        "accessionB": accession_b,
        "status": status,
        "error": error,
        "rank": "",
        "assembly_identifier": (
            f"{entry_id}-{assembly_id}" if entry_id and assembly_id else ""
        ),
        "entry_id": entry_id,
        "assembly_id": assembly_id,
        "download_url_pdb": (
            _download_url(entry_id, assembly_id, "pdb")
            if entry_id and assembly_id
            else ""
        ),
        "download_url_cif": (
            _download_url(entry_id, assembly_id, "cif")
            if entry_id and assembly_id
            else ""
        ),
    }
    return row


def _assembly_record(
    query_record: dict[str, str],
    accession_a: str,
    accession_b: str,
    entry_id: str,
    assembly_id: str,
    *,
    entry_cache: dict[str, dict],
    assembly_cache: dict[tuple[str, str], dict],
    entity_cache: dict[tuple[str, str], dict],
    instance_cache: dict[tuple[str, str], dict],
) -> dict[str, object]:
    """Build one candidate row with metadata for a single assembly."""
    if entry_id not in entry_cache:
        entry_cache[entry_id] = _get_json(f"{_RCSB_CORE_URL}/entry/{entry_id}")
    entry_data = entry_cache[entry_id]

    assembly_key = (entry_id, assembly_id)
    if assembly_key not in assembly_cache:
        assembly_cache[assembly_key] = _get_json(
            f"{_RCSB_CORE_URL}/assembly/{entry_id}/{assembly_id}"
        )
    assembly_data = assembly_cache[assembly_key]

    assembly_asym_ids, oper_expressions = _assembly_asym_ids(assembly_data)
    instance_rows: list[dict[str, object]] = []
    assembly_accessions: set[str] = set()

    for asym_id in assembly_asym_ids:
        instance_key = (entry_id, asym_id)
        if instance_key not in instance_cache:
            instance_cache[instance_key] = _get_json(
                f"{_RCSB_CORE_URL}/polymer_entity_instance/{entry_id}/{asym_id}"
            )
        instance_data = instance_cache[instance_key]
        instance_ids = (
            instance_data.get("rcsb_polymer_entity_instance_container_identifiers", {})
            or {}
        )
        entity_id = _stringify(instance_ids.get("entity_id"))
        auth_asym_id = _stringify(instance_ids.get("auth_asym_id"))
        if not entity_id:
            continue

        entity_key = (entry_id, entity_id)
        if entity_key not in entity_cache:
            entity_cache[entity_key] = _get_json(
                f"{_RCSB_CORE_URL}/polymer_entity/{entry_id}/{entity_id}"
            )
        entity_data = entity_cache[entity_key]
        accessions = _extract_uniprot_ids(entity_data)
        assembly_accessions.update(accessions)

        instance_rows.append(
            {
                "entity_id": entity_id,
                "asym_id": asym_id,
                "auth_asym_id": auth_asym_id,
                "uniprot_ids": accessions,
            }
        )

    def mapping_for(accession: str) -> tuple[str, str, str]:
        matched = [row for row in instance_rows if accession in row["uniprot_ids"]]
        entity_ids = _semicolon_join([str(row["entity_id"]) for row in matched])
        asym_ids = _semicolon_join([str(row["asym_id"]) for row in matched])
        auth_ids = _semicolon_join([str(row["auth_asym_id"]) for row in matched])
        return entity_ids, asym_ids, auth_ids

    a_entities, a_asym_ids, a_auth_ids = mapping_for(accession_a)
    b_entities, b_asym_ids, b_auth_ids = mapping_for(accession_b)

    requested = {accession_a.upper(), accession_b.upper()}
    extra_accessions = sorted(
        accession for accession in assembly_accessions if accession not in requested
    )

    assembly_info = assembly_data.get("rcsb_assembly_info", {}) or {}
    struct_assembly = assembly_data.get("pdbx_struct_assembly", {}) or {}
    resolution = _resolution_value(entry_data)

    return {
        "query_proteinA": query_record["query_proteinA"],
        "query_proteinB": query_record["query_proteinB"],
        "label": query_record.get("label", ""),
        "family": query_record.get("family", ""),
        "references": query_record.get("references", ""),
        "accessionA": accession_a,
        "accessionB": accession_b,
        "status": "ok",
        "error": "",
        "rank": "",
        "assembly_identifier": f"{entry_id}-{assembly_id}",
        "entry_id": entry_id,
        "assembly_id": assembly_id,
        "experimental_method": _experimental_method(entry_data),
        "resolution_A": resolution if resolution is not None else "",
        "release_date": _stringify(
            (entry_data.get("rcsb_accession_info", {}) or {}).get(
                "initial_release_date"
            )
        ),
        "polymer_composition": _stringify(assembly_info.get("polymer_composition")),
        "oligomeric_details": _stringify(struct_assembly.get("oligomeric_details")),
        "oligomeric_count": _stringify(struct_assembly.get("oligomeric_count")),
        "protein_entity_count": _stringify(
            assembly_info.get("polymer_entity_count_protein")
        ),
        "protein_instance_count": _stringify(
            assembly_info.get("polymer_entity_instance_count_protein")
            or assembly_info.get("polymer_entity_instance_count")
        ),
        "proteinA_entity_ids": a_entities,
        "proteinA_asym_ids": a_asym_ids,
        "proteinA_auth_chains": a_auth_ids,
        "proteinB_entity_ids": b_entities,
        "proteinB_asym_ids": b_asym_ids,
        "proteinB_auth_chains": b_auth_ids,
        "extra_uniprot_ids": _semicolon_join(extra_accessions),
        "assembly_asym_ids": _semicolon_join(assembly_asym_ids),
        "assembly_oper_expression": _semicolon_join(oper_expressions),
        "title": _stringify((entry_data.get("struct", {}) or {}).get("title")),
        "download_url_pdb": _download_url(entry_id, assembly_id, "pdb"),
        "download_url_cif": _download_url(entry_id, assembly_id, "cif"),
    }


def _reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Apply a stable, user-facing column order."""
    ordered = [column for column in _COLUMN_ORDER if column in df.columns]
    remaining = [column for column in df.columns if column not in ordered]
    return df[ordered + remaining]


def find_native_complexes(
    query_records: list[dict[str, str]],
    *,
    max_hits: int = 10,
) -> pd.DataFrame:
    """Resolve pairs, search RCSB assemblies, and return a candidate table."""
    rows: list[dict[str, object]] = []
    entry_cache: dict[str, dict] = {}
    assembly_cache: dict[tuple[str, str], dict] = {}
    entity_cache: dict[tuple[str, str], dict] = {}
    instance_cache: dict[tuple[str, str], dict] = {}

    for query_record in query_records:
        accession_a = ""
        accession_b = ""

        try:
            accession_a = _resolve_to_accession(query_record["query_proteinA"])
            accession_b = _resolve_to_accession(query_record["query_proteinB"])
        except ValueError as exc:
            rows.append(
                _status_row(
                    query_record,
                    status="resolve_error",
                    error=str(exc),
                )
            )
            continue

        try:
            assembly_pairs = search_native_assemblies(
                accession_a,
                accession_b,
                max_hits=max_hits,
            )
        except (requests.RequestException, RuntimeError, ValueError) as exc:
            rows.append(
                _status_row(
                    query_record,
                    accession_a=accession_a,
                    accession_b=accession_b,
                    status="search_error",
                    error=str(exc),
                )
            )
            continue

        if not assembly_pairs:
            rows.append(
                _status_row(
                    query_record,
                    accession_a=accession_a,
                    accession_b=accession_b,
                    status="no_hit",
                    error="No experimental assembly matches found.",
                )
            )
            continue

        ok_rows: list[dict[str, object]] = []
        error_rows: list[dict[str, object]] = []
        for entry_id, assembly_id in assembly_pairs:
            try:
                ok_rows.append(
                    _assembly_record(
                        query_record,
                        accession_a,
                        accession_b,
                        entry_id,
                        assembly_id,
                        entry_cache=entry_cache,
                        assembly_cache=assembly_cache,
                        entity_cache=entity_cache,
                        instance_cache=instance_cache,
                    )
                )
            except (requests.RequestException, RuntimeError, ValueError) as exc:
                error_rows.append(
                    _status_row(
                        query_record,
                        accession_a=accession_a,
                        accession_b=accession_b,
                        status="metadata_error",
                        error=f"{entry_id}-{assembly_id}: {exc}",
                        entry_id=entry_id,
                        assembly_id=assembly_id,
                    )
                )

        ok_rows.sort(key=_candidate_rank_key)
        for rank, row in enumerate(ok_rows, start=1):
            row["rank"] = rank
        rows.extend(ok_rows)
        rows.extend(error_rows)

    if not rows:
        return pd.DataFrame(columns=_COLUMN_ORDER)
    return _reorder_columns(pd.DataFrame(rows))


def download_assembly(
    entry_id: str,
    assembly_id: str,
    *,
    accession_a: str,
    accession_b: str,
    download_dir: str | Path,
    file_format: str = "pdb",
) -> str:
    """Download and decompress one biological assembly into ``download_dir``."""
    if file_format not in {"pdb", "cif"}:
        raise ValueError("file_format must be 'pdb' or 'cif'")

    pair_dir = Path(download_dir) / (
        f"{_sanitize_alias_stem(accession_a)}_vs_{_sanitize_alias_stem(accession_b)}"
    )
    pair_dir.mkdir(parents=True, exist_ok=True)

    suffix = ".pdb" if file_format == "pdb" else ".cif"
    output_path = pair_dir / f"{entry_id}_assembly{assembly_id}{suffix}"
    if output_path.exists():
        return str(output_path)

    response = requests.get(
        _download_url(entry_id, assembly_id, file_format),
        timeout=_REQUEST_TIMEOUT,
    )
    response.raise_for_status()
    output_path.write_bytes(gzip.decompress(response.content))
    return str(output_path)


def download_selected_candidates(
    candidates_df: pd.DataFrame,
    *,
    download_dir: str | Path,
    file_format: str = "pdb",
    selected_assembly: str,
) -> pd.DataFrame:
    """Download one explicitly selected assembly from the candidate table."""
    if not selected_assembly:
        raise ValueError("selected_assembly is required")

    if candidates_df.empty:
        return pd.DataFrame(columns=list(candidates_df.columns) + ["downloaded_path"])

    ok_rows = candidates_df[candidates_df["status"] == "ok"].copy()
    if ok_rows.empty:
        return pd.DataFrame(columns=list(candidates_df.columns) + ["downloaded_path"])

    ok_rows = ok_rows[ok_rows["assembly_identifier"] == selected_assembly]
    if ok_rows.empty:
        raise ValueError(
            f"Assembly not found in candidate table: {selected_assembly}"
        )
    ok_rows = ok_rows.head(1)

    downloads: list[dict[str, object]] = []
    for _, row in ok_rows.iterrows():
        downloaded_path = download_assembly(
            _stringify(row.get("entry_id")),
            _stringify(row.get("assembly_id")),
            accession_a=_stringify(row.get("accessionA")),
            accession_b=_stringify(row.get("accessionB")),
            download_dir=download_dir,
            file_format=file_format,
        )
        record = row.to_dict()
        record["downloaded_path"] = downloaded_path
        downloads.append(record)

    if not downloads:
        return pd.DataFrame(columns=list(candidates_df.columns) + ["downloaded_path"])
    return _reorder_columns(pd.DataFrame(downloads))


def _build_parser() -> argparse.ArgumentParser:
    """Construct the CLI parser for ``ppinsight fetch-native``."""
    parser = argparse.ArgumentParser(
        prog="ppinsight fetch-native",
        description=(
            "Find experimental native complexes for one protein pair or a "
            "pairs file by searching RCSB biological assemblies."
        ),
    )
    parser.add_argument(
        "identifiers",
        nargs="*",
        help=(
            "Two protein identifiers in single-pair mode. Accepts UniProt "
            "accessions, FASTA-style IDs, and the same human gene/name forms "
            "as 'ppinsight fetch'."
        ),
    )
    parser.add_argument(
        "--pairs",
        default=None,
        help=(
            "Path to a flat pairs CSV/TSV with proteinA and proteinB columns. "
            "Optional label/family/references columns are carried into the "
            "candidate table."
        ),
    )
    parser.add_argument(
        "--max-hits",
        type=int,
        default=10,
        metavar="N",
        help="Maximum number of assembly hits to keep per pair (default: 10).",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help=(
            "Write the candidate table to CSV/TSV (extension decides format). "
            "If omitted, the table is written to stdout as TSV."
        ),
    )
    parser.add_argument(
        "--download-dir",
        default=os.path.join(_project_root(), "data", "input", "native_complexes"),
        help=(
            "Directory for downloaded native assemblies when --select is used."
        ),
    )
    parser.add_argument(
        "--format",
        choices=["pdb", "cif"],
        default="pdb",
        help=(
            "Structure format to download (default: pdb). Use cif when the "
            "assembly is too large or you want the mmCIF source file."
        ),
    )
    parser.add_argument(
        "--select",
        default=None,
        help=(
            "Download one specific assembly identifier from the candidate "
            "table, for example '3V2A-1'. Single-pair mode only."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    """CLI entry point for ``ppinsight fetch-native``."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.pairs and args.identifiers:
        parser.error("Use either --pairs or two positional identifiers, not both")
    if not args.pairs and len(args.identifiers) != 2:
        parser.error("Provide either --pairs FILE or exactly two identifiers")
    if args.select and args.pairs:
        parser.error("--select is only supported in single-pair mode")
    if args.max_hits < 1:
        parser.error("--max-hits must be at least 1")

    try:
        if args.pairs:
            query_records = _load_pairs_file(args.pairs)
        else:
            query_records = [
                {
                    "query_proteinA": args.identifiers[0],
                    "query_proteinB": args.identifiers[1],
                    "label": "",
                    "family": "",
                    "references": "",
                }
            ]
    except ValueError as exc:
        parser.error(str(exc))

    candidates_df = find_native_complexes(query_records, max_hits=args.max_hits)

    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        candidates_df.to_csv(output_path, sep=_guess_sep(output_path), index=False)
        print(
            f"Wrote {len(candidates_df)} candidate rows to {output_path}",
            file=sys.stderr,
        )
    else:
        candidates_df.to_csv(sys.stdout, sep="\t", index=False)

    if args.select:
        try:
            downloads_df = download_selected_candidates(
                candidates_df,
                download_dir=args.download_dir,
                file_format=args.format,
                selected_assembly=args.select,
            )
        except ValueError as exc:
            parser.error(str(exc))
        if downloads_df.empty:
            print("No candidate assemblies were downloaded.", file=sys.stderr)
        else:
            for _, row in downloads_df.iterrows():
                print(
                    (
                        f"Downloaded {row['assembly_identifier']} -> "
                        f"{row['downloaded_path']}"
                    ),
                    file=sys.stderr,
                )


if __name__ == "__main__":
    main()
