"""Tests for the native-complex search and download workflow."""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

from ppinsight import fetch_native


class DummyResponse:
    """Minimal requests.Response stand-in for unit tests."""

    def __init__(self, *, json_data=None, content=b"", status_code=200):
        self._json_data = json_data
        self.content = content
        self.status_code = status_code

    def raise_for_status(self):
        return None

    def json(self):
        return self._json_data


def test_find_native_complexes_builds_ranked_candidate_rows(monkeypatch):
    """Assembly search should return ranked rows with chain mappings."""

    accession_map = {"VEGFR2": "P35968", "VEGFA": "P15692"}

    def fake_resolve(identifier):
        return accession_map.get(identifier, identifier)

    def fake_post(url, json, timeout):
        assert url == fetch_native._RCSB_SEARCH_URL
        assert json["return_type"] == "assembly"
        return DummyResponse(
            json_data={
                "result_set": [
                    {"identifier": "4ABC-1"},
                    {"identifier": "3V2A-1"},
                ]
            }
        )

    entry_payloads = {
        "3V2A": {
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"resolution_combined": [2.4]},
            "rcsb_accession_info": {"initial_release_date": "2012-01-04"},
            "struct": {"title": "VEGFR2 VEGFA complex"},
        },
        "4ABC": {
            "exptl": [{"method": "ELECTRON MICROSCOPY"}],
            "rcsb_entry_info": {"resolution_combined": [4.1]},
            "rcsb_accession_info": {"initial_release_date": "2018-06-01"},
            "struct": {"title": "Broader signaling complex"},
        },
    }

    assembly_payloads = {
        ("3V2A", "1"): {
            "pdbx_struct_assembly": {
                "oligomeric_count": 2,
                "oligomeric_details": "heterodimeric",
            },
            "pdbx_struct_assembly_gen": [
                {"asym_id_list": ["A", "B"], "oper_expression": "1"}
            ],
            "rcsb_assembly_info": {
                "polymer_composition": "heteromeric protein",
                "polymer_entity_count_protein": 2,
                "polymer_entity_instance_count_protein": 2,
            },
        },
        ("4ABC", "1"): {
            "pdbx_struct_assembly": {
                "oligomeric_count": 3,
                "oligomeric_details": "heterotrimeric",
            },
            "pdbx_struct_assembly_gen": [
                {"asym_id_list": ["A", "B", "C"], "oper_expression": "1"}
            ],
            "rcsb_assembly_info": {
                "polymer_composition": "heteromeric protein",
                "polymer_entity_count_protein": 3,
                "polymer_entity_instance_count_protein": 3,
            },
        },
    }

    instance_payloads = {
        ("3V2A", "A"): {
            "rcsb_polymer_entity_instance_container_identifiers": {
                "entity_id": "1",
                "auth_asym_id": "R",
            }
        },
        ("3V2A", "B"): {
            "rcsb_polymer_entity_instance_container_identifiers": {
                "entity_id": "2",
                "auth_asym_id": "A",
            }
        },
        ("4ABC", "A"): {
            "rcsb_polymer_entity_instance_container_identifiers": {
                "entity_id": "1",
                "auth_asym_id": "R",
            }
        },
        ("4ABC", "B"): {
            "rcsb_polymer_entity_instance_container_identifiers": {
                "entity_id": "2",
                "auth_asym_id": "A",
            }
        },
        ("4ABC", "C"): {
            "rcsb_polymer_entity_instance_container_identifiers": {
                "entity_id": "3",
                "auth_asym_id": "X",
            }
        },
    }

    entity_payloads = {
        ("3V2A", "1"): {
            "rcsb_polymer_entity_container_identifiers": {
                "uniprot_ids": ["P35968"]
            }
        },
        ("3V2A", "2"): {
            "rcsb_polymer_entity_container_identifiers": {
                "uniprot_ids": ["P15692"]
            }
        },
        ("4ABC", "1"): {
            "rcsb_polymer_entity_container_identifiers": {
                "uniprot_ids": ["P35968"]
            }
        },
        ("4ABC", "2"): {
            "rcsb_polymer_entity_container_identifiers": {
                "uniprot_ids": ["P15692"]
            }
        },
        ("4ABC", "3"): {
            "rcsb_polymer_entity_container_identifiers": {
                "uniprot_ids": ["Q99999"]
            }
        },
    }

    def fake_get(url, timeout):
        if "/entry/" in url:
            entry_id = url.rsplit("/", 1)[-1]
            return DummyResponse(json_data=entry_payloads[entry_id])
        if "/assembly/" in url:
            entry_id, assembly_id = url.rsplit("/", 2)[-2:]
            return DummyResponse(json_data=assembly_payloads[(entry_id, assembly_id)])
        if "/polymer_entity_instance/" in url:
            entry_id, asym_id = url.rsplit("/", 2)[-2:]
            return DummyResponse(json_data=instance_payloads[(entry_id, asym_id)])
        if "/polymer_entity/" in url:
            entry_id, entity_id = url.rsplit("/", 2)[-2:]
            return DummyResponse(json_data=entity_payloads[(entry_id, entity_id)])
        raise AssertionError(f"Unexpected URL: {url}")

    monkeypatch.setattr(fetch_native, "_resolve_to_accession", fake_resolve)
    monkeypatch.setattr(fetch_native.requests, "post", fake_post)
    monkeypatch.setattr(fetch_native.requests, "get", fake_get)

    candidates = fetch_native.find_native_complexes(
        [
            {
                "query_proteinA": "VEGFR2",
                "query_proteinB": "VEGFA",
                "label": "interaction",
                "family": "RTK",
                "references": "6,7",
            }
        ],
        max_hits=5,
    )

    assert list(candidates["assembly_identifier"])[:2] == ["3V2A-1", "4ABC-1"]
    assert list(candidates["rank"])[:2] == [1, 2]
    assert candidates.loc[0, "proteinA_auth_chains"] == "R"
    assert candidates.loc[0, "proteinB_auth_chains"] == "A"
    assert candidates.loc[0, "extra_uniprot_ids"] == ""
    assert candidates.loc[1, "extra_uniprot_ids"] == "Q99999"
    assert candidates.loc[0, "label"] == "interaction"
    assert candidates.loc[0, "family"] == "RTK"


def test_find_native_complexes_records_no_hit_rows(monkeypatch):
    """Pairs with no assembly hits should stay visible in the output table."""

    monkeypatch.setattr(fetch_native, "_resolve_to_accession", lambda token: token)
    monkeypatch.setattr(
        fetch_native.requests,
        "post",
        lambda url, json, timeout: DummyResponse(json_data={"result_set": []}),
    )

    candidates = fetch_native.find_native_complexes(
        [{"query_proteinA": "P11111", "query_proteinB": "P22222"}],
        max_hits=3,
    )

    assert len(candidates) == 1
    assert candidates.loc[0, "status"] == "no_hit"
    assert candidates.loc[0, "accessionA"] == "P11111"
    assert candidates.loc[0, "accessionB"] == "P22222"


def test_download_selected_candidates_writes_decompressed_structure(
    tmp_path,
    monkeypatch,
):
    """Selected downloads should be materialized as plain structure files."""

    payload = gzip.compress(b"HEADER TEST\n")

    def fake_get(url, timeout):
        assert url.endswith("3V2A.pdb1.gz")
        return DummyResponse(content=payload)

    monkeypatch.setattr(fetch_native.requests, "get", fake_get)

    candidates = pd.DataFrame(
        [
            {
                "query_proteinA": "VEGFR2",
                "query_proteinB": "VEGFA",
                "accessionA": "P35968",
                "accessionB": "P15692",
                "status": "ok",
                "rank": 1,
                "assembly_identifier": "3V2A-1",
                "entry_id": "3V2A",
                "assembly_id": "1",
            }
        ]
    )

    downloaded = fetch_native.download_selected_candidates(
        candidates,
        download_dir=tmp_path,
        file_format="pdb",
        selected_assembly="3V2A-1",
    )

    output_path = Path(downloaded.loc[0, "downloaded_path"])
    assert output_path.exists()
    assert output_path.read_text(encoding="utf-8") == "HEADER TEST\n"
    assert output_path.name == "3V2A_assembly1.pdb"
