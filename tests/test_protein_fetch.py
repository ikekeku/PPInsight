"""
Tests for the protein_fetch module.

Consolidated from test_protein_fetch.py (Walter) and
test_protein_fetch_1.py (Rita).
"""

from pathlib import Path

import pytest

import ppinsight
from ppinsight import protein_fetch

# ── Smoke tests ──────────────────────────────────────────────────

def test_smoke_walter(tmp_path, monkeypatch):
    """
    author: wmavila
    reviewer: mkamenetskiy
    category: smoke test

    Ensure that get_uniprot_data() runs without raising exceptions.
    """
    monkeypatch.chdir(tmp_path)
    accession_ids = ["P69905", "P68871"]
    ppinsight.protein_fetch.get_uniprot_data(
        accession_ids,
        fasta_file="proteins.fasta",
        csv_file="protein",
    )


def test_smoke_rita(tmp_path, monkeypatch):
    """
    author: mkamenetskiy
    reviewer: wmavila
    category: smoke test

    Simple smoke test with a single accession.
    """
    monkeypatch.chdir(tmp_path)
    protein_fetch.get_uniprot_data(
        ["P15692"],
        fasta_file="proteins.fasta",
        csv_file="protein",
    )


# ── One-shot tests ───────────────────────────────────────────────

def test_oneshot_walter(tmp_path, monkeypatch):
    """
    author: wmavila
    reviewer: mkamenetskiy
    category: one-shot test

    Test that get_uniprot_data() returns the expected protein data for P69905
    (human hemoglobin alpha chain — 142 amino acids, at least one PDB).
    """
    monkeypatch.chdir(tmp_path)
    structured_data, pdb_info = ppinsight.protein_fetch.get_uniprot_data(
        ["P69905"],
    )
    expected_length = 142

    assert structured_data[0]["Sequence Length"] == expected_length, (
        f"Expected {expected_length}, got {structured_data[0]['Sequence Length']}"
    )
    assert pdb_info["P69905"] is not None, "Expected a PDB ID for P69905"


def test_get_uniprot_data_p15692_structured(tmp_path, monkeypatch):
    """
    author: mkamenetskiy
    reviewer: wmavila
    category: one-shot test

    One-shot test with a known UniProt accession P15692 (VEGFA_HUMAN).
    """
    monkeypatch.chdir(tmp_path)
    expected_structured_data = [{
        'ID': 'sp|P15692|VEGFA_HUMAN',
        'Name': 'sp|P15692|VEGFA_HUMAN',
        'Description': 'sp|P15692|VEGFA_HUMAN Vascular endothelial growth '
        'factor A, long form OS=Homo sapiens OX=9606 GN=VEGFA PE=1 SV=3',
        'Sequence Length': 395,
        'Sequence': (
            'MTDRQTDTAPSPSYHLLPGRRRTVDAAASRGQGPEPAPGGGVEGVGARGVALKLFVQLLG'
            'CSRFGGAVVRAGEAEPSGAARSASSGREEPQPEEGEEEEEKEEERGPQWRLGARKPGSWT'
            'GEAAVCADSAPAARAPQALARASGRGGRVARRGAEESGPPHSPSRRGSASRAGPGRASET'
            'MNFLLSWVHWSLALLLYLHHAKWSQAAPMAEGGGQNHHEVVKFMDVYQRSYCHPIETLV'
            'DIFQEYPDEIEYIFKPSCVPLMRCGGCCNDEGLECVPTEESNITMQIMRIKPHQGQHIG'
            'EMSFLQHNKCECRPKKDRARQEKKSVRGKGKGQKRKRKKSRYKSWSVPCGPCSERRKHL'
            'FVQDPQTCKCSCKNTDSRCKARQLELNERTCRCDKPRR'
        ),
    }]
    actual_structured_data, _ = protein_fetch.get_uniprot_data(["P15692"])
    assert actual_structured_data == expected_structured_data


def test_get_uniprot_data_p15692_pdbinfo(tmp_path, monkeypatch):
    """
    author: mkamenetskiy
    reviewer: wmavila
    category: one-shot test

    One-shot test using known UniProt accession P15692 (VEGFA_HUMAN).
    """
    monkeypatch.chdir(tmp_path)
    expected_pdb_info = {'P15692': '1BJ1'}
    _, actual_pdb_info = protein_fetch.get_uniprot_data(["P15692"])
    assert actual_pdb_info == expected_pdb_info


# ── Edge-case tests ──────────────────────────────────────────────

def test_edgecase_walter(tmp_path, monkeypatch):
    """
    author: wmavila
    reviewer: mkamenetskiy
    category: edge-case test

    Fake UniProt accession → should raise ValueError.
    """
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match="Invalid UniProt accession ID"):
        ppinsight.protein_fetch.get_uniprot_data(["fake_id"])


def test_nonexistent_protein(tmp_path, monkeypatch):
    """
    author: mkamenetskiy
    reviewer: wmavila
    category: edge test

    Non-existent protein should raise ValueError even if other IDs are valid.
    """
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match="Invalid UniProt accession ID"):
        protein_fetch.get_uniprot_data(["P000000000"])


# ── Pattern tests ────────────────────────────────────────────────

def test_pattern_walter(tmp_path, monkeypatch):
    """
    author: wmavila
    reviewer: mkamenetskiy
    category: pattern test

    Validate known patterns:
    1. The sequence must not be empty.
    2. The ID should appear in the Description.
    """
    monkeypatch.chdir(tmp_path)
    accession_ids = ["P69905", "P68871"]
    sequence_data, pdb_info = ppinsight.protein_fetch.get_uniprot_data(
        accession_ids,
    )
    for entry in sequence_data:
        assert entry["Sequence"], f"Empty sequence for {entry['ID']}"
        assert entry["ID"] in entry["Description"], (
            f"ID {entry['ID']} not found in Description"
        )


def test_uniprot_fasta_parsing_pattern(tmp_path, monkeypatch):
    """
    author: mkamenetskiy
    reviewer: wmavila
    category: pattern test

    Number of parsed records should equal the number of IDs provided.
    """
    monkeypatch.chdir(tmp_path)
    for accession_ids in [
        ["P15692"],
        ["P15692", "P69905"],
        ["P15692", "P69905", "P68871"],
    ]:
        structured_data, _ = protein_fetch.get_uniprot_data(accession_ids)
        assert len(structured_data) == len(accession_ids)
        for record in structured_data:
            assert record["Sequence Length"] > 0


def test_fetch_writes_accession_named_pdb_aliases(tmp_path, monkeypatch):
    """Fetch should save pair-ready accession stems instead of raw .ent names."""

    class DummyResponse:
        def __init__(self, text="", json_data=None):
            self.text = text
            self._json_data = json_data

        def raise_for_status(self):
            return None

        def json(self):
            return self._json_data

    fasta_map = {
        "P69905": (
            ">sp|P69905|HBA_HUMAN Hemoglobin subunit alpha\n"
            "MVLSPADKTN\n"
        ),
        "P68871": (
            ">sp|P68871|HBB_HUMAN Hemoglobin subunit beta\n"
            "MVHLTPEEKS\n"
        ),
    }

    json_map = {
        "P69905": {"uniProtKBCrossReferences": [{"database": "PDB", "id": "1A00"}]},
        "P68871": {"uniProtKBCrossReferences": [{"database": "PDB", "id": "1A00"}]},
    }

    def fake_get(url, timeout=10):
        accession = url.rsplit("/", 1)[-1].split(".")[0]
        if url.endswith(".fasta"):
            return DummyResponse(text=fasta_map[accession])
        return DummyResponse(json_data=json_map[accession])

    class DummyPDBList:
        def retrieve_pdb_file(self, pdb_id, pdir, file_format="pdb"):
            raw_path = Path(pdir) / f"pdb{pdb_id.lower()}.ent"
            raw_path.write_text(f"HEADER {pdb_id}\n", encoding="utf-8")
            return str(raw_path)

    monkeypatch.setattr(protein_fetch.requests, "get", fake_get)
    monkeypatch.setattr(protein_fetch, "PDBList", DummyPDBList)

    _, pdb_info = protein_fetch.get_uniprot_data(
        ["P69905", "P68871"],
        pdb_dir=str(tmp_path),
    )

    assert pdb_info == {"P69905": "1A00", "P68871": "1A00"}
    assert (tmp_path / "P69905.pdb").exists()
    assert (tmp_path / "P68871.pdb").exists()
    assert not (tmp_path / "pdb1a00.ent").exists()


def test_fetch_tolerates_missing_raw_pdb_artifact(tmp_path, monkeypatch):
    """If PDB retrieval reports an ID but no local file exists, fetch continues."""

    class DummyResponse:
        def __init__(self, text="", json_data=None):
            self.text = text
            self._json_data = json_data

        def raise_for_status(self):
            return None

        def json(self):
            return self._json_data

    fasta_map = {
        "P69905": (
            ">sp|P69905|HBA_HUMAN Hemoglobin subunit alpha\n"
            "MVLSPADKTN\n"
        ),
    }

    json_map = {
        "P69905": {"uniProtKBCrossReferences": [{"database": "PDB", "id": "1A00"}]},
    }

    def fake_get(url, timeout=10):
        accession = url.rsplit("/", 1)[-1].split(".")[0]
        if url.endswith(".fasta"):
            return DummyResponse(text=fasta_map[accession])
        return DummyResponse(json_data=json_map[accession])

    class DummyPDBList:
        def retrieve_pdb_file(self, pdb_id, pdir, file_format="pdb"):
            # Return a canonical path but intentionally do not create the file.
            return str(Path(pdir) / f"pdb{pdb_id.lower()}.ent")

    monkeypatch.setattr(protein_fetch.requests, "get", fake_get)
    monkeypatch.setattr(protein_fetch, "PDBList", DummyPDBList)

    structured_data, pdb_info = protein_fetch.get_uniprot_data(
        ["P69905"],
        pdb_dir=str(tmp_path),
    )

    assert pdb_info == {"P69905": "1A00"}
    assert structured_data and structured_data[0]["ID"].startswith("sp|P69905|")
    assert not (tmp_path / "P69905.pdb").exists()
