import ppinsight
from ppinsight import protein_fetch
import pytest

def test_smoke_walter():
    """
    author: wmavila
    reviewer: mkamenetskiy
    category: smoke test

    This test ensures that the get_uniprot_data() function runs without raising exceptions.
    """
    accession_ids = ["P69905", "P68871"] # human hemoglobin alpha and beta chains
    ppinsight.protein_fetch.get_uniprot_data(accession_ids, 
                                             fasta_file = "proteins.fasta",
                                             csv_file = "protein")
    return

def test_oneshot_walter():
    """
    author: wmavila
    reviewer: mkamenetskiy
    category: one-shot test

    Here, we test if get_uniprot_data() returns the expected protein data for P69905.
    This the human hemoglobin alpha chain, which we know houses 142 amino acids in its 
    sequence and corresponds to at least one structure (PDB) file.
    """
    accession_ids = ["P69905"]
    structured_data, pdb_info = ppinsight.protein_fetch.get_uniprot_data(accession_ids)
    expected_length = 142

    assert structured_data[0]["Sequence Length"] == expected_length, \
        f"Expected {expected_length}, got {structured_data[0]['Sequence Length']}"
    
    assert pdb_info["P69905"] is not None, "Expected a PDB ID for P69905"

def test_edgecase_walter():
    """
    author: wmavila
    reviewer: mkamenetskiy
    category: edge-case test

    Next, we pass a fake UniProt accession ID to get_uniprot_data(). The function should
    raise a ValueError instead of continuing silently and creating empty output folders.

    """
    accession_ids = ["fake_id"]
    with pytest.raises(ValueError) as excinfo:
        ppinsight.protein_fetch.get_uniprot_data(accession_ids)
    assert "Invalid UniProt accession ID" in str(excinfo.value)

def test_pattern_walter():
    """
    author: wmavila
    reviewer: mkamenetskiy
    category: pattern test

    Finally, we validate the following known patterns in the results returned by 
    get_uniprot_data():
    
    1. The sequence must not be empty.
    2. The ID should appear in the Description.
    """
    accession_ids = ["P69905", "P68871"]
    sequence_data, pdb_info = ppinsight.protein_fetch.get_uniprot_data(accession_ids)

    for entry in sequence_data:

        # Pattern 1: non-empty sequence
        assert entry["Sequence"], f"Empty sequence for {entry['ID']}"

        # Pattern 2: ID in description
        assert entry["ID"] in entry["Description"], \
            f"ID not found in description for {entry['ID']}"