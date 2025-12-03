"""
Tests for the test_protein_fetch function
"""
import pytest


from ppinsight import protein_fetch

# Smoke test- ensure that the function runs

accession_list_simple = ['P15692']
accession_list_3 = ['P04637', 'P68871', 'Q8NEC1']
# This is the function to call:
# protein_fetch.get_uniprot_data(accession_list, fasta_file="proteins.fasta",
# csv_file="protein")


def test_smoke_rita():
    """
    author: mkamenetskiy
    reviewer: wmavila
    category: smoke test

    Simple smoke test to ensure function runs
    """
    protein_fetch.get_uniprot_data(accession_list_simple,
                                   fasta_file="proteins.fasta",
                                   csv_file="protein")
    return


def test_get_uniprot_data_p15692_structured():
    """
    author: mkamenetskiy
    reviewer: wmavila
    category: one-shot test

    One-shot test with a known UniProt accession P15692 (VEGFA_HUMAN).
    Expected structured data should match exactly!
    """
    expected_structured_data = [{
        'ID': 'sp|P15692|VEGFA_HUMAN',
        'Name': 'sp|P15692|VEGFA_HUMAN',
        'Description': 'sp|P15692|VEGFA_HUMAN Vascular endothelial growth '
        'factor A, long form OS=Homo sapiens OX=9606 GN=VEGFA PE=1 SV=3',
        'Sequence Length': 395,
        'Sequence': 'MTDRQTDTAPSPSYHLLPGRRRTVDAAASRGQGPEPAPGGGVEGVGARGV'
        'ALKLFVQLLGCSRFGGAVVRAGEAEPSGAARSASSGREEPQPEEGEEEEEKEEERGPQWRLG'
        'ARKPGSWTGEAAVCADSAPAARAPQALARASGRGGRVARRGAEESGPPHSPSRRGSASRAGPG'
        'RASETMNFLLSVHWSLALLLYLHHAKWSQAAPMAEGGGQNHHEVVKFMDVYQRSYCHPIETLV'
        'DIFQEYPDEIEYIFKPSCVPLMRCGGCCNDEGLECVPTEESNITMQIMRIKPHQGQHIGEMSFL'
        'QHNKCECRPKKDRARQEKKSVRGKGKGQKRKRKKSRYKSWSVPCGPCSERRKHLFVQDPQTCKC'
        'SCKNTDSRCKARQLELNERTCRCDKPRR'
    }]
    actual_structured_data, actual_pdb_info = protein_fetch.get_uniprot_data(
        ["P15692"])

    assert actual_structured_data == expected_structured_data


def test_get_uniprot_data_p15692_pdbinfo():
    """
    author: mkamenetskiy
    reviewer: wmavila
    category: one-shot test

    One-shot test using known UniProt accession P15692 (VEGFA_HUMAN).
    Expected PDB info should match exactly!
    """
    expected_pdb_info = {'P15692': '1BJ1'}
    actual_structured_data, actual_pdb_info = protein_fetch.get_uniprot_data(
        ["P15692"])

    assert actual_pdb_info == expected_pdb_info


def test_nonexistent_protein():
    """
    author: mkamenetskiy
    reviewer: wmavila
    category: edge test

    When get_uniprot_data is called with a protein that doesn't exist or that
    can't be downloaded, I want get_uniprot_data to do:
    - Say that that protein ID was not correct within the Uniprot DataBase
    - I will want the function to stop and error out
        - raise some type of error saying that the value of the ID was bad
            - ValueError makes sense
    - Even if there are other protein IDs that are correct, I want the
    function to stop
    """
    with pytest.raises(ValueError, match="Could not fetch FASTA for"):
        protein_fetch.get_uniprot_data(["P000000000"])


def test_uniprot_fasta_parsing_pattern():
    """
    author: mkamenetskiy
    reviewer: wmavila
    category: pattern test

    This will be a pattern test to ensure the number of parsed records equals
    the number of IDs that were provided.
    """
    def test(accession_ids):
        structured_data, _ = protein_fetch.get_uniprot_data(accession_ids)
        assert len(structured_data) == len(accession_ids)
        for record in structured_data:
            assert record["Sequence Length"] > 0  # Not empty seq
        return

    # Run for different sizes of accession lists
    for accession_ids in [["P15692"], ["P15692", "P69905"],
                          ["P15692", "P69905", "P68871"]]:
        test(accession_ids)
