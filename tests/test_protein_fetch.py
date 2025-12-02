import ppinsight
from ppinsight import protein_fetch


def test_smoke_walter():
    """
    Walter's smoke test: ensure the get_uniprot_data() function runs w/o raising exceptions.
    """
    accession_ids = ["P69905", "P68871"]
    try:
        structured_data, pdb_info = ppinsight.protein_fetch.get_uniprot_data(
            accession_ids,
            fasta_file="test_sequences.fasta",
            csv_file="test_sequences.csv",
            pdb_dir="test_pdb_files"
        )
        print("Smoke test passed: function ran w/o exceptions.")
    except Exception as e:
        print(f"Smoke test failed: unexpected exception {e}")

#def test_oneshot_walter():
 #   """
  #  Walter's oneshot test: 