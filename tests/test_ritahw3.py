"""
Tests for the test_protein_fetch function
"""
import pytest
import numpy as np

from ppinsight import protein_fetch

# Smoke test- ensure that the function runs 

accession_list = ['P04637', 'P68871', 'Q8NEC1']
ppinsight.protein_fetch.get_uniprot_data(accession_list, fasta_file="proteins.fasta", csv_file="protein")

# one-shot test (1 point)
#edge test (1 point)
#pattern test (1 point)