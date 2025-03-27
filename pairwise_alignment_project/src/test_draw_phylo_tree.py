# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 00:28:55 2025

@author: james
"""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from phylo_tree_functions import remove_inner_labels, plot_tree  

# Mock sequences for alignment
@pytest.fixture
def mock_alignment():
    seqs = [
        SeqRecord(Seq("ACTG"), id="A"),
        SeqRecord(Seq("ACTA"), id="B"),
        SeqRecord(Seq("ACTC"), id="C")
    ]
    return MultipleSeqAlignment(seqs)

# Test alignment creation
def test_alignment_shape(mock_alignment):
    assert isinstance(mock_alignment, MultipleSeqAlignment)
    assert len(mock_alignment) == 3  # 3 sequences
    assert mock_alignment.get_alignment_length() == 4  # each is length 4

# Test distance matrix generation
def test_distance_matrix(mock_alignment):
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(mock_alignment)
    assert dm is not None
    assert dm["A", "B"] >= 0

# Test tree construction
def test_tree_construction(mock_alignment):
    calculator = DistanceCalculator("identity")
    
    dm = calculator.get_distance(mock_alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    assert isinstance(tree, Phylo.BaseTree.Tree)
    assert len(tree.get_terminals()) == 3

# Test internal label removal
def test_remove_inner_labels(mock_alignment):
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(mock_alignment)
    tree = DistanceTreeConstructor().nj(dm)
    tree_with_removed = remove_inner_labels(tree)

    for clade in tree_with_removed.get_nonterminals():
        assert clade.name is None

