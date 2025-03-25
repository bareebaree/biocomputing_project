# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 10:51:08 2025

@author: james
"""

from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
import matplotlib.pyplot as plt
import matplotlib
from Bio import SeqIO


def remove_inner_labels(tree) -> Phylo.BaseTree.Tree:
    """ 
    Function to remove inner labels for nodes in tree, by selecting labels, then replacing value with None to make tree clearer to read.
    
    PARAMETERS: tree -> (Phylo.BaseTree.Tree) tree object
    
    RETURNS: tree -> tree object with non terminal labels removed.
    
    """
    for clade in tree.find_clades():
        if not clade.is_terminal():  # Internal node
            clade.name = None  # Remove "Inner" labels
    return tree



def plot_tree(tree) -> Phylo.BaseTree.Tree:
    """
    "Function that calls matplotlib and Phylo.draw to create a phylogenetic tree plot"

    PARAMETERS: tree -> (Phylo.BaseTree.Tree) tree object

    RETURNS: 
    None.
    Plots the tree.

    """
    
    matplotlib.rc("font", size=6)
    fig, ax = plt.subplots(figsize=(10, 20), dpi=100)
    Phylo.draw(tree, axes=ax)  # ✅ Pass tree object directly
    plt.show()
    
### append query sequence to database


# Define file paths
database_path = "C:\\Users\\james\\Masters_Degree\\biocomputing\\research_project\\pairwise_alignment_project\\data\\dog_breeds.fa"
query_path = "C:\\Users\\james\\Masters_Degree\\biocomputing\\research_project\\pairwise_alignment_project\\data\\mystery.fa"




# Read sequences from both files
database_seqs = list(SeqIO.parse(database_path, "fasta"))
query_seqs = list(SeqIO.parse(query_path, "fasta"))

# Append sequences to form a single alignment list
alignment_seqs = database_seqs + query_seqs

# Create MultipleSeqAlignment object
alignment = MultipleSeqAlignment(alignment_seqs)

# Compute distance matrix
calculator = DistanceCalculator("identity")
dm = calculator.get_distance(alignment)
print(dm)

# Construct tree
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)

# Save tree in Newick format
Phylo.write(tree, "phylogenetic_tree.nwk", "newick")



# Apply pruning function
tree = remove_inner_labels(tree)


# Display ASCII representation
Phylo.draw_ascii(tree)

# Plot the tree
plot_tree(tree)
