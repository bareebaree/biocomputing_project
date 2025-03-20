# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 10:51:08 2025

@author: james
"""
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
import matplotlib.pyplot as plt
import matplotlib

# Load alignment
alignment = AlignIO.read(
    "C:\\Users\\james\\Masters_Degree\\biocomputing\\research_project\\pairwise_alignment_project\\data\\20250318_232936_msa_output.aln",
    "clustal",
)
print(alignment)

# Compute distance matrix
calculator = DistanceCalculator("identity")
dm = calculator.get_distance(alignment)
print(dm)

# Construct tree
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)  # Use UPGMA instead: constructor.upgma(dm)

# Save tree in Newick format
Phylo.write(tree, "phylogenetic_tree.nwk", "newick")

def remove_inner_labels(tree):
    for clade in tree.find_clades():
        if not clade.is_terminal():  # Internal node
            clade.name = None  # Remove "Inner" labels
    return tree

# Apply pruning function
tree = remove_inner_labels(tree)

# Fix the plot function (Removed incorrect StringIO conversion)
def plot_tree(tree):
    matplotlib.rc("font", size=6)
    fig, ax = plt.subplots(figsize=(10, 20), dpi=100)
    Phylo.draw(tree, axes=ax)  # ✅ Pass tree object directly
    plt.show()

# Display ASCII representation
Phylo.draw_ascii(tree)

# Plot the tree
plot_tree(tree)
