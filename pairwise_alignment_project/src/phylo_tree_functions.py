from Bio import Phylo
import matplotlib.pyplot as plt
from Bio.Phylo.BaseTree import Tree


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




def plot_tree(tree: Tree, font_size: int = 6, figsize=(10, 20), dpi: int = 100) -> None:
    """
    Plots a phylogenetic tree using matplotlib.

    Parameters:
        tree (Tree): A Biopython Tree object.
        font_size (int): Font size for labels.
        figsize (tuple): Figure size in inches.
        dpi (int): Dots per inch for rendering.

    Returns:
        None. Displays the plot.
    """
    plt.rc("font", size=font_size)
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    Phylo.draw(tree, axes=ax)
    plt.show()### append query sequence to database