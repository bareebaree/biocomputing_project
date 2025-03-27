
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from Bio import SeqIO
from phylo_tree_functions import remove_inner_labels, plot_tree





# Define file paths. These are default filepaths. If you need to do a different query or dataset, replace the filepaths with correct directories.
# The database uses iterative pairwise alignment so does not scale particularly well for larger tasks. If larger databases, or longer sequences are desired,
# it is recommended to use a MSA tool such as MAFFT or ClustalOmega.
database_path = "./data/dog_breeds.fa"
query_path = "./data/mystery.fa"




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
