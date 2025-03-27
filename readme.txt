# This project sets out to:

Identify the most similar sequence in a dataset to an input sequence.

# It functions by:

1. Iteratively performing pairwise sequence analysis using Biopython, and choosing the highest scoring sequence.

2. Computing probabilities across dataset.

3. Reconstructing a phylogenetic tree with the dataset and input sequence.

To run:

1. Download entire project folder 
2. Ensure python 3.12 is installed
3. ensure following libraries are installed:
  -BioPython
  -Matplotlib
  -Pytest
  -datetime
  -tqdm
4. chmod +x run_pairwise_and_tree.sh to 
5. (optional if on windows using WSL) - may be necessary to use
- dos2unix run_pairwise_and_tree.sh 
to convert it to unix.

run the shell script  run_pairwise_and_tree.sh

If this doesn't print:

1. Best match

2. P-value of best match

3. Save the phylogenetic tree

Try running in an IDE, or the command line the two python scripts sequentially:

1. pairwise_alignment.py
2. draw_phylo_tree.py

Filepaths given are relative to their position in this folder, etc, ./data/...

But if they do not path on your system for whatever reason, it should work if you use the absolute filepath.



