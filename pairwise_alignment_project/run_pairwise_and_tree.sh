#!/bin/bash

# Exit script on any error


echo "Running pairwise alignment..."
python3 src/pairwise_alignment.py

echo "Drawing phylogenetic tree..."
python3 src/draw_phylo_tree.py

echo "All tasks completed successfully!"
