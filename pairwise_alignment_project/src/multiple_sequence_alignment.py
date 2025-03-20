import datetime
from Bio import SeqIO
from tqdm import tqdm  # For progress bars
from Bio.Align.Applications import ClustalOmegaCommandline


### append query sequence to database

# Define file paths
database_path = "./data/dog_breeds.fa"
query_path = "./data/mystery.fa"
combined_path= "./data/combined_sequences.fa"



# Read sequences from both files
database_seqs = list(SeqIO.parse(database_path, "fasta"))
query_seqs = list(SeqIO.parse(query_path, "fasta"))

# Append the sequences
combined_seqs = database_seqs + query_seqs

# Write the combined sequences to a new file
SeqIO.write(combined_seqs, combined_path, "fasta")

print(f"Sequences combined and saved to: {combined_path}")

### Perform multiple sequence alignment on database

current_datetime = datetime.datetime.now()
timestamp = current_datetime.strftime("%Y%m%d_%H%M%S")


# Define Clustal Omega command
clustalomega_cline = ClustalOmegaCommandline(infile=combined_path,
                                             outfile=f"./data/{timestamp}_msa_output.aln",
                                             seqtype="dna",
                                             verbose=True, outfmt="clustal", auto=True)

# force the skibidi Clustal Omega sub to work for you
stdout, stderr = clustalomega_cline()
print("Alignment completed!")
