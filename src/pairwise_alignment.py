from Bio import SeqIO
from Bio import pairwise2

# File paths
database_path = "./data/dog_breeds.fa"
query_path = "./data/mystery.fa"

# Read the FASTA database into a dictionary
database = {record.id: str(record.seq) for record in SeqIO.parse(database_path, "fasta")}

# Read the query sequence (assuming only one sequence in the file)
query_record = next(SeqIO.parse(query_path, "fasta"))  

# Convert to string for alignment
query_seq = str(query_record.seq)  
best_match = None
best_score = float('-inf')

# Find the most similar sequence
for name, seq in database.items():
    score = pairwise2.align.globalxx(query_seq, seq, score_only=True)
    
    if score > best_score:
        best_score = score
        best_match = name

# Retrieve the best matching sequence
best_match_seq = database[best_match]

print(f"Best match: {best_match} with score: {best_score}")

# Perform full sequence alignment
alignments = pairwise2.align.globalxx(query_seq, best_match_seq)

# Print the best alignment
print("Best alignment:")
print(pairwise2.format_alignment(*alignments[0]))

# Compute percent identity
top_alignment = alignments[0]
num_matches = sum(1 for a, b in zip(top_alignment.seqA, top_alignment.seqB) if a == b)
total_length = max(len(top_alignment.seqA), len(top_alignment.seqB))
percent_identity = (num_matches / total_length) * 100

print(f"Percent Identity: {percent_identity:.2f}%")

# Save results to a file
output_file = "alignment_results.txt"
with open(output_file, "w") as f:
    f.write(f"Query ID: {query_record.id}\n")
    f.write(f"Query sequence: {query_seq}\n")
    f.write(f"Best match: {best_match}\n")
    f.write(f"Alignment Score: {best_score}\n")
    f.write(f"Percent Identity: {percent_identity:.2f}%\n\n")
    f.write("Alignment:\n")
    f.write(pairwise2.format_alignment(*alignments[0]))

print(f"Results saved to {output_file}")
