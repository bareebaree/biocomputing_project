import random
import datetime
import multiprocessing
from Bio import SeqIO, pairwise2
from tqdm import tqdm  # For progress bars
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord





# File paths
database_path= "C:\\Users\\james\\Masters_Degree\\biocomputing\\research_project\\pairwise_alignment_project\\data\\dog_breeds.fa"
query_path = "C:\\Users\\james\\Masters_Degree\\biocomputing\\research_project\\pairwise_alignment_project\\data\\mystery.fa"

# Support for parallelization
def get_num_cores(user_input):
    """Validates and returns the number of cores to use."""
    max_cores = multiprocessing.cpu_count()
    if user_input > max_cores:
        raise ValueError(f"Error: Requested {user_input} cores, but only {max_cores} cores available.")
    return user_input

# Set number of cores (change as needed)
USER_NUM_CORES = 2  
num_cores = get_num_cores(USER_NUM_CORES)  

print(f"Using {num_cores} CPU cores for parallel computation...")

# Read the FASTA database into a dictionary
database = {record.id: str(record.seq) for record in SeqIO.parse(database_path, "fasta")}

# Read the query sequence (assuming only one sequence in the file)
query_record = next(SeqIO.parse(query_path, "fasta"))  
query_seq = str(query_record.seq)  

best_match = None
best_score = float('-inf')

print("Searching for the best match in the database...")

# Find the most similar sequence (with progress bar)
for name, seq in tqdm(database.items(), desc="Aligning sequences", unit="seq"):
    score = pairwise2.align.globalxx(query_seq, seq, score_only=True)
    
    if score > best_score:
        best_score = score
        best_match = name

# Retrieve the best matching sequence
best_match_seq = database[best_match]

print(f"\nBest match: {best_match} with score: {best_score}")

# Perform full sequence alignment
alignments = pairwise2.align.globalxx(query_seq, best_match_seq)

# Get the best alignment
best_alignment = alignments[0]
aligned_query_seq = best_alignment.seqA
aligned_best_match_seq = best_alignment.seqB

# Print the best alignment
print("\nBest alignment:")
print(pairwise2.format_alignment(*best_alignment))

# Compute percent identity
num_matches = sum(1 for a, b in zip(aligned_query_seq, aligned_best_match_seq) if a == b)
total_length = max(len(aligned_query_seq), len(aligned_best_match_seq))
percent_identity = (num_matches / total_length) * 100

print(f"\nPercent Identity: {percent_identity:.2f}%")

# ---------------- Empirical P-value Calculation ----------------

def shuffle_sequence(seq):
    """Shuffle a sequence while preserving its composition."""
    seq_list = list(seq)
    random.shuffle(seq_list)
    return "".join(seq_list)

def shuffle_and_align(_):
    """Shuffles the query sequence and computes an alignment score."""
    shuffled_query = shuffle_sequence(query_seq)
    return pairwise2.align.globalxx(shuffled_query, best_match_seq, score_only=True)

num_shuffles = 100  # Number of shuffled alignments
num_cores = min(4, multiprocessing.cpu_count())  

print("\nComputing empirical P-value with shuffled sequences...")

# Progress bar for shuffled sequences
random_scores = []
for _ in tqdm(range(num_shuffles), desc="Shuffling and aligning", unit="shuffle"):
    shuffled_query = shuffle_sequence(query_seq)
    score = pairwise2.align.globalxx(shuffled_query, best_match_seq, score_only=True)
    random_scores.append(score)

# Compute empirical P-value
empirical_p_value = sum(s >= best_score for s in random_scores) / num_shuffles

print(f"\nEmpirical P-value: {empirical_p_value:.8f}")

# ---------------- Save results ----------------

current_datetime = datetime.datetime.now()
timestamp = current_datetime.strftime("%Y%m%d_%H%M%S")

# Create SeqRecords for the shuffled aligned query and the best-match sequence
shuffled_query_record_aligned = SeqRecord(Seq(shuffled_query), id=f"{query_record.id}_shuffled", description="Shuffled aligned query sequence")
best_match_record_aligned = SeqRecord(Seq(aligned_best_match_seq), id=best_match, description="Aligned best match sequence")

# Generate filename for the shuffled alignment FASTA file
shuffled_alignment_fasta_file = f"C:\\Users\\james\\Masters_Degree\\biocomputing\\research_project\\pairwise_alignment_project\\results\\{timestamp}_shuffled_alignment.fasta"

# Write shuffled aligned query and aligned best-match sequences to FASTA
with open(shuffled_alignment_fasta_file, "w") as fasta_out:
    SeqIO.write([shuffled_query_record_aligned, best_match_record_aligned], fasta_out, "fasta")

print(f"\nShuffled alignment saved as FASTA file: {shuffled_alignment_fasta_file}")

output_file = f"C:\\Users\\james\\Masters_Degree\\biocomputing\\research_project\\pairwise_alignment_project\\results\\{timestamp}_alignment_results.txt"
alignment_fasta_file = f"C:\\Users\\james\\Masters_Degree\\biocomputing\\research_project\\pairwise_alignment_project\\results\\{timestamp}_alignment.fasta"

# Save results in a text file
with open(output_file, "w") as f:
    f.write(f"Query ID: {query_record.id}\n")
    f.write(f"Query sequence: {query_seq}\n")
    f.write(f"Best match: {best_match}\n")
    f.write(f"Alignment Score: {best_score}\n")
    f.write(f"Percent Identity: {percent_identity:.2f}%\n")
    f.write(f"Empirical P-value: {empirical_p_value:.20f}\n\n")
    f.write("Alignment:\n")
    f.write(pairwise2.format_alignment(*best_alignment))

print(f"\nResults saved to {output_file}")

# ---------------- Save Alignment as FASTA ----------------

# Create SeqRecords for aligned sequences
query_record_aligned = SeqRecord(Seq(aligned_query_seq), id=query_record.id, description="Aligned query sequence")
best_match_record_aligned = SeqRecord(Seq(aligned_best_match_seq), id=best_match, description="Aligned best match sequence")

# Write aligned sequences to FASTA
with open(alignment_fasta_file, "w") as fasta_out:
    SeqIO.write([query_record_aligned, best_match_record_aligned], fasta_out, "fasta")

print(f"\nAlignment saved as FASTA file: {alignment_fasta_file}")
