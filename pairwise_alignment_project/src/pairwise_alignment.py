import random
import datetime
import multiprocessing
from Bio import SeqIO, pairwise2
from tqdm import tqdm # For progress bars

# File paths
database_path = "./data/dog_breeds.fa"
query_path = "./data/mystery.fa"

# Support for parallisation

def get_num_cores(user_input):
    """Validates and returns the number of cores to use."""
    max_cores = multiprocessing.cpu_count()
    if user_input > max_cores:
        raise ValueError(f"Error: Requested {user_input} cores, but only {max_cores} cores available.")
    return user_input

# Set number of cores (change as needed)
USER_NUM_CORES = 4  
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

# Print the best alignment
print("\nBest alignment:")
print(pairwise2.format_alignment(*alignments[0]))

# Compute percent identity
top_alignment = alignments[0]
num_matches = sum(1 for a, b in zip(top_alignment.seqA, top_alignment.seqB) if a == b)
total_length = max(len(top_alignment.seqA), len(top_alignment.seqB))
percent_identity = (num_matches / total_length) * 100

print(f"\nPercent Identity: {percent_identity:.2f}%")

# ---------------- Empirical p value calculation ----------------

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
num_cores = min(2, multiprocessing.cpu_count())  

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
timestamp = current_datetime.timestamp()

output_file = f"./results/{timestamp}_alignment_results.txt"
with open(output_file, "w") as f:
    f.write(f"Query ID: {query_record.id}\n")
    f.write(f"Query sequence: {query_seq}\n")
    f.write(f"Best match: {best_match}\n")
    f.write(f"Alignment Score: {best_score}\n")
    f.write(f"Percent Identity: {percent_identity:.2f}%\n")
    f.write(f"Empirical P-value: {empirical_p_value:.20f}\n\n")
    f.write("Alignment:\n")
    f.write(pairwise2.format_alignment(*alignments[0]))

print(f"\nResults saved to {output_file}")
