import os
import platform
import subprocess

# Define paths
DB_NAME = "dna_database"
DB_FASTA = "./dog_breeds.fa"
QUERY_FASTA = "./mystery.fa"
OUTPUT_FILE = "./results.txt"

# Detect the OS
is_windows = platform.system() == "Windows"

# Set BLAST command format
makeblastdb_cmd = ["makeblastdb", "-in", DB_FASTA, "-dbtype", "nucl", "-out", DB_NAME]
blastn_cmd = ["blastn", "-query", QUERY_FASTA, "-db", DB_NAME, "-out", OUTPUT_FILE, "-outfmt", "6"]

def run_command(command):
    """Run a system command and handle errors."""
    try:
        result = subprocess.run(command, shell=is_windows, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(result.stdout.decode())
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {command}\n{e.stderr.decode()}")
        exit(1)

# Run BLAST database creation
print("Creating BLAST database...")
run_command(makeblastdb_cmd)

# Run BLAST search
print("Running BLAST search...")
run_command(blastn_cmd)

print(f"BLAST search complete. Results saved in {OUTPUT_FILE}")
