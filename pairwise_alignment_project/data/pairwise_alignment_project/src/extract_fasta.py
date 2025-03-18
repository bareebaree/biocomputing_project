from Bio import SeqIO


database = "C:\\Users\\james\\Masters_Degree\\biocomputing\\research_project\\msms\\data\\dog_breeds.fa"


# Open and read the FASTA file


for record in SeqIO.parse(database, "fasta"):
    print(f"ID: {record.id}")
    print(f"Sequence: {record.seq}\n")
