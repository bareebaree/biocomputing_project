from Bio.Blast.Applications import NcbiblastnCommandline

blastn_cline = NcbiblastnCommandline(query="query.fasta", db="mydb",
                                     out="blast_results.xml", outfmt=5)
stdout, stderr = blastn_cline()
