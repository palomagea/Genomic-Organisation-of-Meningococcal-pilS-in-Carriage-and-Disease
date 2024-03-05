#This script takes the whole genome fasta file and creates a query.bed tab delimited file with the start and end positions of the sliding windows 
#It works with the deletions_across_genome.sh script

from Bio import SeqIO
import os
import sys

isolate = sys.argv[1]
whole_genome_fasta = sys.argv[2]
query_bed = sys.argv[3] 


window_size = sys.argv[4]
step_size = sys.argv[5]


with open (query_bed, 'w') as bed_file:
    for record in SeqIO.parse(open(whole_genome_fasta, 'fasta'):
        seq_name = record.id
        seq_length = len(record)
        
        for start in range(0, seq_length - window_size +1, step_size):
            end = start + window_size - 1
            if end >= seq_length:
                end = seq_length - 1
                
            bed_line = f"{seq_name}\t{start}\t{end}\n"
            bed_file.write(bed_line)
            
