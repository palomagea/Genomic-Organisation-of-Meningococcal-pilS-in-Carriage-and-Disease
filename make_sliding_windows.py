#This script takes the whole genome fasta file and creates a query.bed tab delimited file with the start and end positions of the sliding windows 

from Bio import SeqIO
import os
import sys

isolate = sys.argv[1]
whole_genome_fasta = sys.argv[2]
query_bed = sys.argv[3] 
path = "/NGS/active/IPL/MENINGO/analysis/paloma/2023/"
filename = "_corrected_consensus.fasta"
deletions_folder = "deletions_across_genome"

window_size = 10000
step_size = 500

bed_file_path = os.path.join(path, deletions_folder, isolate, isolate + "_WG_windows_query.bed")

with open (bed_file_path, 'w') as bed_file:
    for record in SeqIO.parse(open(os.path.join(path, isolate, isolate + filename)), 'fasta'):
        seq_name = record.id
        seq_length = len(record)
        
        for start in range(0, seq_length - window_size +1, step_size):
            end = start + window_size - 1
            if end >= seq_length:
                end = seq_length - 1
                
            bed_line = f"{seq_name}\t{start}\t{end}\n"
            bed_file.write(bed_line)
            
