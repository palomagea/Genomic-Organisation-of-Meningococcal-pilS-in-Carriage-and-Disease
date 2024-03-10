#This script takes the whole genome fasta file and creates a query.bed tab delimited file with the start and end positions of the sliding windows 
#It works with the deletions_across_genome.sh script

from Bio import SeqIO
import os
import sys

#These variables are set in the deletions_across_genome.sh script
isolate = sys.argv[1]
whole_genome_fasta = sys.argv[2]
query_bed = sys.argv[3] 
window_size = int(sys.argv[4])
step_size = int(sys.argv[5])

#A whole genome query bed file is made that has a list of windows, each window has the contig name and start and end position of the window
#The windows span across the entire genome 
with open (query_bed, 'w') as bed_file:
    for record in SeqIO.parse(open(whole_genome_fasta, 'r'),'fasta'): #The whole genome fasta file is opened and the contig name and length of each contig is recorded 
        seq_name = record.id
        seq_length = len(record)
        
        for start in range(0, seq_length - window_size +1, step_size): #Each contig is iterated over to make windows of the set length and step size and the start and end positions of each window are recorded
            end = start + window_size - 1
            if end >= seq_length:
                end = seq_length - 1
                
            bed_line = f"{seq_name}\t{start}\t{end}\n" #For each window, the contig name, start and end positon are recorded and written to a new line of the whole genome query bed file
            bed_file.write(bed_line)
            
