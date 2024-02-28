## this script should take the whole genome fasta file and create a query.bed tab delimited file with the start and end positions of the sliding windows 
## specifically for genomes with many contigs

from Bio import SeqIO
import os

isolate = input("Input isolate name: ")
contig_to_process = input("Contig to process: ")  # Specify the contig you want to process
path = "/NGS/active/IPL/MENINGO/analysis/paloma/2023/"
filename = "_corrected_consensus.fasta"
deletions_folder = "deletions_across_genome"

window_size = 10000
step_size = 500

bed_file_path = os.path.join(path, deletions_folder, isolate, isolate + "_WG_windows_query.bed")

with open (bed_file_path, 'w') as bed_file:
    for record in SeqIO.parse(open(os.path.join(path, isolate, isolate + filename)), 'fasta'):
        seq_name = record.id

        # Only process the specified contig
        if seq_name != contig_to_process:
            continue

        seq_length = len(record)
        
        for start in range(0, seq_length - window_size +1, step_size):
            end = start + window_size - 1
            if end >= seq_length:
                end = seq_length - 1
                
            bed_line = f"{seq_name}\t{start}\t{end}\n"
            bed_file.write(bed_line)
