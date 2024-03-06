#This script works with the extract_pilS_annotation.sh pipeline to create a genome file containing the length of the whole genome

import sys 

#These variables are set in the extract_pilS_annotation.sh pipeline
genome_length = sys.argv[1]
output_filename = sys.argv[2] 

with open(genome_length, 'r') as openfile: #The txt file containing the genome length is opened
    rfile = openfile.readlines()
    length = rfile[1].strip()  #The length of the genome is saved as the length variable 

with open(output_filename, 'w') as f: #A genome file is made with the length of the genome 
    f.write(f"contig_1\t{length}\n") 
