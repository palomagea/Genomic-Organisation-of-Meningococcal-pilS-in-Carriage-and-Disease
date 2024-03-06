#This script works with the extract_pilS_annotation.sh pipeline to create a genome file containing the length of the whole genome

import sys 

#These variables are set in the extract_pilS_annotation.sh pipeline
genome_length = sys.argv[1]

with open(genome_length, 'r') as openfile:
    rfile = openfile.readlines()
    length = rfile[1].strip()  

output_filename = path + isolate + "/" + isolate + ".genome"

with open(output_filename, 'w') as f:
    f.write(f"contig_1\t{length}\n") 
