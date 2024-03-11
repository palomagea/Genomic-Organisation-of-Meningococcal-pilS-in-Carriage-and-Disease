#This script works with the class_II_pilE_extraction.sh pipeline to create a query file for the pilE region of isolates with a class II pilE
#The query file contains the contig pilE is on and the start and end positions of pilE (pilS includes both katA and prlC +1kb on either side as a buffer)

import sys 

#These variables are set in the extract_pilS_annotation.sh pipeline
isolate = sys.argv[1]
pilE_start = sys.argv[2]
pilE_end = sys.argv[3] 
pilE_contig = sys.argv[4]
query_bed = sys.argv[5]

original_stdout = sys.stdout

query_str = pilE_contig +"\t" + pilE_start + "\t" + pilE_end
fixed_query_str = query_str.strip()

with open(query_bed, 'w') as f:
        sys.stdout = f
        print (fixed_query_str)
        sys.stdout = original_stdout

