#This script works with the extract_pilS_annotation.sh pipeline to create a query file for the pilS region
#The query file contains the contig pilS is on and the start and end positions of pilS (pilS includes both lpxC and fbp +1kb on either side as a buffer)

import sys 

#These variables are set in the extract_pilS_annotation.sh pipeline
isolate = sys.argv[1]
pilS_start = sys.argv[2]
pilS_end = sys.argv[3] 
pilS_contig = sys.argv[4]
query_bed = sys.argv[5]

original_stdout = sys.stdout

path = str("/NGS/active/IPL/MENINGO/analysis/paloma/2023/")

str1 = path+isolate+ "/query.bed"
query_str = "contig_1\t" + start_pilS + "\t" + end_pilS
fixed_query_str = query_str.strip()

with open(str1, 'w') as f:
        sys.stdout = f
        print (fixed_query_str)
        sys.stdout = original_stdout
