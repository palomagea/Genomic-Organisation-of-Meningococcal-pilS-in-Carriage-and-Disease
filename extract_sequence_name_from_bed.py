#This script works with the pilS_bed.sh pipeline to search the pilS.bed file that contains the reads mapping across pilS and save their sequence names into a pilS_read_names.txt file

import sys 

#These variables are set in the extract_pilS_annotation.sh pipeline
pilS_bed = sys.argv[1]
pilS_names_txt = sys.argv[2]

with open(pilS_bed, 'r') as bed_file:
    lines = bed_file.readlines()

sequence_names = [line.split()[3].strip() for line in lines]

with open(pilS_names_txt, 'w') as txt_file:
    txt_file.write('\n'.join(sequence_names))
