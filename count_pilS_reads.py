#This script works with the extract_pilS_annotation.sh pipeline to count the number of reads overlapping the whole pilS region 
#The pilS bed file only contains reads that map across the entire pilS region 

import sys 

#These variables are set in the extract_pilS_annotation.sh pipeline
pilS_bed = sys.argv[1]

pilS_list = open (pilS_bed, 'r') #open the pilS bed file

count = 0 

for line in pilS_list: #each line in the pilS bed file is a read that maps across pilS so for each line the count increases by 1 

        if line != "\n": 

                count += 1 

print ("The total number of reads fully mapping across PilS are ", count) #the total number of reads mapping the whole way across pilS is printed
