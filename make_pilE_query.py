
isolate = input("input isolate name ")

start_pilS = input(str("input the start of pilE "))
end_pilS = input(str("input the end of pilE "))

import sys

print ("contig_1\t", start_pilS,"\t", end_pilS)

original_stdout = sys.stdout

path = str("/NGS/active/IPL/MENINGO/analysis/paloma/2023/")

str1 = path+isolate+ "/pilE_query.bed"
query_str = "contig_1\t" + start_pilS + "\t" + end_pilS
fixed_query_str = query_str.strip()

with open(str1, 'w') as f:
        sys.stdout = f
        print (fixed_query_str)
        sys.stdout = original_stdout

