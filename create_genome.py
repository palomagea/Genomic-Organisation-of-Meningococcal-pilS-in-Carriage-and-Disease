#This script works with the extract_pilS_annotation.sh pipeline to create a genome file containing the length of the whole genome

isolate = input("input isolate name ")
path = "/NGS/active/IPL/MENINGO/analysis/paloma/2023/"
filename = "genome_length.txt"

with open(path + isolate + "/" + filename, 'r') as openfile:
    rfile = openfile.readlines()
    length = rfile[1].strip()  

output_filename = path + isolate + "/" + isolate + ".genome"

with open(output_filename, 'w') as f:
    f.write(f"contig_1\t{length}\n") 
