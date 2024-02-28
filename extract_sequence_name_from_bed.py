isolate = input("input isolate name ")
path = "/NGS/active/IPL/MENINGO/analysis/paloma/2023/"
filename = "_pilS.bed"


input_bed_file = path + isolate + "/" + isolate + filename
output_txt_file = 'pilS_read_names.txt'

with open(input_bed_file, 'r') as bed_file:
    lines = bed_file.readlines()

sequence_names = [line.split()[3].strip() for line in lines]

with open(output_txt_file, 'w') as txt_file:
    txt_file.write('\n'.join(sequence_names))

