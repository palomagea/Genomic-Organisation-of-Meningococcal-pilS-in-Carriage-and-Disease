path = "/NGS/active/IPL/MENINGO/analysis/paloma/2023/"
isolate = input("input isolate name ")
filename = "/"+isolate+"_pilS_nt_full_blast_sequences"

import pandas as pd
df = pd.read_csv(path+isolate+filename, delimiter = "\t", header=None)
df.columns = ["sseqid", "qstart", "qend", "sseq"]


import sys
original_stdout = sys.stdout
str1 = path+isolate+ "/" + isolate + "_pilS_nt_full_blast_seq.fasta"

with open(str1, 'w') as f:
        sys.stdout = f
        for index, row in df.iterrows():
            seq_id = row["sseqid"]
            sequence = row["sseq"]
            print (">", seq_id, "\n", sequence, sep = '')
        sys.stdout = original_stdout
