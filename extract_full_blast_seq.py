#This script works with the extract_pilS_annotation.sh pipeline to extract the sequences of the PubMLST alleles from the database that align to the pilS region
#The fasta file of these sequences will later be uploaded to geneious to annotate the pilS region

import sys 

#These variables are set in the extract_pilS_annotation.sh pipeline
blast_seq_file = sys.argv[1]
output_file = sys.argv[2]

import pandas as pd
df = pd.read_csv(blast_seq_file, delimiter = "\t", header=None)
df.columns = ["sseqid", "qstart", "qend", "sseq"]


import sys
original_stdout = sys.stdout

with open(output_file, 'w') as f:
        sys.stdout = f
        for index, row in df.iterrows():
            seq_id = row["sseqid"]
            sequence = row["sseq"]
            print (">", seq_id, "\n", sequence, sep = '')
        sys.stdout = original_stdout
