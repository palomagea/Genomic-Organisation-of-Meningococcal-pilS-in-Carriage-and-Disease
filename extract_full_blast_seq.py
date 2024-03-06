#This script works with the extract_pilS_annotation.sh pipeline to extract the sequences of the PubMLST alleles from the database that align to the pilS region
#The fasta file of these sequences will later be uploaded to geneious to annotate the pilS region

import sys 
import pandas as pd

#These variables are set in the extract_pilS_annotation.sh pipeline
blast_seq_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_csv(blast_seq_file, delimiter = "\t", header=None) #read the blast output into a dataframe
df.columns = ["sseqid", "qstart", "qend", "sseq"] #name the columns of the dataframe so can call the sequence id names and the sequences below

original_stdout = sys.stdout

with open(output_file, 'w') as f: #create the output file
        sys.stdout = f
        for index, row in df.iterrows(): #print out all of the sequence id and sequences of the alleles that align to pilS in fasta format 
            seq_id = row["sseqid"]
            sequence = row["sseq"]
            print (">", seq_id, "\n", sequence, sep = '')
        sys.stdout = original_stdout
