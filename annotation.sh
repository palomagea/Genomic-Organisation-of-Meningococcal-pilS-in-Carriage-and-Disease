#!/bin/bash
#this pipeline will do all the annotation steps you need to upload into geneious to make the figures 
#this includes annotating pilS with bakta and making the blast db to annotate with blast in geneious

#start in correct folder 
echo Enter pwd for the folder you want to be working in 
read path_dir

#set isolate name as variable 
echo Enter isolate name
read isolate

#do bakta run on pilS fasta
cd ${path_dir}/${isolate}
mkdir bakta
cd bakta
module load bakta/1.5.0
bakta --output ${isolate}_pilS --genus Neisseria ${path_dir}/${isolate}/${isolate}_pilS.fasta 

#do the limited blast for nucleotide sequences with the longest matches 
module load ncbi-blast/2.12.0
cd ${path_dir}/${isolate}
blastn -db /NGS/active/IPL/MENINGO/analysis/paloma/pilS_nt_db/pilS_nt_db -query ${isolate}_pilS.fasta -outfmt "6 sseqid qstart qend sseq" -out ${isolate}_pilS_nt_limited_blast_sequences -max_target_seqs 6000 -max_hsps 1  -culling_limit 1 

#full blast to get all sequences that align  
blastn -db /NGS/active/IPL/MENINGO/analysis/paloma/pilS_nt_db/pilS_nt_db -query ${isolate}_pilS.fasta -outfmt "6 sseqid qstart qend sseq" -out ${isolate}_pilS_nt_full_blast_sequences -max_target_seqs 6000

#run python script to extract the sequences from the limited blast  
module load python/3.7.3
python ${path_dir}/extract_limited_blast_seq.py

#run python script to get the sequences from the full blast   
python ${path_dir}/extract_full_blast_seq.py

