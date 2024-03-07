#This pipeline is for the isolates with a class II pilE at a distinct genomic locus
#This pipeline will annotate the pilE region with Bakta
#This pipeline will use the blast database of pilE alleles and look for alignments of alleles from this database to pilE, saving the alleles with sequence similarity to upload into Geneious to annotate pilE
#For this pipeline need to set the genomic region of pilE from the prokka annotation and add 1kb to each side

#!/bin/bash

######### Set the variables #########
isolate="the_name_of_the_isolate_or_experiment" #Use the same isolate name as in the seq_analysis.sh pipeline so that all of the files are named consistently
pilE_start="the_genomic_start_positon_of_pilS" #The positon in bp of the start of pilS (including fbp and lpxC +1kb on either end as a buffer)
pilE_end="the_genomic_end_position_of_pilS" #The positon in bp of the end of pilS (including fbp and lpxC +1kb on either end as a buffer)
pilE_contig="the_contig_that_pilS_is_on" #Very few of the isolates had more than one contig. If there was only one contig set contig_1. In isolates with more than one contig set the contig pilS was on. 

#Start in correct folder
path_dir=$(pwd)
db_name="pilS_nt_db" #This only needs editing if using a different database
blast_db="${path_dir}/${db_name}" #This only needs editing if using a different database

#Load modules
module load python/3.7.3
module load bedtools/2.19.1
module load samtools/1.9
module load bakta/1.5.0
module load ncbi-blast/2.12.0

#Run python script to make a query bed file that contains the contig name and start and end positions of pilE
cd ${isolate}
python make_pilE_query.py ${isolate} ${pilE_start} ${pilE_end} ${pilE_contig} ${path_dir}/${isolate}/pilE_query.bed

#Extract the pilE sequence from the whole genome fasta
bedtools getfasta -fi $path_dir/$isolate/${isolate}_corrected_consensus.fasta -bed pilE_query.bed -fo ${isolate}_class_II_pilE.fasta 

#Use Bakta to annotate the pilE region this will confirm the prokka annotation and will be uploaded to geneious to annotate the pilE region
cd bakta
bakta --output ${isolate}_class_II_pilE --genus Neisseria ${path_dir}/${isolate}/${isolate}_class_II_pilE.fasta 
cd ${path_dir}/${isolate}

#Use BLAST to look for alleles in the database with sequence similarity to the pilS region 
blastn -db blast_db -query ${isolate}_class_II_pilE.fasta -outfmt "6 sseqid qstart qend sseq" -out ${isolate}_class_II_pilE_nt_full_blast_sequences -max_target_seqs 6000

#Run python script to extract just the allele sequences that matched to pilE from the BLAST result
python ${path_dir}/extract_pilE_full_blast_seq.py ${path_dir}/${isolate}/${isolate}__class_II_pilE_nt_full_blast_sequences ${path_dir}/${isolate}/${isolate}_class_II_pilE_nt_full_blast_seq.fasta
