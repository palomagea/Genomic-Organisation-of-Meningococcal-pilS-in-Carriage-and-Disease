#This pipeline will extract the reads that map fully across the pilS region and count them
#This pipeline will then annotate pilS with bakta
#This pipeline will use the blast database of pilS alleles and look for alignments of alleles from this database to pilS, saving the alleles with sequence similarity to upload into Geneious to annotate pilS
#For this pipeline need to set the genomic region of pilS from the prokka annotation of fbp and lpxC and add 1kb to each side

#!/bin/bash

######### Set the variables #########
isolate="the_name_of_the_isolate_or_experiment" #Use the same isolate name as in the seq_analysis.sh pipeline so that all of the files are named consistently
pilS_start="the_genomic_start_positon_of_pilS" #The positon in bp of the start of pilS (including fbp and lpxC +1kb on either end as a buffer)
pilS_end="the_genomic_end_position_of_pilS" #The positon in bp of the end of pilS (including fbp and lpxC +1kb on either end as a buffer)
pilS_contig="the_contig_that_pilS_is_on" #Very few of the isolates had more than one contig. If there was only one contig set "contig_1". In isolates with more than one contig set the contig pilS was on. 
#Set the working directory with the scripts saved in it 
path_dir=$(pwd)
db_name="pilS_nt_database

#Load Modules 
module load python/3.7.3
module load bedtools/2.31.0
module load samtools/1.9
module load bakta/1.5.0
module load ncbi-blast/2.6.0

#Run python script to make a query bed file that contains the contig name and start and end positions of pilS
cd ${isolate}
python ${path_dir}/make_query.py ${isolate} ${pilS_start} ${pilS_end} ${pilS_contig} ${path_dir}/${isolate}/query.bed

#Make a genome file that contains the length of the genome
awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' $path_dir/$isolate/${isolate}_corrected_consensus.fasta > $path_dir/$isolate/genome_length.txt #Get the length of genome and save to a txt file
python ${path_dir}/create_genome.py ${path_dir}/${isolate}/genome_length.txt #Run python script to make a genome file from the txt file with the genome length

#Extract only the reads that map the whole way across the pilS region
bedtools intersect -a ${isolate}_allreads.bed -b query.bed -F 1 > ${isolate}_pilS.bed #the query bed file contains the pilS region and the intersect function extracts reads from the allreads.bed file that overlap the whole pilS region

#Run python script that counts reads mapping pilS from the pilS.bed file 
python ${path_dir}/count_pilS_reads.py ${path_dir}/${isolate}/${isolate}_pilS.bed #this script prints the number of reads mapping across the entire pilS 

#Extract the pilS sequence from the whole genome fasta
bedtools getfasta -fi $path_dir/$isolate/${isolate}_corrected_consensus.fasta -bed query.bed -fo ${isolate}_pilS.fasta #The query bed file containing the pilS region is used to extract the pilS fasta from the whole genome fasta

#Use Bakta to annotate the pilS region
mkdir bakta #make a directory for the Bakta annotation 
cd bakta
bakta --output ${isolate}_pilS --genus Neisseria ${path_dir}/${isolate}/${isolate}_pilS.fasta #Set the genus as Neisseria

#Use BLAST to look for alleles in the database with sequence similarity to the pilS region 
blastn -db blast_db -query ${isolate}_pilS.fasta -outfmt "6 sseqid qstart qend sseq" -out ${isolate}_pilS_nt_full_blast_sequences -max_target_seqs 6000 #Outfmt includes the sequence ID and the sequence of the alleles that matched to pilS

#Run python script to extract just the allele sequences that matched to pilS from the BLAST result
python ${path_dir}/extract_full_blast_seq.py ${path_dir}/${isolate}/${isolate}_pilS_nt_full_blast_sequences ${path_dir}/${isolate}/${isolate}_pilS_nt_full_blast_seq.fasta #the sequence of the aligning alleles is saved to a fasta file that gets uploaded to geneious prime to annotate the genome

