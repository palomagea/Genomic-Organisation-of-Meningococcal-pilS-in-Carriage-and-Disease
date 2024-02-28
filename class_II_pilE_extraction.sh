#!/bin/bash
#this pipeline will extract class II pilE it is different because class II pilE are not located inside the pilS region

#start in correct folder
echo Enter pwd for the folder you want to be working in
read path_dir

#check which isolate you are analysing and set variable name

echo Enter isolate name
read isolate

module load python/3.7.3
module load bedtools2/2.19.1
module load samtools/1.9

#make query.bed for the pilE region from input 
python make_pilE_query.py

#create genome file from input 
python create_genome.py

cd $path_dir/$isolate
bedtools intersect -a ${isolate}_allreads.bed -b pilE_query.bed -F 1 > ${isolate}_class_II_pilE.bed 
bedToBam -i ${isolate}_class_II_pilE.bed -g $path_dir/$isolate > ${isolate}_class_II_pilE.bam 
samtools index ${isolate}_class_II_pilE.bam

cd $path_dir


#make pilE fasta
cd $isolate 
bedtools getfasta -fi $path_dir/$isolate/${isolate}_corrected_consensus.fasta -bed pilE_query.bed -fo ${isolate}_class_II_pilE.fasta 

#do bakta run on pilS fasta
cd ${path_dir}/${isolate}
cd bakta
module load bakta/1.5.0
bakta --output ${isolate}_class_II_pilE --genus Neisseria ${path_dir}/${isolate}/${isolate}_class_II_pilE.fasta 

cd ${path_dir}/${isolate}

#full blast to get all sequences that align  
module load ncbi-blast/2.12.0
blastn -db /NGS/active/IPL/MENINGO/analysis/paloma/pilS_nt_db/pilS_nt_db -query ${isolate}_class_II_pilE.fasta -outfmt "6 sseqid qstart qend sseq" -out ${isolate}_class_II_pilE_nt_full_blast_sequences -max_target_seqs 6000

#run python script to get the sequences from the full blast   
python ${path_dir}/extract_pilE_full_blast_seq.py