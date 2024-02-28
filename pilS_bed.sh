#!/bin/bash
#this pipeline will create an alignment file of all the reads mapping to pilS and then count the depth at each pilS position

#start in correct folder
path_dir="/NGS/active/IPL/MENINGO/analysis/paloma/2023"

#check which isolate you are analysing and set variable name

echo Enter isolate name
read isolate

module load python/3.7.3
module load bedtools/2.31.0
module load samtools/1.9

python create_genome.py

cd ${path_dir}/${isolate}

#get sequences that span whole pilS region
bedtools intersect -a ${isolate}_allreads.bed -b query.bed -F 1 > ${isolate}_pilS.bed

#get txt of all the sequence names of sequences spanning pilS into a file called pilS_read_names.txt
python ${path_dir}/extract_sequence_name_from_bed.py

#use the list of sequence names to get a sam file containing the pilS mapping reads from the original allreads bam file 
samtools view -H ${isolate}_allreads.bam > header_for_pilS_sam.txt
samtools view ${isolate}_allreads.bam | grep -f 'pilS_read_names.txt' > ${isolate}_pilS.sam
cat header_for_pilS_sam.txt ${isolate}_pilS.sam > ${isolate}_pilS_with_header.sam

#chnage sam file to bam file and index 
samtools view -bS ${isolate}_pilS_with_header.sam | samtools sort -o ${isolate}_pilS.bam
samtools index ${isolate}_pilS.bam

#get depths at all the positions in the pilS region in the pilS bam file of only uninterrupted reads mapping pilS
samtools depth -b query.bed ${isolate}_pilS.bam > ${isolate}_pilS_depths.txt
