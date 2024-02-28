#!/bin/bash 
#this pipeline will align RNA seq reads with the ONT reference genomes

#set the path and the name of the RNA seq stranded data 
path_dir="/NGS/active/IPL/MENINGO/analysis/paloma/RNA_seq"

echo Enter isolate name
read isolate 

echo Enter strand one data
read strand_1

echo Enter strand two data 
read strand_2

#make folder for rna seq of isolate 
cd ${path_dir}
mkdir $isolate
cd ${isolate}

#copy reference genome to folder and write index file for reference genome
module load bwa/0.7.17
cp /NGS/active/IPL/MENINGO/analysis/paloma/2023/${isolate}/${isolate}_corrected_consensus.fasta ${path_dir}/${isolate}
bwa index -a is ${isolate}_corrected_consensus.fasta

#align fastq files 
module load samtools/1.9
bwa mem -t 4 ${isolate}_corrected_consensus.fasta /NGS/active/IPL/MENINGO/analysis/transcriptome/${strand_1} /NGS/active/IPL/MENINGO/analysis/transcriptome/${strand_2} > ${isolate}_RNA_seq_aligned_to_reference.sam 
samtools view -S -b ${isolate}_RNA_seq_aligned_to_reference.sam > ${isolate}_RNA_seq_aligned_to_reference.bam

#sort bam file 
samtools sort ${isolate}_RNA_seq_aligned_to_reference.bam > ${isolate}_RNA_seq_aligned_to_reference.sorted.bam

#clean bam file
#remove unmapped reads, include only reads that are properly paired 
# -F 1024 removes PCR or optical duplicates 
#remove reads with secondary alignments
#remove aligned reads that have failed instrument QC
#add unas filter to only take reads where the 150 paired reads are 300-700bp 
samtools view -f 3 -F 4 -F 1024 -b -o ${isolate}_RNA_seq_aligned_to_reference.sorted.temp1.bam ${isolate}_RNA_seq_aligned_to_reference.sorted.bam
samtools view -b -F 0x100  ${isolate}_RNA_seq_aligned_to_reference.sorted.temp1.bam >${isolate}_RNA_seq_aligned_to_reference.sorted.temp2.bam
samtools rmdup ${isolate}_RNA_seq_aligned_to_reference.sorted.temp2.bam ${isolate}_RNA_seq_aligned_to_reference.temp3.sorted.bam
samtools view -h ${isolate}_RNA_seq_aligned_to_reference.temp3.sorted.bam | awk 'substr($0,1,1)=="@" || ($9>= 300 && $9<=700) || ($9<=-300 && $9>=-700)' | samtools view -b > ${isolate}_RNA_seq_aligned_to_reference.clean.sorted.bam
rm ${isolate}_RNA_seq_aligned_to_reference.sorted.temp1.bam ${isolate}_RNA_seq_aligned_to_reference.sorted.temp2.bam ${isolate}_RNA_seq_aligned_to_reference.temp3.sorted.bam

#index bam file 
samtools index ${isolate}_RNA_seq_aligned_to_reference.clean.sorted.bam



