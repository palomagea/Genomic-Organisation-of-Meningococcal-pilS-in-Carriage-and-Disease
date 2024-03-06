#This pipeline will align the Illumina RNA seq reads to the polished genomes. 
#To run, this pipeline required paired end RNA seq reads

#!/bin/bash 

######### Set the variables #########
isolate="the_name_of_the_isolate_or_experiment" Use the same isolate name as in the seq_analysis.sh pipeline so that all of the files are named consistently
pair_1="the_name_and_path_of_RNA_seq_reads_pair_1" #Path to and name of the RN Aseq reads pair 1
pair_2="the_name_and_path_of_RNA_seq_reads_pair_2" #Path to and name of the RNA seq reads pair 2

#Set the working directory with the scripts saved in it 
path_dir=$(pwd)

#Load modules
module load bwa/0.7.17
module load samtools/1.9


#Make folder for RNA seq Data 
cd ${isolate}
mkdir RNA_seq
cd RNA_seq
cp ${path_dir}/${isolate}/${isolate}_corrected_consensus.fasta ${path_dir}/${isolate}/RNA_seq #copy the genome assembly file into the RNA seq folder and index it
bwa index -a is ${isolate}_corrected_consensus.fasta

#Align the RNA seq fastq files to the genome assembly using BWA 
bwa mem -t 4 ${isolate}_corrected_consensus.fasta /NGS/active/IPL/MENINGO/analysis/transcriptome/${pair_1} /NGS/active/IPL/MENINGO/analysis/transcriptome/${pair} > ${isolate}_RNA_seq_aligned_to_reference.sam #both paired ends are aligned to the genome assembly
samtools view -S -b ${isolate}_RNA_seq_aligned_to_reference.sam > ${isolate}_RNA_seq_aligned_to_reference.bam #sam file is changed to a bam file
samtools sort ${isolate}_RNA_seq_aligned_to_reference.bam > ${isolate}_RNA_seq_aligned_to_reference.sorted.bam #bam file is sorted 

#Filter and clean the alignment file
#Remove unmapped reads, include only reads that are properly paired 
# -F 1024 removes PCR or optical duplicates 
#Remove reads with secondary alignments
#Remove aligned reads that have failed instrument QC
#Add paired end filter to only retain reads with an insert size of 300-700bp to reflect the size of the pilE gene
samtools view -f 3 -F 4 -F 1024 -b -o ${isolate}_RNA_seq_aligned_to_reference.sorted.temp1.bam ${isolate}_RNA_seq_aligned_to_reference.sorted.bam
samtools view -b -F 0x100  ${isolate}_RNA_seq_aligned_to_reference.sorted.temp1.bam >${isolate}_RNA_seq_aligned_to_reference.sorted.temp2.bam
samtools rmdup ${isolate}_RNA_seq_aligned_to_reference.sorted.temp2.bam ${isolate}_RNA_seq_aligned_to_reference.temp3.sorted.bam
samtools view -h ${isolate}_RNA_seq_aligned_to_reference.temp3.sorted.bam | awk 'substr($0,1,1)=="@" || ($9>= 300 && $9<=700) || ($9<=-300 && $9>=-700)' | samtools view -b > ${isolate}_RNA_seq_aligned_to_reference.clean.sorted.bam
rm ${isolate}_RNA_seq_aligned_to_reference.sorted.temp1.bam ${isolate}_RNA_seq_aligned_to_reference.sorted.temp2.bam ${isolate}_RNA_seq_aligned_to_reference.temp3.sorted.bam #Remove temporary files

#Index the final alignment file 
samtools index ${isolate}_RNA_seq_aligned_to_reference.clean.sorted.bam



