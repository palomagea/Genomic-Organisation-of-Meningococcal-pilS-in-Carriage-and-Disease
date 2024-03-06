#This script is used on the corrected consensus fasta file, which is the Flye assembly that has been corrected using Unicycler with Illumina reads 
#Running this script assumes that the filtered fastq reads and the corrected consensus assembly for the isolate are saved in the isolate directory and everything else has been moved to a new "not Illumina corrected" sub directory
#The steps in the seq_analysis directory are repeated with the corrected genome assembly

#!/bin/bash

######### Set the variables #########
isolate="the_name_of_the_isolate_or_experiment" #Use the same isolate name as in the seq_analysis.sh pipeline so that all of the files are named consistently

#Set the working directory with the scripts saved in it 
path_dir=$(pwd)


#Load modules 
module load prokka/1.14.5
module load samtools/1.9
module load minimap2/2.24
module load python/3.7.3
module load bedtools2/2.19.1

cd $path_dir/$isolate

#Prokka is used to annotate the consensus genome assembly. This is used to identify the fbp and lpxC genes that flank the pilS region.
mkdir annotation
cp ${isolate}_corrected_consensus.fasta $path_dir/$isolate/annotation
cd annotation
prokka ${isolate}_corrected_consensus.fasta --genus neisseria -- prefix $isolate

#An alignment file is generated of the filtered fastq reads mapped to the consensus assembly 
cd $path_dir/$isolate
minimap2 -ax map-ont $path_dir/$isolate/${isolate}_corrected_consensus.fasta $path_dir/$isolate/${isolate}_filt.fastq | samtools view -bS | samtools sort -o ${isolate}_allreads.bam #Minimap was used to map the fastq reads back onto the assembly
samtools index ${isolate}_allreads.bam 

#Generate a txt file from the alignment file that contains the depth at each genomic position. 
#The txt file has three columns: contig name, position and depth at that positon
samtools depth -aa ${isolate}_allreads.bam > ${isolate}_coverage.txt

#Generate txt files with the lengths of each read mapped to the genome
samtools view -F 4 ${isolate}_allreads.bam | cut -f 10 | perl -ne 'chomp;print length ($_). "\n"' | sort -n | uniq -c > ${isolate}_numberreads_readlength.txt #a txt file is made with two columns, the length of the read (right column) and the number of reads that are that length (left column)
samtools view -F 4 ${isolate}_allreads.bam | awk '{print length($10)}' | sort -n > listreadlengths.txt #a txt file is generated that has the length of every mapped read on a new line

#Run python scripts that will print the average depth and standard deviation and average length of mapped reads 
python ${path_dir}/average_depth_std.py ${path_dir}/${isolate}/${isolate}_coverage.txt
python ${path_dir}/average_read_length.py ${path_dir}/${isolate}/listreadlengths.txt

#Print the length of the longest read
echo the longest read length was 
tail -1 $path_dir/$isolate/listreadlengths.txt
