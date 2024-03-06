#This pipeline will create an alignment file of all the reads mapping to pilS and then count the depth at each pilS position
#To run, this pipeline requires the pilS bed file of the reads mapping across pilS 

#!/bin/bash 

######### Set the variables #########
isolate="the_name_of_the_isolate_or_experiment" #Use the same isolate name as in the seq_analysis.sh pipeline so that all of the files are named consistently

#Set the working directory with the scripts saved in it 
path_dir=$(pwd)

#Load modules
module load python/3.7.3
module load bedtools/2.31.0
module load samtools/1.9

cd ${path_dir}/${isolate}

#Get txt of all the sequence names of sequences spanning pilS from the pilS.bed file and then add to a pilS_read_names.txt file 
python ${path_dir}/extract_sequence_name_from_bed.py ${path_dir}/${isolate}/${isolate}_pilS.bed ${path_dir}/${isolate}/pilS_read_names.txt

#Use the list of sequence names to get a sam file containing the reads that map pilS from the original allreads bam file 
samtools view -H ${isolate}_allreads.bam > header_for_pilS_sam.txt #Need to make a header first for formatting 
samtools view ${isolate}_allreads.bam | grep -f 'pilS_read_names.txt' > ${isolate}_pilS.sam #Search the bam file for reads that are in the pilS_read_names.txt file
cat header_for_pilS_sam.txt ${isolate}_pilS.sam > ${isolate}_pilS_with_header.sam #Add the header to the file

#Convert the sam file to a bam alignment file and index
samtools view -bS ${isolate}_pilS_with_header.sam | samtools sort -o ${isolate}_pilS.bam
samtools index ${isolate}_pilS.bam

#Get depths at each positions in the pilS region in the pilS bam file
#Since the bam file is reads mapping the whole way across pilS, expect the depth to be the same at every position
#When there is a change in depth of coverage, this signals a deletion 
samtools depth -b query.bed ${isolate}_pilS.bam > ${isolate}_pilS_depths.txt #save the position and the depth to a pilS_depths.txt file 
