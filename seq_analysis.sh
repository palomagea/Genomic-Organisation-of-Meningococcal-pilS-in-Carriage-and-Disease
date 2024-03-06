#This is the first pipeline used to analyse the sequencing run and make a genome assembly from the fastq files, annotate with prokka and basic assembly statistics, average depth, read length and longest read
#This pipeline was used before the issues with medaka were discovered

#!/bin/bash

######### Set the variables #########
run="set_the_name_of_the_run" #This name is used to name the pycoQC summary file for this run and distinguish it from the other sequencing runs
isolate="the_name_of_the_isolate_or_experiment" #Multiple isolates were sequenced in the same run with different barcodes. This attaches the isolate name to the sequences for that barcode
seq_run_path="the_name_and_path_where_the_fastq_runs_are_saved" #This is the path to folder of fastq files with the barcode corresponding to this isolate
seq_summary="the_name_and_path_of_the_sequencing_summary_file_for_this_run" #The path to the sequencing summary file for this run and the name of the summary file

#Load modules 
module load pycoqc/2.5.0.3
module load nanopack/1.1.1
module load flye/2.8-1
module load minimap2/2.24
module load racon/1.4.3
module load medaka/1.7.1
module load prokka/1.14.5
module load samtools/1.9
module load python/3.7.3
module load bedtools2/2.19.1

#Make a summary file of the sequencing run
pycoQC --summary_file seq_summary -o ${run}_qc.html

#Make a directory for the isolate 
mkdir ${isolate} 

#Join all of the fastq files for the isolate into a single fastq and move into the isolate directory
cat $seq_run_path/*.fastq.gz > $path_dir/$isolate/$isolate.fastq.gz

#Filter the fastq files to remove any with a q score below 12 or a length below 5kb
cd ${path_dir}/$isolate
gunzip $isolate.fastq.gz
NanoFilt -q 12 --length 5000 $isolate.fastq > ${isolate}_filt.fastq

#Use to assemble genome and polish with Medaka and Racon 
flye --nano-raw ${isolate}_filt.fastq --threads 20 --iterations 3 --out-dir ${isolate}_flye #Flye is run on the filtered fastq reads and saved in an output directory
minimap2 -x map-ont -t 20 ${isolate}_flye/assembly.fasta ${isolate}_filt.fastq > ${isolate}_flye/assembly_racon.paf #Minimap is used to map the filtered fastq onto the flye assembled genome and saved in paf format to use with Racon
racon -m 8 -x -6 -g -8 -w 500 -t 20 ${isolate}_filt.fastq ${isolate}_flye/assembly_racon.paf ${isolate}_flye/assembly.fasta > ${isolate}_flye/racon_on_flye.fasta #Racon is used to polish the flye assembly
medaka_consensus -i ${isolate}_filt.fastq -d ${isolate}_flye/racon_on_flye.fasta -o ${isolate}_flye/medakon_on_racon -t 3 -m r941_min_hac_g507 #Medaka is used to polish the Racon on Flye assembly

#Copy polished assembly out of folder and into main directory
cp $path_dir/$isolate/*_flye/medakon_on_racon/consensus.fasta $path_dir/$isolate/${isolate}_consensus.fasta

#Give some stats about flye assembly like coverage and N50/N90 
cd $path_dir/$isolate/*_flye 
grep 'Mean coverage\|Reads N50\|Total read length\|Total length' flye.log > $path_dir/$isolate/genome_assembly_stats.txt #the Flye log is searched for the mean coverage, N50, total read length and genome length and this is recorded in a txt file
cd $path_dir/$isolate

#Prokka is used to annotate the consensus genome assembly. This is used to identify the fbp and lpxC genes that flank the pilS region.
mkdir annotation
cp ${isolate}_consensus.fasta $path_dir/$isolate/annotation #the consensus genome is added to the annotation directory to run Prokka in
cd annotation
prokka ${isolate}_consensus.fasta --genus neisseria -- prefix $isolate

#An alignment file is generated of the filtered fastq reads mapped to the consensus assembly 
cd $path_dir/$isolate
minimap2 -ax map-ont $path_dir/$isolate/${isolate}_consensus.fasta $path_dir/$isolate/${isolate}_filt.fastq | samtools view -bS | samtools sort -o ${isolate}_allreads.bam #Minimap was used to map the fastq reads back onto the assembly
samtools index ${isolate}_allreads.bam 

#depth of coverage all bps txt file ref seq, base index, coverage as headers 
samtools depth -aa ${isolate}_allreads.bam > ${isolate}_coverage.txt

#get txt file with the numbers of reads of each length mapped to genome
samtools view -F 4 ${isolate}_allreads.bam | cut -f 10 | perl -ne 'chomp;print length ($_). "\n"' | sort -n | uniq -c > ${isolate}_numberreads_readlength.txt

#get average reads mapped to genome 
samtools view -F 4 ${isolate}_allreads.bam | awk '{print length($10)}' | sort -n > listreadlengths.txt 

#run python scripts for analysis need to be in folder with scripts 
#these scripts give the average depth of coverage and standard deviation, the average read length and the longest read length 
cd $path_dir

python average_depth_std.py
python average_read_length.py
echo the longest read length was 
tail -1 $path_dir/$isolate/listreadlengths.txt


