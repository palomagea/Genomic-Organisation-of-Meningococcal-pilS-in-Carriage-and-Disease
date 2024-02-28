#!/bin/bash
#this was the first pipeline used to analyse the sequencing run and make a genome assembly from the fastq files
#this pipeline was before the issues with medaka were discovered

#start in correct folder
echo Enter pwd for the folder you want to be working in 
read path_dir

#check which isolate you are analysing and set variable name 
echo Enter name of run 
read run 

echo Enter isolate name
read isolate
echo Analysing sequence run from isolate $isolate 


#make folder for this isolate to store analysis
mkdir $isolate

#locate path where the sequencing run is stored
echo Enter path of isolate sequencing run
read seq_run_path

#make summary file and save in home directory
cd $seq_run_path
cd ..
cd ..
module load pycoqc/2.5.0.3
pycoQC --summary_file sequencing_summary* -o /home/pgea/${run}_qc.html
 
#cat together fastq files from run and move into new directory
cd $path_dir
cat $seq_run_path/*.fastq.gz > $path_dir/$isolate/$isolate.fastq.gz

#filter the fastq files to remove any with a q score below 12 or a length below 5kb
module load nanopack/1.1.1
cd $isolate
gunzip $isolate.fastq.gz
NanoFilt -q 12 --length 5000 $isolate.fastq > ${isolate}_filt.fastq

#run flye to assemble genome and correction step with medaka and racon 
module load flye/2.8-1
module load minimap2/2.24
module load racon/1.4.3
module load medaka/1.7.1

flye --nano-raw ${isolate}_filt.fastq --threads 20 --iterations 3 --out-dir ${isolate}_flye
minimap2 -x map-ont -t 20 ${isolate}_flye/assembly.fasta ${isolate}_filt.fastq > ${isolate}_flye/assembly_racon.paf
racon -m 8 -x -6 -g -8 -w 500 -t 20 ${isolate}_filt.fastq ${isolate}_flye/assembly_racon.paf ${isolate}_flye/assembly.fasta > ${isolate}_flye/racon_on_flye.fasta
medaka_consensus -i ${isolate}_filt.fastq -d ${isolate}_flye/racon_on_flye.fasta -o ${isolate}_flye/medakon_on_racon -t 3 -m r941_min_hac_g507


#copy polished assembly out of folder and into main directory
cp $path_dir/$isolate/*_flye/medakon_on_racon/consensus.fasta $path_dir/$isolate/${isolate}_consensus.fasta

#give some stats about flye assembly like coverage and N50/N90 
cd $path_dir/$isolate/*_flye
grep 'Mean coverage\|Reads N50\|Total read length\|Total length' flye.log > $path_dir/$isolate/genome_assembly_stats.txt

cd $path_dir/$isolate

#annotate with prokka 
mkdir annotation
cp ${isolate}_consensus.fasta $path_dir/$isolate/annotation
cd annotation

module load prokka/1.14.5
prokka ${isolate}_consensus.fasta --genus neisseria -- prefix $isolate

#make sam bam bai folders index
#use minimap to map the fastq reads onto the flye assembly just built to make a bam file and index it
module load samtools/1.9
module load minimap2/2.24
module load python/3.7.3
module load bedtools2/2.19.1

cd $path_dir/$isolate
minimap2 -ax map-ont $path_dir/$isolate/${isolate}_consensus.fasta $path_dir/$isolate/${isolate}_filt.fastq | samtools view -bS | samtools sort -o ${isolate}_allreads.bam
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


