#This script is used on the corrected consensus fasta file, which is the Flye assembly that has been corrected using Unicycler with Illumina reads 
#Running this script assumes that the filtered fastq reads and the corrected consensus assembly for the isolate are saved in the isolate directory and everything else has been moved to a new "not Illumina corrected" sub directory
#The steps in the seq_analysis directory are repeated with the corrected genome assembly

#!/bin/bash



cd $path_dir/$isolate

#prokka annotation
mkdir annotation
cp ${isolate}_corrected_consensus.fasta $path_dir/$isolate/annotation
cd annotation

module load prokka/1.14.5
prokka ${isolate}_corrected_consensus.fasta --genus neisseria -- prefix $isolate

#make sam bam bai folders index
#map the filtered reads onto the consensus assembly to make a bam file 
module load samtools/1.9
module load minimap2/2.24
module load python/3.7.3
module load bedtools2/2.19.1

cd $path_dir/$isolate
minimap2 -ax map-ont $path_dir/$isolate/${isolate}_corrected_consensus.fasta $path_dir/$isolate/${isolate}_filt.fastq | samtools view -bS | samtools sort -o ${isolate}_allreads.bam
samtools index ${isolate}_allreads.bam 

#depth of coverage all bps txt file ref seq, base index, coverage as headers 
samtools depth -aa ${isolate}_allreads.bam > ${isolate}_coverage.txt

#get txt file with frequency and num of reads mapped to genome
samtools view -F 4 ${isolate}_allreads.bam | cut -f 10 | perl -ne 'chomp;print length ($_). "\n"' | sort -n | uniq -c > ${isolate}_numberreads_readlength.txt

#get average reads mapped to genome 
samtools view -F 4 ${isolate}_allreads.bam | awk '{print length($10)}' | sort -n > listreadlengths.txt 

#run python scripts for analysis need to be in folder with scripts 
cd $path_dir

python average_depth_std.py
python average_read_length.py
echo the longest read length was 
tail -1 $path_dir/$isolate/listreadlengths.txt


