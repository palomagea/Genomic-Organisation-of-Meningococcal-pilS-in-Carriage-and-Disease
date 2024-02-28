
#!/bin/bash
#this pipeline will extract the pilS region from the consensus assembly. need to know the boundaries of where pilS is from the prokka annotation of fbp and lpxC and then add 1kb to each side
#this pipeline will then annotate pilS with bakta and make a bam file of the reads that map the whole way across pilS and count the reads that map the whole way across pilS
#it will also get the blast results from the database that have sequence similarity to pilS

#start in correct folder
echo Enter pwd for the folder you want to be working in
read path_dir

#check which isolate you are analysing and set variable name

echo Enter isolate name
read isolate
echo Analysing sequence run from isolate $isolate

module load python/3.7.3
module load bedtools/2.31.0
module load samtools/1.9

#make query.bed for the pilS region from input 
python make_query.py

#get length of genome and save to txt file 
awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' $path_dir/$isolate/${isolate}_corrected_consensus.fasta > $path_dir/$isolate/genome_length.txt 

#create genome file from input 
python create_genome.py

#extract only the reads that map the whole way across the pilS region
cd $path_dir/$isolate
bedtools intersect -a ${isolate}_allreads.bed -b query.bed -F 1 > ${isolate}_pilS.bed 

cd $path_dir

#count reads mapping pilS
python count_pilS_reads.py

#make pilS fasta
cd $isolate 
bedtools getfasta -fi $path_dir/$isolate/${isolate}_corrected_consensus.fasta -bed query.bed -fo ${isolate}_pilS.fasta 

#do bakta run on pilS fasta
cd ${path_dir}/${isolate}
mkdir bakta
cd bakta
module load bakta/1.5.0
bakta --output ${isolate}_pilS --genus Neisseria ${path_dir}/${isolate}/${isolate}_pilS.fasta 

#full blast to get all sequences that align  
blastn -db /NGS/active/IPL/MENINGO/analysis/paloma/pilS_nt_db/pilS_nt_db -query ${isolate}_pilS.fasta -outfmt "6 sseqid qstart qend sseq" -out ${isolate}_pilS_nt_full_blast_sequences -max_target_seqs 6000

#run python script to get the sequences from the full blast   
python ${path_dir}/extract_full_blast_seq.py

