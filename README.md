# Genomic Organisation of Meningococcal pilS in Carriage and Disease

**Author:** Paloma Gea

These are the scripts used in the thesis Genomic Organisation of Meningococcal pilS in Carriage and Disease

#### Requirements
\- Oxford Nanopore sequence data for meningococcal isolates including the fastq sequences for each barcode and the sequencing summary file \
\- Illumina sequencing data for each of the isolates\
\- Illumina RNA sequencing reads

###### To complete the analysis outlined in the thesis, the following steps have been carried out:

#### 1. Run the seq_analysis.sh pipeline
Set the following variables within the pipeline before use
- run="set_the_name_of_the_run" This name is used to name the pycoQC summary file for this run and distinguish it from the other sequencing runs
- isolate="the_name_of_the_isolate_or_experiment" Multiple isolates were sequenced in the same run with different barcodes. This attaches the isolate name to the sequences for that barcode
- seq_run_path="the_name_and_path_where_the_fastq_runs_are_saved" This is the path to folder of fastq files with the barcode corresponding to this isolate
- seq_summary="the_name_and_path_of_the_sequencing_summary_file_for_this_run" The path to the sequencing summary file for this run and the name of the summary file

This pipeline will
- Generate a pycoQC sequencing summary report
- Filter the fastq reads and remove any reads below q score 12 and less than 5kb in length
- Assemble the reads into a genome using Flye
- Polish the genome using Racon and Medaka
- Annotate the genome using Prokka 
- Create an alignment of fastq reads mapped onto the assembly genome
- Print the average depth across the genome and the standard deviation
- Print the average length of mapped reads
- Print the length of the longest mapped read

Despite having to repeat many of these steps after correcting with Illumina, this pipeline has been left to show how the genomes were polished with Medaka and Racon and what the flye assemblies looked like before and after Illumina correction. 

#### 2. Illumina Corrections
As outlined in the thesis, there were issues encountered when polishing with Medaka. The filtered fastq reads were retained in the directory for that isolat,e and the rest of the output from the first pipeline was into a new sub-directory "not_Illumina_corrected". This was later used to compare the Illumina corrected sequence with the not Illumina corrected sequence.
The following steps were carried out to use the Illumina reads to correct the flye assembly.

```bash
#unicycler Version 0.4.8
unicycler_polish --threads 20 --ale [path_to_ALE_executable_file] -1 [path_to_illumina_reads_first_reads_in_pair] -2 [path_to_illumina_reads_second_reads_in_pair] -a [path_to_the_flye_genome_assembly]
```

The final genome generated by Unicycler was saved as "isolate_corrected_consensus.fasta" where isolate is the same isolate name set in the seq_analysis.sh pipeline
The corrected genome assembly is moved into the directory for the isolate and is used for all subsequent analysis. The other files generated by Unicycler were deleted.

#### 3. Run the after_correction.sh pipeline
This pipeline repeats the same steps as the seq_analysis.sh pipeline but uses the genome assembly file that has been corrected with Illumina reads.
Set the following variables within the pipeline before use 
- isolate="the_name_of_the_isolate_or_experiment" Use the same isolate name as in the seq_analysis.sh pipeline so that all of the files are named consistently

This pipeline will
- Annotate the genome using Prokka 
- Create an alignment of fastq reads mapped onto the assembly genome
- Print the average depth across the genome and the standard deviation
- Print the average length of mapped reads
- Print the length of the longest mapped read

#### 4. Define the pilS region
As outlined in the thesis, the pilS region is flanked by fbp and lpxC genes. The prokka annotation file was used to identify the genomic location of these genes.
The region between and including fbp and lpxC, with an additional 1kb buffer on each side was defined as the pilS region. This was recorded to enter into the next pipeline.

#### 5. BLAST database
All alleles for fbp, lpxC, pilE and pilS were downloaded from pubMLST and joined into one fasta file called pilS_nt_db.fas with a unique identifier for each allele
These alleles were used to make a BLAST database
The following steps were carried out to make the database

```bash
#ncbi-blast Version 2.6.0
makeblastdb -in [pilS_nt_db.fas] -dbtype nucl -input_type fasta -out [pilS_nt_db]
```

The pilS nucleotide database generated was saved and is used in the extract_pilS_annotation.sh pipeline below
The database used in this thesis has been uploaded into the repository.

#### 6. Run the extract_pilS_annotation.sh pipeline
For this pipeline to run need to know the boundaries of the pilS region and have created the BLAST database from the PubMLST alleles and saved it as pilS_nt_db
Set the following variables within the pipeline before use
- isolate="the_name_of_the_isolate_or_experiment" Use the same isolate name as in the seq_analysis.sh pipeline so that all of the files are named consistently
- pilS_start="the_genomic_start_positon_of_pilS" The positon in bp of the start of pilS (including fbp and lpxC +1kb on either end as a buffer)
- pilS_end="the_genomic_end_position_of_pilS" The positon in bp of the end of pilS (including fbp and lpxC +1kb on either end as a buffer)
- pilS_contig="the_contig_that_pilS_is_on" Very few of the isolates had more than one contig. If there was only one contig set "contig_1". In isolates with more than one contig set the contig pilS was on.

This pipeline will 
- Extract all of the reads that map the whole way across pilS and count them
- Extract the fasta sequence of the pilS region
- Annotate the pilS region with Bakta
- Perform a BLAST alignment of the PubMLST alleles in the database to the pilS region
- Save all of the matched alleles into a fasta file

#### 7. Run the class_II_pilE_extraction.sh pipeline
This pipeline only needs to be run for isolates with a class II pilE (NZ97/052, NZCM111 and NZCM112)
For this pipeline to run need to know the genomic location of the pilE gene in the class II isolates, located using the flanking katA and prlC genes and then include a buffer of +1kb on either side
Set the following variables within the pipeline before use
- isolate="the_name_of_the_isolate_or_experiment" Use the same isolate name as in the seq_analysis.sh pipeline so that all of the files are named consistently
- pilE_start="the_genomic_start_positon_of_pilS" The positon in bp of the start of pilS (including fbp and lpxC +1kb on either end as a buffer)
- pilE_end="the_genomic_end_position_of_pilS" The positon in bp of the end of pilS (including fbp and lpxC +1kb on either end as a buffer)
- pilE_contig="the_contig_that_pilS_is_on" Very few of the isolates had more than one contig. If there was only one contig set "contig_1". In isolates with more than one contig set the contig pilS was on.

This pipeline will
- Extract the reads that map the whole way across pilE
- Extract the fasta sequence of pilE
- Annotate the pilE region with Bakta
- Perform a BLAST alignment of the PubMLST alleles in the database to the pilS region
- Save all of the matched alleles into a fasta file

#### 8. Annotate in Geneious
The pilS regions of all isolates and the pilE region of the class II isolates were annotated using Geneious prime with the methods outlined in the thesis.
All of the files needed for this annotation were generated in the pipelines above

#### 9. Run the align_RNA_seq.sh pipeline
This pipeline only needs to be run for the isolates with RNA seq data
Set the following variables within the pipeline before use
- isolate="the_name_of_the_isolate_or_experiment" Use the same isolate name as in the seq_analysis.sh pipeline so that all of the files are named consistently
- pair_1="the_name_and_path_of_RNA_seq_reads_pair_1" 
- pair_2="the_name_and_path_of_RNA_seq_reads_pair_2" 

This pipeline will
- Align the RNA seq reads to the genome assembly using BWA
- Filter the RNA seq reads in the alignment file

#### 10. FeatureCounts
The consensus annotations of pilS were downloaded from geneious prime and reformatted into a gtf tab delimited file
An example of the file is:
```plaintext
##gff-version 2
##sequence-region contig_1 1 2209937
contig_1 Geneious gene 291821 292150 . + . gene_id "fbp"
contig_1 Geneious gene 292175 292544 . + . gene_id "pilS_6"
contig_1 Geneious gene 292590 292949 . + . gene_id "pilS_5"
contig_1 Geneious gene 293175 293560 . + . gene_id "pilS_4"
contig_1 Geneious gene 293624 293983 . + . gene_id "pilS_3"
contig_1 Geneious gene 294029 294394 . + . gene_id "pilS_2"
contig_1 Geneious gene 294729 295149 . + . gene_id "pilS_1"
contig_1 Geneious gene 295965 296498 . + . gene_id "pilE"
contig_1 Geneious gene 297674 298597 . + . gene_id "lpxC"
```
FeatureCounts was then used to count the number of reads mapping to the genomic features 

```bash
# featureCounts Version 2.0.1
#For non stranded
featureCounts -p -t gene -a [path_to_gtf_file.gtf] -o [isolate_RNA_counts.tx]t [isolate_RNA_seq_aligned_to_reference.clean.sorted.bam]

#For stranded 
featureCounts -p -s 2 -t gene -a [path_to_gtf_file.gtf] -o [isolate_RNA_counts.txt] [isolate_RNA_seq_aligned_to_reference.clean.sorted.bam]
```

#### 11. Run the pilS_bed.sh pipeline
Set the following variables within the pipeline before use
- isolate="the_name_of_the_isolate_or_experiment" Use the same isolate name as in the seq_analysis.sh pipeline so that all of the files are named consistently

This pipeline will
- Convert the pilS bed file of reads mapping across pilS to a bam alignment file
- Calculate the depth of coverage at each position in pilS to look for deletions in the pilS region

#### 12. Run the deletions_across_genome.sh pipeline
Set the following variables within the pipeline before use
- isolate="name_of_the_isolate_or_sample" Use the same isolate name as in the seq_analysis.sh pipeline so that all of the files are named consistently
- WG_fasta="name_and_path_of_the_whole_genome_fasta_file" The genome assembly for the specified isolate
- allreads_bam="path_and_name_of_whole_genome_alignment_file" This is the alignment file of all the nanopore fastq reads aligned to the consensus genome assembly
- window_size="set_the_window_size" The size of the window. In this thesis I used 10000 however this can be changed
- window_step="set_the_window_step" The step between windows. In this thesis I used 500 however this can be changed

This pipeline will
- Make sliding windows across the genome of the set size and step
- Extract all of the reads mapping across the whole window and make an alignment file
- Count the depth of coverage at each position in the window
- Calculate the variance in depth at each window, the signal for a deletion is high variance in depth
- Save the averages for each window into a txt file that can be graphed to look at which genomic windows have high variance




