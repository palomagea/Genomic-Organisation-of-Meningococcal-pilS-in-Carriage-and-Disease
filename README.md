These scripts have been used in the thesis "Genomic Organisation of Meningococcal pilS in Carriage and Disease" by Paloma Gea 

**The following is required to run this analysis**
- Oxford Nanopore sequence data for meningococcal isolates including the fastq sequences for each barcode in the run that corresponds to the different isolates and the sequencing summary file
- Illumina sequencing data for each of the isolates

To complete the analysis outlined in the thesis, the following steps have been carried out:

1. Run the seq_analysis.sh pipeline
Set the following variables within the pipeline before use
run="set_the_name_of_the_run" #This name is used to name the pycoQC summary file for this run and distinguish it from the other sequencing runs
isolate="the_name_of_the_isolate_or_experiment" #Multiple isolates were sequenced in the same run with different barcodes. This attaches the isolate name to the sequences for that barcode
seq_run_path="the_name_and_path_where_the_fastq_runs_are_saved" #This is the path to folder of fastq files with the barcode corresponding to this isolate
seq_summary="the_name_and_path_of_the_sequencing_summary_file_for_this_run" #The path to the sequencing summary file for this run and the name of the summary file

This pipeline will
- Generate a pycoQC sequencing summary report
- Filter the fastq reads and remove any reads below q score 12 and less than 5kb in length
- Assemble the reads into a genome using Flye
- Polish the genome using Racon and Medaka
- Create an alignment of fastq reads mapped onto the assembly genome
- Print the average depth across the genome and the standard deviation
- Print the average length of mapped reads
- Print the length of the longest mapped read

Despite having to repeat many of these steps after correcting with Illumina, this pipeline has been left to show how the genomes were polished with Medaka and Racon and what the flye assemblies looked like before and after Illumina correction. 

2. Illumina Corrections
As outlined in the thesis, there were issues encountered when polishing with Medaka. The filtered fastq reads were retained in the directory for that isolat,e and the rest of the output from the first pipeline was into a new sub-directory "not_Illumina_corrected". This was later used to compare the Illumina corrected sequence with the not Illumina corrected sequence.
The following steps were carried out to use the Illumina reads to correct the flye assembly.

module load unicycler/0.4.8
unicycler_polish --threads 20 --ale path_to_ALE_executable_file -1 path_to_illumina_reads_first_reads_in_pair -2 path_to_illumina_reads_second_reads_in_pair -a path_to_the_flye_genome_assembly

The final genome generated by Unicycler was saved as "isolate_corrected_consensus.fasta" and moved into the directory for the isolate. The other files were deleted.

3. Run the after_correction.sh pipeline
