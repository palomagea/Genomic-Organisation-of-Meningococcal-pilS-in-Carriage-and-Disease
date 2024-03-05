#This script breaks the genome into windows, extracts the reads that map across the whole window and calculates the depth at each position. 
#The window name, mean, mode and median coverage, sample_variance and standard_deviation of coverage per window and minimum and maximum covergae are calculated for each window and recorded in the output file "averages per window".
#The sample variance of each window across the genome can be plotted to look at the variance across the genome.
#Windows with high variance signal that there are deletions within that region. 

#This script works with the "make_sliding_windows.py" script to generate a list of the windows for the genome to be broken into
#This can be edited to change the window size and step across the genome.

#!/bin/bash

#Variables to set:
allreads_bam="path_and_name_of_whole_genome_alignment_file" #this is the alignment file of all the nanopore fastq reads aligned to the consensus genome assembly
isolate="name_of_the_isolate or sample" #this will be the name of the folder that all the output goes into

mkdir ${isolate}
cd ${isolate}

#create a WG_sliding widowns query.bed file across the whole genome
module load python/3.7.3
module load bedtools/2.31.0
module load samtools/1.9

python make_sliding_windows.py

#use the WG_sliding windows query file and iterate through each window in the file to pull out all of the reads that map completely across that window
#make separate bam files for each window of only the reads mapping fully across and then get the depth of coverage across each windows bam file
#because the bam files only contain reads mapping the whole way across, regions of low depth of coverage correspond to a deletion

output_file="averages_per_window.txt"
echo -e "window\tmean\tmode\tmedian\tsample_variance\tstandard_deviation\tminimum\tmaximum" > "$output_file"

while IFS=$'\t' read -r chrom start end; do
    window_name="window_${start}_${end}"
    window_name_clean=$(echo "$window_name" | sed 's/[^a-zA-Z0-9]/_/g')
    window_bed="${window_name_clean}.bed"
    window_sam="${window_name_clean}.sam"
    window_bam="${window_name_clean}.bam"
    query_bed="${window_name_clean}_query.bed"  

    echo -e "$chrom\t$start\t$end" > "$query_bed"  # Create the query BED file
    bedtools intersect -a "${path_dir}/${isolate}/${isolate}_allreads.bed" -b "$query_bed" -F 1 > "${window_name_clean}.bed"

    if [ -s "${window_name_clean}.bed" ]; then
        python ${path_dir}/deletions_across_genome/extract_seq_names_from_bed.py ${window_name_clean}.bed ${window_name_clean}_read_names.txt
        samtools view -H "$allreads_bam" > "header_${window_name_clean}.sam"
        samtools view "$allreads_bam" | grep -f "${window_name_clean}_read_names.txt" > "$window_sam"
        cat "header_${window_name_clean}.sam" "$window_sam" > "${window_name_clean}_with_header.sam"
        samtools view -bS "${window_name_clean}_with_header.sam" | samtools sort -o "$window_bam"
        samtools index "$window_bam"
        samtools depth -b "$window_bed" "$window_bam" > "${window_name_clean}_depths.txt"
    
#Calculate averages for each window and print to averages per window txt file 
        values=()
        while read -r _ _ depth; do
            values+=("$depth")
        done < "${window_name_clean}_depths.txt"

        mean=$(printf '%s\n' "${values[@]}" | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }')
        mode=$(printf '%s\n' "${values[@]}" | sort | uniq -c | sort -rn | head -1 | awk '{ print $2 }')
   
        n=${#values[@]}
        sorted_values=($(printf '%s\n' "${values[@]}" | sort -n))

        median_index=$((n / 2))
        if ((n % 2 == 0)); then
            median=$(( (sorted_values[median_index - 1] + sorted_values[median_index]) / 2 ))
        else
            median="${sorted_values[median_index]}"
        fi

        sum=0
        for val in "${values[@]}"; do
            sum=$((sum + val))
        done
        avg=$((sum / ${#values[@]}))

        variance=0
        for val in "${values[@]}"; do
            diff=$((val - avg))
            variance=$((variance + (diff * diff)))
        done
        sample_variance=$(awk "BEGIN {print $variance / (${#values[@]} - 1)}")
        standard_deviation=$(awk "BEGIN {print sqrt($variance / ${#values[@]})}")

        min=$(printf '%s\n' "${values[@]}" | sort -n | head -1)
        max=$(printf '%s\n' "${values[@]}" | sort -n | tail -1)

        echo -e "$window_name\t$mean\t$mode\t$median\t$sample_variance\t$standard_deviation\t$min\t$max" >> "$output_file"

        rm "$window_bed" "${window_name_clean}_read_names.txt" "header_${window_name_clean}.sam"
        
    else
        echo "Skipping empty window: $window_name"
    fi

done < "${isolate}_WG_windows_query.bed"

