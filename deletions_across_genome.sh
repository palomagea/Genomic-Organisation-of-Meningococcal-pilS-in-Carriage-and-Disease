#This script breaks the genome into windows, extracts the reads that map across the whole window and calculates the depth at each position. 
#The window name, and averages including variance are calculated for each window and recorded in the output file "averages per window".
#The sample variance of each window across the genome can be plotted to look at the variance across the genome.
#Windows with high variance signal that there are deletions within that region.
#Need to set the variables isolate name, whole genome fasta, whole genome alignment file

#This script works with the "make_sliding_windows.py" script to generate a list of the windows for the genome to be broken into
#This can be edited to change the window size and step across the genome.

#!/bin/bash

####################Variables to set: ######################
isolate="name_of_the_isolate_or_sample" #Use the same isolate name as in the seq_analysis.sh pipeline so that all of the files are named consistently
WG_fasta="name_and_path_of_the_whole_genome_fasta_file" #The genome assembly for the specified isolate
allreads_bam="path_and_name_of_whole_genome_alignment_file" #This is the alignment file of all the nanopore fastq reads aligned to the consensus genome assembly
window_size="set_the_window_size" #The size of the window. In this thesis I used 10000 however this can be changed
window_step="set_the_window_step" #The step between windows. In this thesis I used 500 however this can be changed

#Set the working directory with the scripts saved in it 
path_dir=$(pwd)

#Load modules
module load python/3.7.3
module load bedtools/2.31.0
module load samtools/1.9

#Make directory for everything to save in
mkdir ${isolate}
cd ${isolate}

#Use python script to create a WG_sliding windows query.bed file of the windows that will span across the whole genome, this can be edited to change window step and size
${path_dir}/python make_sliding_windows.py ${isolate} ${WG_fasta} ${path_dir}/${isolate}/${isolate}_WG_windows_query.bed ${window_size} ${window_step}

#Use the WG_sliding windows query file and iterate through each window defined in the query file 
#For each window pull out all of the reads that map completely across that window and make separate bam alignment files for each window with the reads mapping fully across
#Get the depth of coverage across each windows bam file and calculate the mean, mode and median coverage, sample_variance and standard_deviation of coverage and minimum and maximum covergae per window
#Because the bam files only contain reads mapping the whole way across, regions of low depth of coverage correspond to a deletion

#Create the output file for averages per window to be saved in
output_file="averages_per_window.txt"
echo -e "window\tmean\tmode\tmedian\tsample_variance\tstandard_deviation\tminimum\tmaximum" > "$output_file"

#The following steps are repeated for each window defined in the WG_query bed file
while IFS=$'\t' read -r chrom start end; do
    window_name="window_${start}_${end}" #Name the window based on the contig it came from and the start and end positon
    window_name_clean=$(echo "$window_name" | sed 's/[^a-zA-Z0-9]/_/g') #Remove extra characters from the window name to clean up
    window_bed="${window_name_clean}.bed" #Define the bed file name for each window
    window_sam="${window_name_clean}.sam" #Define the sam file name for each window
    window_bam="${window_name_clean}.bam" #Define the bam file name for each window
    query_bed="${window_name_clean}_query.bed" #Define the query file name for each window

    echo -e "$chrom\t$start\t$end" > "$query_bed"  #Create the query bed file for each window based on what was listed in the WG_query file
    bedtools intersect -a "${path_dir}/${isolate}/${isolate}_allreads.bed" -b "$query_bed" -F 1 > "${window_name_clean}.bed" #Search the WG alignment file for the reads that map only across the specified window 

    if [ -s "${window_name_clean}.bed" ]; then #This if statement makes sure there is content in the bed file of reads mapping across the window by checking if the file size is greater than 0 
        python ${path_dir}/deletions_across_genome/extract_seq_names_from_bed.py ${window_name_clean}.bed ${window_name_clean}_read_names.txt #This python script extracts the sequence names from the whole genome alignment file that map across the specified window and adds them to a new file
        samtools view -H "$allreads_bam" > "header_${window_name_clean}.sam" #This takes the header from the whole genome alignment file and saves it so it can be added to the window sam file later
        samtools view "$allreads_bam" | grep -f "${window_name_clean}_read_names.txt" > "$window_sam" #The reads in the whole genome alignment file with sequence names that map across the specified window are added to the window sam file
        cat "header_${window_name_clean}.sam" "$window_sam" > "${window_name_clean}_with_header.sam" #The header is added to the sam file for the specified window so the format is correct 
        samtools view -bS "${window_name_clean}_with_header.sam" | samtools sort -o "$window_bam" #The sam file is converted to a bam alignment file for the specified window
        samtools index "$window_bam" #The bam file for the specified window is indexed 
        samtools depth -b "$window_bed" "$window_bam" > "${window_name_clean}_depths.txt" #The depth at each position in the window is calculated from the bam alignment file for the window and saved into a text file 
    
#Calculate averages for each window based on the depth at each position and print to averages_per_window.txt file 
#The depth at each position is added to an array named values 
        values=()
        while read -r _ _ depth; do
            values+=("$depth")
        done < "${window_name_clean}_depths.txt"

        mean=$(printf '%s\n' "${values[@]}" | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }') #The mean depth for the window is calculated 
        mode=$(printf '%s\n' "${values[@]}" | sort | uniq -c | sort -rn | head -1 | awk '{ print $2 }') #The mode depth for the window is calculated 
   
        n=${#values[@]} 
        sorted_values=($(printf '%s\n' "${values[@]}" | sort -n)) #The depths get sorted and put in order, saved in sorted_values

        median_index=$((n / 2)) #The median depth for the window is calculated
        if ((n % 2 == 0)); then
            median=$(( (sorted_values[median_index - 1] + sorted_values[median_index]) / 2 )) #If the total number is even, then it takes the average of the middle two depths in sorted_values
        else
            median="${sorted_values[median_index]}" #If the total number is odd it takes the middle element directly
        fi

        sum=0 #Calculates the mean/average value to use in the variance calculation below 
        for val in "${values[@]}"; do #Adds all of the depth values together 
            sum=$((sum + val))
        done
        avg=$((sum / ${#values[@]})) #Average is the sum divided by the number of values

        variance=0 #To calculate the variance, the difference between each depth value is compared to the average and the difference is squared, the squared differences are all added to the variance variable to accumulate the sum of squared differences
        for val in "${values[@]}"; do
            diff=$((val - avg))
            variance=$((variance + (diff * diff)))
        done
        sample_variance=$(awk "BEGIN {print $variance / (${#values[@]} - 1)}") #The sample variance is the sum of squared differences divided by the number of values minus 1 
        standard_deviation=$(awk "BEGIN {print sqrt($variance / ${#values[@]})}") #The standard deviation is the square root of the variance which is the sum of squared differences divided by the number of values  

        min=$(printf '%s\n' "${values[@]}" | sort -n | head -1) #The minimum depth is the lowest value when sorted
        max=$(printf '%s\n' "${values[@]}" | sort -n | tail -1) #The maximum depth is the highest value when sorted 

        echo -e "$window_name\t$mean\t$mode\t$median\t$sample_variance\t$standard_deviation\t$min\t$max" >> "$output_file" #For each window the name and all of the averages calculated are recorded on a new line in the output file averages_per_window.txt 

        rm "$window_bed" "${window_name_clean}_read_names.txt" "header_${window_name_clean}.sam" #Unwanted files are removed 
        
    else #This statement occurs if the bed file for the specified window had a size of 0, therefore didnt have any reads mapping across it. The corresponding window is skipped and the script continues. 
        echo "Skipping empty window: $window_name"
    fi

done < "${isolate}_WG_windows_query.bed" #This done statement names the WG_query.bed file that contains all of the window names as the input for the loop. This means that each window specified in this file will have the windows made and averages calculated. 

