#this script gives the average read length of all the fastq files (after filtering)
#this script is used with the seq_analysis.sh and after_correction.sh

import sys

#this variable is set in the seq_analysis.sh or after_corrections.sh scripts
read_lengths_file = sys.argv[1]

#the file that has a list of all the lengths of the reads that mapped to the genome is opened
openfile = open (read_lengths_file, 'r')
rfile = openfile.readlines()
   
count = 0 
total = 0 

#each length in the file is stored on a new line. for each length, the length is added to the total and the count goes up by 1 to keep track of how many lenghts have been added 
for lengths in rfile: 
        lengths.replace('\n','') #formatting changed 
        total += int(lengths) #the length is read as an interger and added to the total
        count += 1 #for each line in the file the count increases by one 
average = total/count #the average is calculated by the total length of all reads divided by the number of reads 
print (str('the average mapped read length for'), isolate, str(' is  '), average) #the average read length is printed 
