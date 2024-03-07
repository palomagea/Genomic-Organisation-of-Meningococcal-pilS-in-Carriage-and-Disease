#this script gives the average depth of coverage across the genome and the standard deviation
#this script is run with the seq_analysis.sh or after_corrections.sh

import sys

#this variable is set in the seq_analysis.sh or after_corrections.sh scripts
coverage_file = sys.argv[1]

#the coverage.txt file which has the depth at each position in the genome is opened. this file has columns for contig, position and depth
openfile = open (coverage_file, 'r') 
rfile = openfile.readlines() 

#set empty variables
count = 0 
total = 0 
coveragelist = [] 

#the contig, positon and depth for each positon are on separate lines. for each line, the depth is added to the total and the count keeps track of how many positons (or lines in the file) have been processed
for lines in rfile: 
        lines.replace('\n','')
        field1,field2,field3 = lines.split('\t') #the formating of the lines is changed to separate them with a tab so they are able to be processed
        total += int (field3) #the depth at each position is read as an interger and added to the total
        count += 1 #for each line the count increases by 1 so the total number of positions is recorded
        coveragelist.append(int(field3)) #the depth at each positon is added to a list called coverage list
result = total /count  #the mean depth is calculated by dividing the total of all depths by the number of genomic positions
print (str('the average depth is '), result) #the mean depth is printed

import statistics 
std = statistics.stdev(coveragelist)  #the standard deviation is calculated from the coverage list that contains all the depths 
print (str('the standard deviation of coverage is '), std) #the standard deviation is printed 
openfile.close() 
