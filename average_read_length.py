#this script gives the average read length of all the fastq files (after filtering)
#this script is used with the seq_analysis.sh and after_correction.sh

path = "/NGS/active/IPL/MENINGO/analysis/paloma/2023/"
isolate = input("input isolate name ")
filename = "/listreadlengths.txt"
openfile = open (path+isolate+filename, 'r')
rfile = openfile.readlines()
   
count = 0 
total = 0 
  
for lengths in rfile: 
        lengths.replace('\n','') 
        total += int(lengths) 
        count += 1 
average = total/count 
print (str('the average mapped read length for'), isolate, str(' is  '), average) 
