#this script gives the average depth of coverage across the genome and the standard deviation
#this script is run with the seq_analysis.sh or after_corrections.sh

path = "/NGS/active/IPL/MENINGO/analysis/paloma/2023/" 
isolate = input("input isolate name ") 
filename = "/" + isolate + "_coverage.txt" 
openfile = open (path+isolate+filename, 'r') 
rfile = openfile.readlines() 
  
count = 0 
total = 0 
coveragelist = [] 

for lines in rfile: 
        lines.replace('\n','') 
        field1,field2,field3 = lines.split('\t') 
        total += int (field3) 
        count += 1 
        coveragelist.append(int(field3))
result = total /count  
print (str('isolate '), isolate, str('the average depth of coverage is '), result) 

import statistics 
std = statistics.stdev(coveragelist) 
print (str('isolate '), isolate, str('the standard deviation of coverage is '), std) 
openfile.close() 
