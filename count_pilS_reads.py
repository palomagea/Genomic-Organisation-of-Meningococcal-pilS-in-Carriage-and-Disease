isolate = input("input isolate name ") 
path = str("/NGS/active/IPL/MENINGO/analysis/paloma/2023/") 
filename = "/" + isolate + "_pilS.bed"

pilS_list = open (path+isolate+filename, 'r') 

count = 0 

for line in pilS_list: 

        if line != "\n": 

                count += 1 

print ("The total number of reads fully mapping across PilS are ", count) 
