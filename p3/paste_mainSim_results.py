#!/usr/bin/env python
#!/usr/bin/python
#!python

from __future__ import print_function
import sys
import os
import re
import getopt
import numpy

outfile = "allResults.csv"
csvList = []
n_proccessed = 1

##Read in list of CSV
for file in os.listdir("./"):
	if file.endswith(".csv"):
		if file != outfile: #exclude allResults.csv
			csvList.append(file)

# print csvList
print(str(len(csvList)) + " files to combine\n")
#what to do if only one file
if len(csvList) == 1:
	print('there is only 1 file no need to combine')
	sys.exit(0)

##read in the first file
infh= open(csvList[0])
print("{0},".format(str(n_proccessed)), end="")
n_proccessed +=1
#keep the name and the dimension
var_names = infh.readline().strip().split(",")
#output to outfile
ofh = open(outfile, "w")
ofh.write(",".join(var_names)+"\n")
for line in infh.xreadlines():
	ofh.write(line)
infh.close()
#read in the rest of the files
for file_name in csvList[1:]:
	infh= open(file_name)
	print("{0},".format(str(n_proccessed)), end="")
	n_proccessed +=1
	var_names_other = infh.readline().strip().split(",")
	#check the name or dimension 
	if len(var_names) != len(var_names_other):
		os.remove(outfile)
		sys.exit("number of variables differs in at least one file")
	if var_names != var_names_other:
		os.remove(outfile)
		sys.exit("colnames differ in at least one file")
	#output to outfile
	for line in infh.xreadlines():
		ofh.write(line)
print("\n")


# ##read in each and output to combined file
# infh= open(inputName)
# for line in infh.xreadlines():
#     fields = line.split()
# dosage=[]
# switch_pos=[]
# for item in fields:
# 	a, b=((item.split(':')))
# 	dosage.append(a)
# 	switch_pos.append(b)
# infh.close()
# ofh = open(outputName, "w")
# start_pos_index = 0
# for idx, value in enumerate(switch_pos):
# 	switch_pos_index =int(value)-1
# 	ancestry = mapping[str(dosage[idx])]
# 	for i in pos[start_pos_index:switch_pos_index]:
# 		ofh.write(ancestry+"\n")
# 	start_pos_index = switch_pos_index;
# ofh.write(mapping[str(dosage[-1])]+"\n") #print the last snp
# ofh.close()
# #check dimension
# #check name

# ofh = open(outputName, "w")
# start_pos_index = 0
# for idx, value in enumerate(switch_pos):
# 	switch_pos_index =int(value)-1
# 	ancestry = mapping[str(dosage[idx])]
# 	for i in pos[start_pos_index:switch_pos_index]:
# 		ofh.write(ancestry+"\n")
# 	start_pos_index = switch_pos_index;
# ofh.write(mapping[str(dosage[-1])]+"\n") #print the last snp
# ofh.close()