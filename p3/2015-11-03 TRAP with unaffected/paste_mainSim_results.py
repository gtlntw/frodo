#!/usr/bin/env python
#!/usr/bin/python
#!python

import sys
import os
import re
import getopt
import numpy

outfile = "allResults.csv"
csvList = []

##Read in list of CSV
for file in os.listdir("./"):
	if file.endswith(".csv"):
		if file != outfile: #exclude allResults.csv
			csvList.append(file)

print csvList
#what to do if only one file
if len(csvList) == 1:
	cat 'there is only 1 file no need to combine'
	return 0

##read in each and output to combined file
infh= open(inputName)
for line in infh.xreadlines():
    fields = line.split()
dosage=[]
switch_pos=[]
for item in fields:
	a, b=((item.split(':')))
	dosage.append(a)
	switch_pos.append(b)
infh.close()
ofh = open(outputName, "w")
start_pos_index = 0
for idx, value in enumerate(switch_pos):
	switch_pos_index =int(value)-1
	ancestry = mapping[str(dosage[idx])]
	for i in pos[start_pos_index:switch_pos_index]:
		ofh.write(ancestry+"\n")
	start_pos_index = switch_pos_index;
ofh.write(mapping[str(dosage[-1])]+"\n") #print the last snp
ofh.close()
#check dimension
#check name

ofh = open(outputName, "w")
start_pos_index = 0
for idx, value in enumerate(switch_pos):
	switch_pos_index =int(value)-1
	ancestry = mapping[str(dosage[idx])]
	for i in pos[start_pos_index:switch_pos_index]:
		ofh.write(ancestry+"\n")
	start_pos_index = switch_pos_index;
ofh.write(mapping[str(dosage[-1])]+"\n") #print the last snp
ofh.close()