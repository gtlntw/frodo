#!/usr/bin/python
import os

#directories to ieterate
dir_list = ["2g2c", "2g2c_null", "2g3c", "2g3c_null", "3g_5ppl", "3g_5ppl_null"]

for i in dir_list:
	dir_path = "./"+i
	os.chdir(dir_path)
	
	#run
	os.system("R --vanilla < paste_mainSim_results.R")
	
	result_name=("2014-10-25 simulation for family data 1kb region more generations and children random risk haplotypes_")
	os.rename("allResults.csv", result_name+i+".csv")
	
	#copy files to upper level
	cmd ="cp \""+result_name+i+".csv\" "+"../"
	os.system(cmd)
	
	os.chdir("../")