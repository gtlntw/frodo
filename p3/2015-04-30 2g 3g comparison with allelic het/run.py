#!/usr/bin/python
import os
import subprocess

#directories to ieterate
dir_list = ["2g2c", "2g2c_null", "2g3c", "2g3c_null", "3g_5ppl", "3g_5ppl_null"]

for i in dir_list:
	dir_path = "./"+i
	os.chdir(dir_path)
	#copy mainSim.R
	cmd ="cp ../\"mainSim.R\"" + " ./"
	os.system(cmd)
	cmd ="cp ../\"run_mainSim_parallel.R\"" + " ./"
	os.system(cmd)
	#clean up
	os.system("rm *.csv")
	os.system("rm *.Rout*")
	
	if(i=="2g2c"):
		n_family="500"
		null="FALSE"
		family_strct="family_strct.2g2c"
	if(i=="2g2c_null"):
		n_family="500"
		null="TRUE"
		family_strct="family_strct.2g2c"
	if(i=="2g3c"):
		n_family="400"
		null="FALSE"
		family_strct="family_strct.2g3c"
	if(i=="2g3c_null"):
		n_family="400"
		null="TRUE"
		family_strct="family_strct.2g3c"
	if(i=="3g_5ppl"):
		n_family="400"
		null="FALSE"
		family_strct="family_strct.3g"
	if(i=="3g_5ppl_null"):
		n_family="400"
		null="TRUE"
		family_strct="family_strct.3g"	
	#run
	cmd = "R --vanilla --args n_family " + n_family + " null " + null + " family_strct " + family_strct + " < run_mainSim_parallel.R > run_mainSim_parallel.txt"
	print(cmd)
	process = subprocess.call(cmd, shell=True)
	
	os.chdir("../")