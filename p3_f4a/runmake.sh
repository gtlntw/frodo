#!/bin/bash
#!/usr/bin/env bash
#!/usr/bin/bash
#!bash
#default job_name = 'test'
launch='slurm';

while getopts i:j:l: opt; do
  case $opt in
 	i)
		d=$OPTARG
		;;
	j)
		job_name=$OPTARG
		;;
	l)
		launch=$OPTARG
		;;
  esac
done

#clean start and end log
rm -rf start*.log end*.log

#command to run
[ -z "$job_name" ] && echo need to give a job name \"./runmake.sh -j xxxxxx\" && exit 1
jobNum=500; #number of parallel jobs
echo "running ${jobNum} parallel jobs"
python ./makefile.py -l ${launch} -j ${job_name};

##run twice to make sure the preempted jobs are rerun
cat <<EOF > runmake_job.sh
#!/bin/bash
make -j ${jobNum} --keep-going -f makefile_${job_name} 
make -j ${jobNum} --keep-going -f makefile_${job_name} 
EOF
chmod +x runmake_job.sh
./runmake_job.sh > "runmake_${job_name}.log" 2>&1 &
