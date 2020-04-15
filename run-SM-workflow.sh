#!/bin/bash

module load viz
module load graphviz
module load python/3.6.1
CF=config.json
CL=cluster.json
SF=cr-snakefile-v2.py

[ ! -d SM ] && mkdir -p SM

snakemake --snakefile $SF --unlock --cores 1 -np

snakemake --snakefile $SF --cluster-config $CL --configfile $CF --cores 256 -np


snakemake --snakefile $SF --cluster-config $CL --configfile $CF --cores 256 --dag | dot -Tsvg > dag-SM-workflow.svg

echo "open svg of DAG in browser and check if workflow looks ok"
read -p "Type Y to run if DAG looks ok, press any other key if not: " key
if [ "$key" = 'Y' ]; then
	nohup snakemake -p --snakefile $SF --cluster-config $CL --configfile $CF --max-jobs-per-second 5 --max-status-checks-per-second 5 -w 15  --cores 256 --cluster "sbatch --mem-per-cpu={cluster.mem} -t {cluster.time} --mail-user={cluster.email} --partition {cluster.partition} --mail-type={cluster.mail_type} --output {cluster.output}" && 
	rm SM/*%j.log &&
	[ ! -d SM/logs ] && mkdir -p SM/logs &&
	mv SM/*.log SM/logs &
else
	echo "not run"
fi
