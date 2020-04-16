#!/bin/bash
###########################################
# (c) Alexandra Jansen van Rensburg
# last modified 12/07/2019 05:49 
###########################################

## Index all bamfiles listed in bamlist

#PBS -N E1.index  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=10:00:00 ##wall time.
#PBS -j oe  #concatenates error and output files (with prefix job1)
#PBS -t 1-48

#run job in working directory
cd $PBS_O_WORKDIR


#load modules

module load apps/samtools-1.8

##Set up array

NAME=$(sed "${PBS_ARRAYID}q;d" bamlist)

##Run script
echo "Indexing ${NAME}"
printf "\n"

echo "time samtools index ${NAME}"
time samtools index ${NAME}
