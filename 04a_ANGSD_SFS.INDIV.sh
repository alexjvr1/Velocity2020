#!/bin/bash
#PBS -N E3_INDIV.MUS.ARRAY  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=20:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)
#PBS -t 1-48 #array job

#Define variable
NAMES="MUS.idx.names"

#run job in working directory
cd $PBS_O_WORKDIR 

#load modules
module load languages/gcc-6.1
realSFS=~/bin/angsd/misc/realSFS

#Define variables
INDIV=$(sed "${PBS_ARRAYID}q;d" MUS.idx.names)

##RUN ANGSD

$realSFS ${INDIV} > ${INDIV}.sfs

