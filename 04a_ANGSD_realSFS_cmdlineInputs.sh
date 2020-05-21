#!/bin/bash
#PBS -N E3_SFSMOD_ANGSD1.ARRAY  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=20:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)

#run job in working directory
cd $PBS_O_WORKDIR 

INPUT=$1
OUTPUT=$2

#load modules
module load languages/gcc-6.1
angsd=~/bin/angsd/angsd
REALSFS=~/bin/angsd/misc/realSFS


time $REALSFS $INPUT -fold 1 > $OUTPUT
