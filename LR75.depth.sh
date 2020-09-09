#!/bin/bash
###########################################
# (c) Alexandra Jansen van Rensburg
# last modified 12/07/2019 05:49 
###########################################

## Index all bamfiles listed in bamlist

#PBS -N LR75  ##job name
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
echo "Output depth stats for ${NAME}"
printf "\n"

echo "time samtools view -b ${NAME} LR761675.1 > ${NAME}.LR75.bam" &&\
time samtools view -b ${NAME} LR761675.1 > ${NAME}.LR75.bam && \
echo "time samtools index ${NAME}.LR75.bam" && \
time samtools index ${NAME}.LR75.bam && \
echo "time samtools depth ${NAME}.LR75.bam > ${NAME}.LR75.depth" &&\
time samtools depth ${NAME} > ${NAME}.LR75.depth
