#!/bin/bash
#PBS -N E3_FoldedSFS  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=20:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)
#PBS -t 1-3 #array job

#run job in working directory
cd $PBS_O_WORKDIR 

#pop1.list
#MODC
#MODC
#MODE

#pop2.list
#MODE
#MUS
#MUS

#load modules
module load languages/gcc-6.1
angsd=~/bin/angsd/angsd

#Define variables
#REGION=$(sed "${PBS_ARRAYID}q;d" regions)
POP1=$(sed "${PBS_ARRAYID}q;d" pop1.list)
POP2=$(sed "${PBS_ARRAYID}q;d" pop2.list)
REGION="LR761675.1"
DATE="MAY26"

#estimate folded SFS for each pop pair using ANGSD

time ~/bin/angsd/misc/realSFS ${POP1}.$REGION.$DATE.saf.idx ${POP2}.$REGION.$DATE.saf.idx -fold 1 > ${POP1}.${POP2}.fold.sfs
