#!/bin/bash
#PBS -N E3_INDIV.MODC.ARRAY  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=20:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)
#PBS -t 1-38 #array job



#Set filters
POP="MODC"
PVAL="0.05"
MINQ="20"
minMAPQ="20"
minDP="2"
REGION="LR761675.1:"
FILEPATH="../02a_modern_mapped/01a_modern_cutadapt_reads/"

#run job in working directory
cd $PBS_O_WORKDIR 

#load modules
module load languages/gcc-6.1
angsd=~/bin/angsd/angsd

#Define variables
INDIV=$(sed "${PBS_ARRAYID}q;d" $POP.poplist.short)


#estimate SFS for modern expanding population using ANGSD

time $angsd -i $FILEPATH${INDIV} -checkBamHeaders 1 -minQ $MINQ -minMapQ $minMAPQ -uniqueOnly 1 -remove_bads 1 \
-only_proper_pairs 1 -r $REGION -GL 1 -out HET.per.INDIV/$POP.${INDIV} \
-doSaf 1 -ref ../RefGenome/*fna -anc ../RefGenome/*fna -rmTriallelic 1 \
-doCounts 1 -dumpCounts 2 -doMajorMinor 4 -doMaf 1\
 -setMinDepthInd $minDP -baq 1 -C 50
