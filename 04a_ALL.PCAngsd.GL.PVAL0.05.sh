#!/bin/bash
#PBS -N E3_ALL.PCAngsd.GL  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=20:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)
#PBS -t 1-87 #array job

#Set filters
MININD=""
MINMAF=""
PVAL="0.05"
MINQ="20"
minMAPQ="20"
minDP="2"
maxDP=""
C="50"
POPLIST="ALL3POPS.poplist"
SPECIESDIR="/newhome/aj18951/E3_Aphantopus_hyperantus_2020"

#run job in working directory
cd $PBS_O_WORKDIR 

#load modules
module load languages/gcc-6.1
angsd=~/bin/angsd/angsd

#Define variables
REGION=$(sed "${PBS_ARRAYID}q;d" regions)


#estimate GL for modern expanding population using ANGSD

time $angsd -b $POPLIST -checkBamHeaders 1 -minQ 20 -minMapQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -r ${REGION} \
-GL 1 -doGlf 2 -out $SPECIESDIR/04a_ANGSD_FINAL/PCAngsd/GL.ALL3POPS.${REGION} -doSaf 1 \
-anc $SPECIESDIR/RefGenome/*fna -ref $SPECIESDIR/RefGenome/*fna \
-rmTriallelic 1 -doCounts 1 -dumpCounts 2 -doMajorMinor 4 -doMaf 1 -setMinDepthInd 2 \
-C $C -baq 1 -doGeno 32 -doPost 1 -SNP_pval $PVAL
