#!/bin/bash
# -----------------------------------------

# (c) Alexandra Jansen van Rensburg
# alexjvr@gmail.com
# Last modified: 03/10/2018 15:39:48

# Description:
# Filter bam files to include only properly paired reads with PHRED > 20. 


#PBS -N E3.Flt.modBam  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=20:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)
#PBS -t 1-78

#run job in working directory
cd $PBS_O_WORKDIR 

# Load modules
module load apps/samtools-1.8

# Define variables
NAME=$(sed "${PBS_ARRAYID}q;d" bamfiles.mod.names)

##Script
echo "Filtering modern bam files: started" >> 02c_FltModBam.log
echo "---------------" >> 02c_FltModBam.log

sample_name=`echo ${NAME} | awk -F "" '{print $2}'`
echo "[mapping running for] $sample_name"
printf "\n"
echo "time samtools view -h -f 2 -q 20 ${NAME} -b -o ${NAME}.flt.bam >> 02c_FltModBam.log"
time samtools view -h -f 2 -q 20 ${NAME} -b -o ${NAME}.flt.bam
