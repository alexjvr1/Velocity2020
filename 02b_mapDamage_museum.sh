#!/bin/bash
# -----------------------------------------

# (c) Alexandra Jansen van Rensburg
# alexjvr@gmail.com
# Last modified: 03/10/2018 15:39:48

# Description:
# MapDamage to assess bias DNA damage in museum vs modern data and rescale mapping quality in .bam files accordingly. 
# Each samples runs in 30min - 5hours. The script completes a species museum and modern bam files in 8 hours on bluecrystal


#PBS -N E3.MapDmg.mus  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=20:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)
#PBS -t 1-48

#run job in working directory
cd $PBS_O_WORKDIR 

# Load modules
module load languages/python-2.7.6
module load languages/R-3.0.2

# Define variables
mapDamage="$HOME/.local/bin/mapDamage"
RefSeq="../RefGenome/*fna"
NAME=$(sed "${PBS_ARRAYID}q;d" bamlist)

##Script
echo "mapDamage2 started" >> map.log
echo "---------------" >> map.log

sample_name=`echo ${NAME} | awk -F ".bam" '{print $1}'`
echo "[mapping running for] $sample_name"
printf "\n"
echo "time $mapDamage --merge-reference-sequences -i ${NAME} -r $RefSeq >> mapdamage.log"
time $mapDamage --merge-reference-sequences -i ${NAME} -r $RefSeq --rescale
