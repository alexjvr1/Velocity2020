#!/bin/bash
################################################################
# (c) Alexandra Jansen van Rensburg
# Last modified 11/07/2019 23:02
################################################################
## Downsample the number of reads in a bam file to the same coverage as in the museum samples
## PROP = mean museum properly paired reads/ number of modern indiv's properply paired reads. This is calculated in the VelocityMappingStatsPerSpecies_AJvR_20190604.xlsx file. Calculate mean number of reads in the museum bam files after filtering for Phred <20. Calculate proportion of modern reads to keep. 
## PROP input = column listing these proportions in the same order as the bam name list. 

#PBS -N E3.DwnSmpl
#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=8G
#PBS -j oe
#PBS -t 1-78

#load modules
module load languages/java-jdk-11.0.3
module load apps/picard-2.20.0
module load apps/samtools-1.8

#Run script in working directory 
cd $PBS_O_WORKDIR


#Filter all bamfiles to remove mapping quality <20 and retain only properly paired reads. 

#samtools view -f 2 -q 20 -b ./rescaled.bam/$RESCALED.BAM.NAME.bam > ./rescaled.bam/$RESCALED.BAM.NAME.flt.bam 


#Calculate proportion reads to remove from each bam
#Define variables

NAME=$(sed "${PBS_ARRAYID}q;d" bamfiles.mod.flt.names)
PROP=$(sed "${PBS_ARRAYID}q;d" Downsample.proportions)

echo "Downsampling started" >> 02c_downsample.log
echo "---------------" >> 02c_downsample.log

sample_name=`echo ${NAME} | awk -F "_" '{print $1}'`
echo "[filtering] $sample_name"
printf "\n"


#Remove reads to that proportion using PicardTools DownsampleSam
java -jar /cm/shared/apps/Picard-2.20.0/picard.jar \
DownsampleSam \
I=${NAME}.bam.flt.bam \
O=${NAME}.Downsampled.bam \
P=${PROP}
