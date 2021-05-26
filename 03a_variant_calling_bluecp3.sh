#!/bin/bash
##################################
# Alexandra Jansen van Rensburg
# alexjvr@gmail.com
# Last modified 23/10/2018 09:35
##################################

# Creates submission script for variant calling using mpileup
# using arrays on BlueCrystal p3

#run job in working directory
cd $PBS_O_WORKDIR

#load your program if it is installed globally or the modules you used to install your program locally (compilers, etc) 
#Specify modules

HOME="/newhome/aj18951"
SPECIES="E3_Aphantopus_hyperantus_2020"
INPUT="02b_museum_mapDamage"
OUTPUT="03.2_Variant_calling_museum"
REF="RefGenome/*fna"
SAMTOOLS="apps/samtools-1.9.1"
JOBNAME="E3_mus_mpileup"
CALLER="bristol-velocity/AJvR_VelocityPipeline/wrapper/03a_call_SNVs_bluecp3.sh"


$HOME/$CALLER -i $HOME/$SPECIES/$INPUT \
-r $HOME/$SPECIES/$REF -o $HOME/$SPECIES/$OUTPUT \
-c c -v 1 -d 0 -s 20 -p 0.05 \
-module1 $SAMTOOLS -job $JOBNAME;
