#!/bin/bash
##################################
# Alexandra Jansen van Rensburg
# alexjvr@gmail.com
# Last modified 05/07/2021 15:56
##################################
#v1 - modified path to wrapper script

# Creates submission script for variant calling using mpileup
# using arrays on BlueCrystal p3

#run job in working directory
cd $SGE_O_WORKDIR

#load your program if it is installed globally or the modules you used to install your program locally (compilers, etc) 
#Specify modules

SHAREDFOLDER="/SAN/ugi/LepGenomics"
SPECIES="E3_Aphantopus_hyperantus"
INPUT="02c_MapDamage_MUS"
OUTPUT="03a_Variant_calling_museum"
REF="RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna"
SAMTOOLS="/share/apps/genomics/samtools-1.9/bin/samtools"
JOBNAME="E3_mod_mpileup"
CALLER="/SAN/ugi/LepGenomics/VelocityPipeline/wrapper/03a_call_SNVs_UCL.sh"

$CALLER -i $SHAREDFOLDER/$SPECIES/$INPUT \
-r $SHAREDFOLDER/$SPECIES/$REF -o $SHAREDFOLDER/$SPECIES/$OUTPUT \
-c c -v 1 -d 0 -s 20 -p 0.05 \
-job $JOBNAME;

