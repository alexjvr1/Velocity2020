#!/bin/bash
#PBS -N D3.BWA_mus  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=10:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)
#PBS -t 1-48

#run job in working directory
cd $PBS_O_WORKDIR 

#Load modules
module load apps/bwa-0.7.15
module load apps/samtools-1.9.1

#Define variables

RefSeq=GCF_905163445.1_ilParAegt1.1_genomic.fna
total_files=`find 01a_mus.concat_cutadapt_reads/ -name '*.fastq.gz' | wc -l`


NAME1=$(sed "${PBS_ARRAYID}q;d" mus.merged.names)

echo "mapping started" >> map.log
echo "---------------" >> map.log

##Check if Ref Genome is indexed by bwa
if [[ ! RefGenome/$RefSeq.fai ]]
then 
	echo $RefSeq" not indexed. Indexing now"
	bwa index RefGenome/$RefSeq
else
	echo $RefSeq" indexed"
fi


##Map with BWA MEM and output sorted bam file

sample_name=`echo ${NAME1} | awk -F ".repaired" '{print $1}'`
echo "[mapping running for] $sample_name"
printf "\n"
echo "time bwa mem RefGenome/$RefSeq 01d_musAll_merged/${NAME1} | samtools sort -o 02a_museum_mapped/${NAME1}.bam" >> map.log
time bwa mem RefGenome/$RefSeq 01d_musAll_merged/${NAME1} | samtools sort -o 02a_museum_mapped/${NAME1}.bam
