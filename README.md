# Velocity2020
NERC Velocity Project analyses using Sanger genomes

Here I'll curate the variant calling pipeline and analyses undertaken using the Lepidoptera genomes generated by Sanger in 2019/2020. 

## Pipeline

### 1. Raw to cleaned and processed data

   1a. [Trim adapter sequence using cutadapt](README.md#1a-demultiplex-and-adapter-trimming)
        
   1b. [Concatenate resequenced museum data](README.md#1b-concatenate-museum-reseq-data) (some individuals have been sequenced >1)
        
   1c. [Repair problems in museum PE data for data from 1.2. (BBrepair)](README.md#1c-prepare-museum-data-for-mapdamage-22-repair-pe-reads)
        
   1d. [Merge overlapping PE reads in museum data (BBmerge)](README.md#1d-prepare-museum-data-for-mapdamage-22-merge-overlapping-pe-reads)
        
#### 2. Map and process

   2a. [Map museum and modern data to Sanger genome](README.md#2a-map-to-reference-genome) 
        
   2b. [Correct museum data for possible deamination](README.md#2b-mapdamage-run-on-museum-data) (MapDamage -> output = corrected bam file) 
        
   2c. [Downsample modern data to the same depth as the museum data](README.md#2c-downsample-modern-data-to-the-same-coverage-as-in-the-museum-samples)
        
#### 3. [SNP discovery and filtering](https://github.com/alexjvr1/Velocity2020#3-angsd)

   3a. [ANGSD filters for SFS](README.md#3a-angsd-filters-for-sfs) (ie. no MAF)
        
   3b. [ANGSD filters for population genomics]
        
#### 4. Analyses: Outliers

   4a. [Outlier analysis in ANGSD]
        
   4b. [Manhattan plot]
        
   4c. [Table of functions]
        
   4d. [Network analysis]

#### 5. Analysis: Genetic diversity and population structure

#### 6. Analysis: LD 

   6a. [ANGSD estimate LD across the genome]


## DATA: Genome

Aphantopus hyperantus (Ringlet) was the first genome available ([NCBI link](https://www.ncbi.nlm.nih.gov/assembly/GCA_902806685.1)), so the pipeline will be set up with this species. 

## DATA: WGS

Whole genome resequencing data was generated for 38 & 40 modern individuals (sampled 2016-2017 & 2019) from a core and expanding population. Museum data was generated from 48 individuals + resequencing of a subset of individuals to increase read coverage. 

### Note on renaming files

Liverpool raw data is named with digits and a dash before the sample names. e.g. 33-AH-01-1900-47_191121_L001_R2.fastq.gz
The easiest way to rename them is with the Perl rename (note that the native linux rename works the same as mv and is not so useful in this case). 
Install the perl script (a version curated [here](https://github.com/subogero/rename)). On bluecp3 I've installed this in my software folder: /newhome/aj18951/software/rename-master/rename.

This works with the sed syntax. e.g. this will remove numbers plus dash from the start of the file names: 

```
../../software/rename-master/rename 's/^[0-9]+-//' *
```

## Pipeline for Velocity project from raw data to mapped reads:

### 1. Raw to cleaned and processed data

#### 1a Demultiplex and Adapter trimming

##### *TIME:*

This runs in 1-2 hours for the full dataset (museum + modern)


##### *METHOD:*

Modern samples arrive demultiplexed by the sequencing facility, but Museum samples need to be demultiplexed. 

We're trimming all adapter sequence from the demultiplexed data. We're also removing all sequences that are shorter than 20bp and 3' quality trimmed to remove bases with PHRED quality score of < 20 with Cutadapt.

If you're running this on BlueCrystal, you'll have to install cutadapt locally first. 
We're using cutadapt version 1.12: 
```
module load languages/python-anaconda3-5.2.0

#install cutadapt in your home directory using the web instructions
pip3 install --user --upgrade cutadapt

#Check that this cutadapt works
~/.local/bin/cutadapt --help


##Check if this directory is in your PATH:
echo $PATH

##And add to PATH if it isn't yet
PATH="$PATH:~/.local/bin/"

##Now you can run cutadapt directly
cutadapt --help

cutadapt version 1.12
```

##### Generate submission script

Run the following scripts to generate the submission script

[01a_museum_cutadapt_filtering_trimming.sh](https://github.com/alexjvr1/UKButterflies/blob/master/01a_museum_cutadapt_filtering_trimming.sh)

[01a_modern_cutadapt_filtering_trimming.sh](https://github.com/alexjvr1/Velocity2020/blob/master/01a_modern_cutadapt_filtering_trimming.sh)


##### Submission script

The above scripts will generate this: 

01a_parallel_cutadapt_bluecp3.sh

Edit this submission script to submit from your home directory:

```
1. Set all paths to your home directory if necessary. 

2. Adjust the number of threads (PBS -t 1-xx) to equal the number of individuals to be analysed. 

3. Check that any empty arguments have been removed from the cutadapt command

4. You might have to set the path to cutadapt to find your local version

```

#### 1b Concatenate museum reseq data

##### *TIME*

~30-40min

##### *METHOD*

A subset of individuals (33 per species) have been sequenced twice to increase mean depth. The data from both sequencing runs need to be concatenated together after adapter trimming. We're using these scripts: 

[1b_concat.fastq.R1.sh](https://github.com/alexjvr1/Velocity2020/blob/master/concat.fastq.R1.sh) and [1b_concat.fastq.R2.sh](https://github.com/alexjvr1/Velocity2020/blob/master/concat.fastq.R2.sh)

Reseq data are kept in the following folders:
```
00_raw_data_museum2

01a_museum2_cutadapt_reads

01a_mus.concat_cutadapt_reads  ## concatenated museum1 and museum2 + all samples that didn't have reseq data added. I'll point to this folder when mapping

02a_museum_mapped  ##see below. This contains all data including concatenated reseq samples. 
```


#### Rename samples

We don't need the really long names given to the samples by the sequencing facilities. We'll rename samples before proceeding further: 

Museum
```
cd 01a_mus.concat_cutadapt_reads/

#move all the log files into a folder
mkdir logfiles
mv *log logfiles

#change the names
~/software/rename_master/rename 's/long_name/new_name/' *gz

##remove all the extra info about lane number and date etc. Final names will bein this format: 

AH-01-1900-01_R1.concat.fastq.gz  AH-01-1900-17_R2.fastq.gz         AH-01-1900-34_R1.concat.fastq.gz
AH-01-1900-01_R2.concat.fastq.gz  AH-01-1900-18_R1.fastq.gz         AH-01-1900-34_R2.concat.fastq.gz

```

Modern
```
I'm mapping simultaneously to mod.core and mod.exp, but I'll keep the cutadapt reads in their different folders. 

cd 01a_modern_cutadapt_reads
#move all the log files into a folder
mkdir logfiles
mv *log logfiles

#change the names

~/software/rename-master/rename 's/R1_001.fastq.gzcutadapt_filtered/mod.core/' *gz
~/software/rename-master/rename 's/R2_001.fastq.gzcutadapt_filtered/mod.core/' *gz

AH-01-2016-01_mod.core_R1.fastq.gz  AH-01-2016-16_mod.core_R1.fastq.gz  AH-01-2017-29_mod.core_R1.fastq.gz
AH-01-2016-01_mod.core_R2.fastq.gz  AH-01-2016-16_mod.core_R2.fastq.gz  AH-01-2017-29_mod.core_R2.fastq.gz


cd  01a_modern.exp_cutadapt_reads
#move all the log files into a folder
mkdir logfiles
mv *log logfiles

~/software/rename-master/rename 's/R1_001.fastq.gzcutadapt_filtered/mod.exp/' *gz
~/software/rename-master/rename 's/R2_001.fastq.gzcutadapt_filtered/mod.exp/' *gz

AH-02-2019-42_mod.exp_R1.fastq.gz  AH-02-2019-57_mod.exp_R1.fastq.gz  AH-02-2019-71_mod.exp_R1.fastq.gz
AH-02-2019-42_mod.exp_R2.fastq.gz  AH-02-2019-57_mod.exp_R2.fastq.gz  AH-02-2019-71_mod.exp_R2.fastq.gz

```



#### 1c Prepare museum data for MapDamage (2.2): Repair PE reads

##### *TIME* 

~1hour

##### *METHOD*

Museum data is prone to post-mortem damage which we will correct for using MapDamage (see 2.2 below). We need to pre-process the museum reads for this. First we need to correct any problems with the PE files. The most common problem we've found is mismatches between R1 and R2 files e.g. not equal in length or mismatches between names. We can correct for this using BBtools's repair.sh script. 

[BBtools](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/) can be installed from the downloaded tarball

We have it installed here
```
/newhome/bzzjrb/Software/bbmap/
```

We'll repair the museum data using this script: [01c_bbtools_repair_museum_ARRAY.sh](https://github.com/alexjvr1/Velocity2020/blob/master/01c_bbtools_repair_museum_ARRAY.sh)

The input files will be all the museum fastq files found in 01a_mus.concat_cutadapt_reads

Create two files with names for the inputs in the home directory for the species (e.g. /E3_Aphantopus_hyperantus_2020/): 
```
ls 01a_mus.concat_cutadapt_reads/*R1*fastq.gz >> R1.museum.names.torepair
ls 01a_mus.concat_cutadapt_reads/*R2*fastq.gz >> R2.museum.names.torepair

##remove the file path from the name files. We need only in the input file names. 
sed -i 's:01a_mus.concat_cutadapt_reads/::g' *torepair
```

This will write all the repaired files to: 01c_musPERepaired

#### 1d Prepare museum data for MapDamage (2.2): Merge overlapping PE reads

##### *TIME*

30-40min

##### *METHOD*

MapDamage requires that overlapping PE reads be merged. Our museum libraries were sequenced with 75bp PE kits, except for the final (museum4) library which was sequenced with 50bp PE. 

We know from the library prep that the museum DNA is quite degraded. Thus we expect that the vast majority of the inserts would have overlapping PE sequences. 

We can count and merge these using bbtools' merge script. 

Edit [01d_bbtools_merge_museum_ARRAY.sh](https://github.com/alexjvr1/Velocity2020/blob/master/01d_bbtools_merge_museum_ARRAY.sh) to set up the array submission script. 

Make the output directory and we need the names files again: 
```
mkdir 01d_musAll_merged

ls 01c_musPERepaired/*R1*gz >> R1.museum.names.repaired
ls 01c_musPERepaired/*R2*gz >> R2.museum.names.repaired

sed -i 's:01c_musPERepaired/::g' *repaired
```

### 2. Map and process

Map museum and modern samples to the Sanger reference genome. Thereafter correct the museum bam files using MapDamage, and downsample the modern data to the same final depth as the corrected museum bam files. 

#### 2a Map to Reference Genome

##### *TIME:*

Museum samples (n=48) ~3 hours 

Modern Core (n=38) ~10 hours 

Modern Exp (n=40) ~10 hours for all but two samples which had to be restarted. They then ran in 4 hours... 

##### *METHOD:*

It is more efficient to run this code in local directory before submitting the mapping script to queue

```
#Index the reference genome if needed. Check if the *fasta.fai* file exists in the SpeciesName/RefGenome/ folder in your local directory. If not, run the indexing code. 

#index reference genome
module load apps/bwa-0.7.15
bwa index RefGenome/*fasta


#Create files with input names
## museum
ls 01a_mus.concat_cutadapt_reads/*R1*fastq.gz >> R1.museum.names
sed -i s:01a_mus.concat_cutadapt_reads/::g R1.museum.names

ls 01a_mus.concat_cutadapt_reads/*R2*fastq.gz >> R2.museum.names
sed -i s:01a_mus.concat_cutadapt_reads/::g R2.museum.names

## modern
#We're pointing to two input folders so I'll leave the path in the sample names folder
ls 01a_modern_cutadapt_reads/*R1* >> R1.modern.names
ls 01d_musAll_merged/*R1* >> R1.museum.names 

ls 01a_modern_cutadapt_reads/*R2* >> R2.modern.names
ls 01d_musAll_merged/*R2* >> R2.museum.names 

sed -i 's:01a_modern_cutadapt_reads/::g' *names
sed -i 's:01d_musAll_merged/::g' *names

#make output directories. 
mkdir 02a_museum_mapped
mkdir 02a_modern_mapped

#Add the additional folders in the 02a_modern_mapped folder when working with an EXPANDING species as the output will be written there. 
mkdir 02a_modern_mapped/01a_modern_cutadapt_reads
mkdir 02a_modern_mapped/01a_modern.exp_cutadapt_reads

#Check that you're pointing to the correct reference genome

#Check that the file separator makes sense: 
##sample_name=`echo ${NAME1} | awk -F "_R" '{print $1}'`
#Change the -F "xxx" according to the file names. 
#e.g the above works for files named as follows: 
#HS-01-2016-26_L007_cutadapt_filtered_R2.fastq.gz
#we want only the first part of this name to carry through. 
```

Run the submission scripts: 

[02a_MapwithBWAmem.ARRAY_museum.sh](https://github.com/alexjvr1/Velocity2020/blob/master/02a_MapwithBWAmem.ARRAY_museum.sh)

[02a_MapwithBWAmem.ARRAY_modern.sh](https://github.com/alexjvr1/Velocity2020/blob/master/02a_MapwithBWAmem.ARRAY_modern.sh)

Check that everything has mapped correctly by checking the file sizes. If the mapping is cut short (e.g. by exceeding the requested walltime) the partial bam file will look complete and can be indexed. But the bam file size will be small (~500kb) and empty when you look at it.
```
#To determine file size

du -sh *bam   

#To see bam file
module load apps/bcftools-1.8
bcftools view file.bam | head
Check the output with samtools flagstat

module load apps/samtools-1.8
samtools flagstat file.bam

#make a flagstat log file for all of the samples
for i in $(ls *bam); do ls $i >>flagstat.log && samtools flagstat $i >> flagstat.log; done
```

Index the bam files with the script [02a_index.bamfiles.sh](https://github.com/alexjvr1/Velocity2020/blob/master/02a_index.bamfiles.sh)



#### 2b. MapDamage run on museum data

##### *TIME*

3-4 hours for 48 samples

##### *METHOD*

[MapDamage2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3694634/) is a package used to estimate and correct for Cytosine deamination (or any other transition/transversion bias in the data). This is a problem anticipated for ancient DNA, and possibly for museum data.

This needs to be locally installed on BlueCrystal. Follow the instructions in the tutorial

MapDamage will be run in the 02a_museum_mapped folder. 


1. Create a file listing all the bamfiles
```

ls *bam > bamfiles.mus.names
```

2. Copy the script [02b_mapDamage_museum.sh](https://github.com/alexjvr1/Velocity2020/blob/master/03a_mapDamage_museum.sh) to the 02a_museum_mapped folder. Change the job name, the number of threads, and check the path to the reference genome.

Submit to queue.


Move all the new rescaled bam files to a new folder: 
```
mkdir 02b_museum_mapdamage.bams
mv 02a_museum_processed.mapped/results*/*bam 02b_museum_mapdamage.bams && cd 02b_museum_mapdamage.bams

module load apps/samtools-1.9
for i in $(ls *bam); do ls $i >> flagstat.log && samtools flagstat $i >> flagstat.log; done
```



#### 2c. Downsample modern data to the same coverage as in the museum samples

Due to the difference in sample quality between museum and modern samples, mean coverage is much higher for the modern data. This may bias the confidence in variant calls downstream. To avoid this problem I will downsample the modern data to the same mean depth as the museum data.

First filter the bam files to include only reads with PHRED quality >20 and properly paired reads using the [02c_Filter_modern_bam_pp.PHRED20.sh](https://github.com/alexjvr1/Velocity2020/blob/master/02c_Filter_modern_bam_pp.PHRED20.sh) script.

We'll need the names file: 
```
ls 02a_modern_mapped/*bam >> bamfiles.mod.names
sed -i 's:02a_modern_mapped/::' bamfiles.mod.names
```

Use samtools flagstat to calculate the number of properly paired reads in the recalibrated and filtered museum files.

```
module load apps/samtools-1.8
for i in $(ls results*/*flt.bam); do ls $i >> mus.flagstat.log && samtools flagstat $i >> mus.flagstat.log; done
```
Do the same for the modern samples.

Enter these data in the "Rescaled.ProperlyPaired.Q20" column in the Velocity_MapingStatsPerSpecies_AJvR_20190604.xlsx sheet on Dropbox. Calculate the mean number of museum reads and the proportion of modern reads to downsample to.

Use the [02c_Downsample_mod_ARRAY.sh](https://github.com/alexjvr1/Velocity2020/blob/master/02c_Downsample_mod_ARRAY.sh) script to downsample the modern bam files. Remember to change the job name and the PROP variables and create the input file listing all the modern bams.


### 3. ANGSD


#### 3a. ANGSD filters for SFS

I'm starting with bam files, so there are already some filters on the mapping quality of the sequences. Prior to that there are some crude filters during the demultiplexing and trimming steps.

To speed up the analysis I will split the ANGSD run across the genome; i.e. all indivs will be analysed for regions 0-x in ARRAY1, x-x1 in ARRAY2, etc. For this we need to split the genome up into regions.

Find all the regions (i.e. all chromosomes and contigs) from the reference index file (.fasta.fai):

```
awk '{print $1}' ../RefGenome/*.fna.fai >> regions

cat regions |wc -l
>87

##We'll run these in an ARRAY. 
```

###### *Filters*

-b[filelist]

-remove_bads 1 : remove reads with 255 flag (not primary, failure and duplicate reads) (1=default)

-uniqueOnly 1 : remove reads with multiple best hits

-minMapQ 20 : PHRED 20. This should already be in place during the mapping.

-minQ 20 : PHRED 20 for individual base score.

-only_proper_pairs 1 : include only properly paired reads (default) and should already have been applied to the museum reads prior to this. 

-trim 0 : We're not trimming any data

-baq 1 : estimate base alignment quality using samtools method.

###ALLELE FREQUENCY ESTIMATION

-doMajorMinor 4 : Force Major allele based on reference. The minor allele is then inferred using doMajorMinor 1. This option needs to be used when calculating SFS for multiple populations as ANGSD otherwise determines a minor allele within each population. I.e. this may not be the same across all the populations.

-ref [..fasta] : For doMM 4 above we need to specify a reference genome.

-doMaf 1 : calculate minor allele frequency

-SNP_pval 0.001 : Only work with SNPs with a p-value above [float]

-GL 1 : I will estimate genotype likelihoods using the SAMtools model

-minInd 18 : I will remove loci where less than 18 individuals have been genotyped. There are 19-45 indivs per group (HOD,FOR, South, New). So this number seems quite high, but this is to be more certain of the 5% MAF.

-setMinDepth : Discard site if total depth (across all indivs) is below [int]. Use -doCounts to determine the distribution of depths

-setMaxDepth : Discard site if total depth (across all indivs) is above [int]

-setMinDepthInd 2 : Minimum depth for a locus for an individual. This is only applicable for analyses using counts (-doCounts obligatory)

-setMaxDepthInd [int]: I'll use a max of meanDP + 2xSD of depth.

-doCounts 1 : Count of the nucleotide bases per individual

-dumpCounts 2 : write a file with all the allele counts per position per individual

-rmTriallelic 1 : include only biallelic loci

-checkBamHeaders 1 : check that the bam headers are compatable for all files.


#### 3b. Call GL

#### 3c. SFS



#### 3d. ISSUES

##### *1. Input bam files*

I've had problems with reading merged bam files created by BBmerge into R. I've reported the issue [here](https://github.com/ANGSD/angsd/issues/260) - it seems to be the same problem encountered by someone using inputs using miniMap2

```
[bammer_main] 1 samples in 1 input files
-> Parsing 1 number of samples
No data for chromoId=0 chromoname=LR761647.1
This could either indicate that there really is no data for this chromosome
Or it could be problem with this program regSize=0 notDone=0

-> Done reading data waiting for calculations to finish
-> Done waiting for threads
-> Output filenames:
	->"angsdput.arg"
-> Fri Apr 17 16:45:14 2020

-> Arguments and parameters for all analysis are located in .arg file
-> Total number of sites analyzed: 0
-> Number of sites retained after filtering: 0 
[ALL done] cpu-time used =  4.24 sec
[ALL done] walltime used =  4.00 sec
```

I thought bbmerge might be the problem, so I will try to read in the example files from: 

1. BBmerge

2. PEAR

3. MapDamage

```
pwd
/newhome/aj18951/E3_Aphantopus_hyperantus_2020/ANGSDinputTests

#convert example sam to bam  
samtools view -S -b ~/software/mapDamage/mapdamage/rescale_test/pe_test/pe_rescaled_correct.sam > mapDamage_pe_rescaled_correct.bam

#move the reference fasta file to the same folder and index
cp ~/software/mapDamage/mapdamage/rescale_test/pe_test/ref.fa .
module load apps/bwa-0.7.15
samtools faidx ref.fa

#submit the following script to bluecrystal p3

> cat angsd.test.sh
#!/bin/bash
#PBS -N angsd.test  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=20:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)


#run job in working directory
cd $PBS_O_WORKDIR 

#load modules
module load languages/gcc-6.1
angsd=~/bin/angsd/angsd


time $angsd -i mapDamage_pe_rescaled_correct.bam -ref ref.fa -GL -out mapdamage.test.out
```





#### 4d. Compare downsampled data to full dataset for modern pops (SFS and GL)

#### 4e. Population structure (PCA)

#### 4f. Diversity stats 

#### 4g. Outliers

#### 4h. LD analyses


### 5. Map outlier loci

Map outlier loci, identify candidates in the area

Synteny between modern and museum samples. 





