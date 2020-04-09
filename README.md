# Velocity2020
NERC Velocity Project analyses using Sanger genomes

Here I'll curate the variant calling pipeline and analyses undertaken using the Lepidoptera genomes generated by Sanger in 2019/2020. 


## Genome

Aphantopus hyperantus (Ringlet) was the first genome available ([NCBI link](https://www.ncbi.nlm.nih.gov/assembly/GCA_902806685.1)), so the pipeline will be set up with this species. 

## WGS

Whole genome resequencing data was generated for 38 & 40 modern individuals (sampled 2016-2017 & 2019) from a core and expanding population. Museum data was generated from 48 individuals + resequencing of a subset of individuals to increase read coverage. 

## Note on renaming files

Liverpool raw data is named with digits and a dash before the sample names. e.g. 33-AH-01-1900-47_191121_L001_R2.fastq.gz
The easiest way to rename them is with the Perl rename (note that the native linux rename works the same as mv and is not so useful in this case). 
Install the perl script (a version curated [here](https://github.com/subogero/rename)). On bluecp3 I've installed this in my software folder: /newhome/aj18951/software/rename-master/rename.

This works with the sed syntax. e.g. this will remove numbers plus dash from the start of the file names: 

```
../../software/rename-master/rename 's/^[0-9]+-//' *
```

## Pipeline for Velocity project from raw data to mapped reads:

The scripts for the final pipeline are outlined below. The scripts are organised in 3 folders within a main folder.

Use this pipeline by running scripts in the pipeline/ folder in the numbered order. These scripts generate submission scripts for each step that can be submitted to the BlueCrystal p3 queue (i.e. qsub script.sh).

Inputs for each step should be submitted via the command line.

```
|
-----> pipeline  

        This contains all the scripts that generate submission scripts for BlueCrystal. Options can be specified in the command line. 

|
-----> wrapper
    
        Scripts called by pipeline scripts. Wrapper scripts for generating Queue request and specifying inputs from command line for                the tools to be called.     
        
|
-----> tools
        
        Scripts of tools or functions used in each step. These are called by the wrapper script
```


### 1. Demultiplex and Adapter trimming

#### *TIME:*

This runs in 1-2 hours for the full dataset (museum + modern)


#### *METHOD:*

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

##### *pipeline

01a_museum_cutadapt_filtering_trimming.sh

01a_modern_cutadapt_filtering_trimming.sh


##### *wrapper

01a_parallel_cutadapt_bluecp3.sh

Edit the generated script above to submit from your home directory:

```
1. Set all paths to your home directory if necessary. 

2. Adjust the number of threads (PBS -t 1-xx) to equal the number of individuals to be analysed. 

3. Check that any empty arguments have been removed from the cutadapt command

4. You might have to set the path to cutadapt to find your local version

```

To incorporate new data (e.g. resequencing of some individuals to increase mean depth), new fastq files need to be adapter trimmed. Fastq files are concatenated after this using the script [concat.fastq.R1.sh](https://github.com/alexjvr1/Velocity2020/blob/master/concat.fastq.R1.sh) and [concat.fastq.R2.sh](https://github.com/alexjvr1/Velocity2020/blob/master/concat.fastq.R2.sh)

Reseq data are kept in the following folders:
```
00_raw_data_museum2

01a_museum2_cutadapt_reads

01a_mus.concat_cutadapt_reads  ## concatenated museum1 and museum2 + all samples that didn't have reseq data added. I'll point to this folder when mapping

02a_museum_mapped  ##see below. This contains all data including concatenated reseq samples. 
```

#### Rename samples before mapping

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

### 2. Map to Reference Genome

#### *TIME:*

Museum samples (n=48) ~6 hours 

Modern Core (n=38) ~10 hours 

Modern Exp (n=40) ~10 hours for all but two samples which had to be restarted. They then ran in 4 hours... 

#### *METHOD:*

It is more efficient to run this code in local directory before submitting to queue

```
#Index the reference genome if needed. Check whether the *fasta.fai* file exists in the SpeciesName/RefGenome/ folder in your local directory. If not, run the indexing code. 

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
ls 01a_modern.exp_cutadapt_reads/*R1* >> R1.modern.names 

ls 01a_modern_cutadapt_reads/*R2* >> R2.modern.names
ls 01a_modern.exp_cutadapt_reads/*R2* >> R2.modern.names 


#make output directories. Remember to add the additional folders in the 02a_modern_mapped folder as the output will be written there. 
mkdir 02a_museum_mapped
mkdir 02a_modern_mapped
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

*pipeline*

[02_MapwithBWAmem.ARRAY_museum.sh](https://github.com/alexjvr1/Velocity2020/blob/master/02_MapwithBWAmem.ARRAY_museum.sh)

[02_MapwithBWAmem.ARRAY_modern.sh](https://github.com/alexjvr1/Velocity2020/blob/master/02_MapwithBWAmem.ARRAY_modern.sh)

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
Index the bam files with the script index.bamfiles.sh
```

Check that everything has mapped correctly by checking the file sizes. If the mapping is cut short (e.g. by exceeding the requested walltime) the partial bam file will look complete and can be indexed. But the bam file size will be small (~500kb) and empty when you look at it.

```
#To determine file size

du -sh *bam   

#To see bam file
module load apps/bcftools-1.8
bcftools view file.bam | head

```


Check the output with samtools flagstat
```
module load apps/samtools-1.8
samtools flagstat file.bam

#make a flagstat log file for all of the samples
for i in $(ls *bam); do ls $i >>flagstat.log && samtools flagstat $i >> flagstat.log; done
```

Index the bam files with the script [index.bamfiles.sh](https://github.com/alexjvr1/Velocity2020/blob/master/index.bamfiles.sh)



### 3. MapDamage: correct for Cytosine deamination in museum data

#### 3a. MapDamage run on museum data

[MapDamage2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3694634/) is a package used to estimate and correct for Cytosine deamination (or any other transition/transversion bias in the data). This is a problem anticipated for ancient DNA, and possibly for museum data.

This needs to be locally installed on BlueCrystal. Follow the instructions in the tutorial

MapDamage will be run in the 02a_museum_mapped folder. 


1. Create a file listing all the bamfiles
```

ls *bam > bamfiles.mus.names
```

2. Copy the script [03a_mapDamage_museum.sh](https://github.com/alexjvr1/Velocity2020/blob/master/03a_mapDamage_museum.sh) to the 02a_museum_mapped folder. Change the job name, the number of threads, and check the path to the reference genome.

Submit to queue.

Analyse output stats

#### 3b. Downsample modern data to the same coverage as in the museum samples


### 4. ANGSD

#### 4a. Set up filters (including min and maxDepth)

#### 4b. Call GL

#### 4c. SFS

#### 4d. Compare downsampled data to full dataset for modern pops (SFS and GL)

#### 4e. Population structure (PCA)

#### 4f. Diversity stats 

#### 4g. Outliers

#### 4h. LD analyses


### 5. Map outlier loci

Map outlier loci, identify candidates in the area

Synteny between modern and museum samples. 





