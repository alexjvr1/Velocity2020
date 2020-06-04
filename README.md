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

There seem to be improperly paired reads in my final dataset. To deal with this I need to change the ANGSD flag -only_proper_pairs to 0. 

Also make sure to use the latest commit of ANGSD. 

I'm using: 
```
angsd version: 0.933-18-gfd1a21a (htslib: 1.10.2-61-g8859b09) build(May  6 2020 14:42:05)
```

#### Set up analysis

To speed up the analysis I will split the ANGSD run across the genome; i.e. all indivs will be analysed for regions 0-x in ARRAY1, x-x1 in ARRAY2, etc. For this we need to split the genome up into regions.

Find all the regions (i.e. all chromosomes and contigs) from the reference index file (.fasta.fai):

```
awk '{print $1}' ../RefGenome/*.fna.fai >> regions

cat regions |wc -l
>87

##We'll run these in an ARRAY. 
```

The shortest contig is 11048 and longest is 18856181 for Ringlet

###### *Filters*

I'm starting with bam files, so there are already some filters on the mapping quality of the sequences. Prior to that there are some crude filters during the demultiplexing and trimming steps.


Other possible filters: 

-b[filelist]

-remove_bads 1 : remove reads with 255 flag (not primary, failure and duplicate reads) (1=default)

-uniqueOnly 1 : remove reads with multiple best hits

-minMapQ 20 : PHRED 20. This should already be in place during the mapping.

-minQ 20 : PHRED 20 for individual base score.

-only_proper_pairs 0 : NBNB THIS flag is changed to 0 because some of the reads in my final mus files are not properly paired! 
***OLD FILTER include only properly paired reads (default) and should already have been applied to the museum reads prior to this.  

-trim 0 : We're not trimming any data

-baq 1 : estimate base alignment quality using samtools method.

###ALLELE FREQUENCY ESTIMATION

-doMajorMinor 4 : Force Major allele based on reference. The minor allele is then inferred using doMajorMinor 1. This option needs to be used when calculating SFS for multiple populations as ANGSD otherwise determines a minor allele within each population. I.e. this may not be the same across all the populations.

-ref [..fasta] : For doMM 4 above we need to specify a reference genome.

-doMaf 1 : calculate minor allele frequency

-SNP_pval 0.001 : Only work with SNPs with a p-value above [float]

-GL 1 : I will estimate genotype likelihoods using the SAMtools model

-minInd 18 : I will remove loci where less than 18 individuals have been genotyped. There are 35-40 indivs per group (MUS, MODC, MODE). This number was chosen to ensure we can find MAF of 2-5% (1/(18*2)=2.8%)

-setMinDepth : Discard site if total depth (across all indivs) is below [int]. Use -doCounts to determine the distribution of depths

-setMaxDepth : Discard site if total depth (across all indivs) is above [int]

-setMinDepthInd 2 : Minimum depth for a locus for an individual. This is only applicable for analyses using counts (-doCounts obligatory)

-setMaxDepthInd [int]: I'll use a max of meanDP + 2xSD of depth.

-doCounts 1 : Count of the nucleotide bases per individual

-dumpCounts 2 : write a file with all the allele counts per position per individual

-rmTriallelic 1 : include only biallelic loci

-checkBamHeaders 1 : check that the bam headers are compatable for all files.



#### 3a. ANGSD filters for SFS

It's important to include all loci in the SFS. Don't include any filters that'll affect the allele frequencies (e.g. MAF or SNP_pval filters). 

It's also important to polarise the Major/Minor alleles so that they're the same between the populations. 

This is done with -doMajorMinor 4 : Force Major allele based on reference. The minor allele is then inferred using doMajorMinor 1 (i.e. from GL). This option needs to be used when calculating SFS for multiple populations as ANGSD otherwise determines a minor allele within each population. I.e. this may not be the same across all the populations.

The order of the filters do not seem to affect the outcome (see https://github.com/alexjvr1/Velocity2020/blob/master/04a_ANGSD_testFilters.md)

To add depth filters we need to add -doCount 1. 

NBNB for MUS data we have to select only_proper_pairs 0 - this is because the processing with MapDamage changes the reads in some way so that they all get filtered out if we select only_proper_pairs 1. They should already be filtered for proper_pairs based on prior mapping filters. 

The [basic process is](https://github.com/ANGSD/angsd/issues/259): 

###### 1. Estimate SAF for each population (unfolded) 

MUS
```
/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04a_ANGSD_TESTS/04_SFS/04a_MUS.SAF.sh

/newhome/aj18951/bin/angsd/angsd -b MUS.poplist -checkBamHeaders 1 -minQ 20 -minMapQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -r LR761675.1: -GL 1 -doSaf 1 -anc ../RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna -ref ../RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna -doCounts 1 -setMinDepthInd 2 -setMaxDepth 144 -doMajorMinor 4 -out MUS -C 50 -baq 1 

-> Total number of sites analyzed: 4786628
-> Number of sites retained after filtering: 4762287 
```

MODC
```
/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04a_ANGSD_TESTS/04_SFS/04a_MODC.SAF.sh

/newhome/aj18951/bin/angsd/angsd -b MODC.poplist -checkBamHeaders 1 -minQ 20 -minMapQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -r LR761675.1: -GL 1 -doSaf 1 -anc ../RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna -ref ../RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna -doCounts 1 -setMinDepthInd 2 -setMaxDepth 328 -doMajorMinor 4 -out MODC -C 50 -baq 1 

-> Total number of sites analyzed: 5715589
-> Number of sites retained after filtering: 5621808 

```

MODE
```
/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04a_ANGSD_TESTS/04_SFS/04a_MODE.SAF.sh

/newhome/aj18951/bin/angsd/angsd -b MODE.poplist -checkBamHeaders 1 -minQ 20 -minMapQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -r LR761675.1: -GL 1 -doSaf 1 -anc ../RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna -ref ../RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna -doCounts 1 -setMinDepthInd 2 -setMaxDepth 621 -doMajorMinor 4 -out MODE -C 50 -baq 1 

-> Total number of sites analyzed: 5688111
-> Number of sites retained after filtering: 5553623 
```


###### 2. unfolded SAF used to produce folded 2D SFS
```
realSFS pop1.unfolded.saf.idx pop2.unfolded.saf.idx -fold 1 >folded.sfs

```

3. unfolded SAF and folded SFS used to generate per-site numerator/denominator of Fst (use -fold in realSFS fst index)
```
realSFS fst index pop1.unfolded.saf.idx pop2.unfolded.saf.idx -sfs folded.sfs -fold 1 -fstout persite

```

4. sum numerator and denominator in windows
```
realSFS fst stat2 persite.fst.idx -win XXXX -step XXXX >window.fst
```


GlobalFST 

1. MODC vs MODE
```
~/bin/angsd/misc/realSFS fst stats MODE.MODC.persite.fst.idx

-> FST.Unweight[nObs:5498007]:0.020322 Fst.Weight:0.126046
0.020322	0.126046
```

2. MODC vs MUS
```
~/bin/angsd/misc/realSFS fst stats MUS.MODC.fstout.fst.idx 
	-> Assuming idxname:MUS.MODC.fstout.fst.idx
	-> Assuming .fst.gz file: MUS.MODC.fstout.fst.gz
	-> FST.Unweight[nObs:4738713]:0.091403 Fst.Weight:0.234494
0.091403	0.234494
```



#### 3b. Call GL

Test dataset: 
GL are called using the samtools method (GL1), and we're only allowing SNPs with a p-value of 0.001 for modern samples, and 0.05 for museum samples. 


Total sites analyzed and kept for each dataset: 

MOD Core p0.001 per SNP
```
/newhome/aj18951/E3_Aphantopus_hyperantus_2020/03a_ANGSD_SFS/logfiles

grep "Total number of sites analyzed" E3_MODC.ANGSD1.ARRAY.o9325111-* | awk -F "analyzed:" '{s+=$2} END {print s}'
381121498

grep "retained" E3_MODC.ANGSD1.ARRAY.o9325111-* | awk -F "filtering:" '{s+=$2} END {print s}'
6645193
```

MOD Exp p0.001 per SNP
```
/newhome/aj18951/E3_Aphantopus_hyperantus_2020/03a_ANGSD_SFS/logfiles

grep "Total number of sites analyzed" E3_MODE.ANGSD1.ARRAY.o9325110-* | awk -F "analyzed:" '{s+=$2} END {print s}'
379718544

grep "retained" E3_MODE.ANGSD1.ARRAY.o9325110-* | awk -F "filtering:" '{s+=$2} END {print s}'
5148017
```

*p0.05 per SNP*

MOD Core p0.05 per SNP
```
/newhome/aj18951/E3_Aphantopus_hyperantus_2020/03a_ANGSD_SFS/logfiles

grep "Total number of sites analyzed" E3_MODC.ANGSD1.ARRAY.o9541530-* | awk -F "analyzed:" '{s+=$2} END {print s}'
381121498 ##OLD p0.001
381142712

grep "retained" E3_MODC.ANGSD1.ARRAY.o9541530-* | awk -F "filtering:" '{s+=$2} END {print s}'
6645193 ##OLD p0.001
7529770
```

MOD Exp p0.05 per SNP
```
/newhome/aj18951/E3_Aphantopus_hyperantus_2020/03a_ANGSD_SFS/logfiles

grep "Total number of sites analyzed" E3_MODE.ANGSD1.ARRAY.o9541531-* | awk -F "analyzed:" '{s+=$2} END {print s}'
379718544 ##old p0.001
379785246

grep "retained" E3_MODE.ANGSD1.ARRAY.o9541531-* |awk -F "filtering:" '{s+=$2} END {print s}'
5148017 ##old p0.001
6265127
```


MUS
```
/newhome/aj18951/E3_Aphantopus_hyperantus_2020/03a_ANGSD_SFS/logfiles

grep "Total number of sites analyzed" E3_ANGSD1.ARRAY.o9403374* | awk -F "analyzed:" '{s+=$2} END {print s}'
346540624

grep "retained" E3_ANGSD1.ARRAY.o9403374* | awk -F "filtering:" '{s+=$2} END {print s}'
23542320
```



ANGSD is performed across all indivs for each contig/scaffold seperately. I then need to merge all of the data into a single file. 

```
module load languages/gcc-6.1

~/bin/angsd/misc/realSFS cat BA.HOD.*idx -outnames MUS.MERGED
```


#### 4d. Compare downsampled data to full dataset for modern pops (SFS and GL)


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

The developers have fixed a bug in the program. In addition I need to change the -only_proper_pairs flag to 0 as there are improperly paired reads in the merged file (i.e. reads located on different scaffolds). 



### 4. ANALYSES

Example papers

[Crested Ibis (Feng et al, Current biology, 2019)](https://www.sciencedirect.com/science/article/pii/S0960982218316099)

[Alpine chipmunts (Bi et al, PLoS Genetics, 2019)](https://journals.plos.org/plosgenetics/article?rev=2&id=10.1371/journal.pgen.1008119)

- Also good example of mapdamage use

- PopGen

- outFlank

[Dryas Monkey (ven der Valk, MBE, 2020)](https://academic.oup.com/mbe/article/37/1/183/5570178)


[Darwins finches (Pachecho, GBE, 2020)](https://academic.oup.com/gbe/article/12/3/136/5735467)

- example of using ANGSD to call GL

- diversity stats

[Lynx (Lucena-Perez, MolEcol, 2020)](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.15366?casa_token=I6nUZVYLpjAAAAAA%3AjGkueGUyt5AP-RdrD2GeejFl-U2BXVlrXqmYLcck4lOfFy4Yb2jMv3r6ZxXZ9WhGYDcmVDM0oJ8l)

- Genotype calls

- Pop structure

- Diversity measures

[Killer whale ROH (BioRxiv, Andy Foote & Excoffier)](https://www.biorxiv.org/content/10.1101/2020.04.08.031344v1.abstract)

- ROH


#### 4a. Diversity stats

##### 1. Tajima's D and Watterson's theta

We can calculate Tajima's D and Watterson's theta as [diversity stats in ANGSD](http://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests) after estimating the SFS for each population. When using folded SFS (as we have here) only the TD and WT will be meaningful despite several stats being calculated by these scripts: 

```
##1. Calculate the SFS from SAF. -fold 1 for folded SFS when reference is unknown

~/bin/angsd/misc/realSFS NAME.saf.idx -fold 1 > NAME.saf.sfs

## or use the script below and input names in the command line. Where $1 = input saf, and $2= output file name: 

qsub 04a_ANGSD_realSFS_cmdlineInputs.sh -F "NAME.saf.idx NAME.saf.sfs"

##2. Estimate thetas from SFS

~/bin/angsd/misc/realSFS saf2theta NAME.saf.idx -sfs NAME.saf.sfs -outname NAME.THETA.theta.gz

##3. Estimate thetas in windows from the above global theta file

~/bin/angsd/misc/thetaStat do_stat NAME.theta.idx -win 50000 -step 10000 -outnames NAME.THETA.window.gz

#4. Estimate theta per chromosome

~/bin/angsd/misc/thetaStat do_stat NAME.theta.idx
```


###### Scripts: 

Step1: [04a_ANGSD_realSFS_cmdlineInputs.sh](https://github.com/alexjvr1/Velocity2020/blob/master/04a_ANGSD_realSFS_cmdlineInputs.sh)


##### 2. observed heterozygosity

Comparison of the number of heterozygous sites in historic vs current datasets

Calculate this in ANGSD

Global estimate
```
/newhome/aj18951/E3_Aphantopus_hyperantus_2020/03a_ANGSD_SFS/test_pval0.001_nomaf

~/bin/angsd/misc/realSFS MODC.MERGED.saf.idx > est.ml

module load languages/R-3.6.2-gcc9.1.0

#in R
a<-scan("est.ml")
(a[2]/sum(a))*100   ##calculates the percentage of heterozygote sites

0.015% for MODC
```

How do we do this per individual??? 






##### 3. nucleotide diversity in windows

ThetaD from ANGSD window-based approach

Fig 3a in Feng et al: Distribution of nucleotide diversity across chromosomes in old vs new in 5Mb windows using nuc.div function in pegas


##### 4. Runs of Homozygosity 

Runs of homozygosity: Dryas monkey MS


Convert output to plink genotype calls: http://www.popgen.dk/angsd/index.php/Plink

Only use indivs with mean depth >> 3x? 

Use Plink to look for ROH in sliding window blocks of 50SNPs with one SNP per xkb (50kb?)
"allow for a maximum of one heterozygous and five missing calls per window before we considered the ROH to be broken"



##### 5. IBD in windows across chromosomes

Fig 3b in Feng et al: Distribution of IBD across chromosomes

##### 6. Deleterious load

Deleterious load: See Feng et al. 2019 method



#### 4b. Population structure

NJ phylogeny (Feng et al)

PCA [using ANGSD](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.5231). 
Use PCAngsd for this. We need Beagle output from ANGSD

??Structure



#### 4c. Outliers, Map, and identify candidates in the area. 

PCAngsd as outlier approach
```


```



##### Fst in windows to find outlying regions in ANGSD

ANGSD generates a ML SFS from genotype likelihoods for each population. It then uses the SFS as a prior in an empirical Bayes approach to estimate the posterior probability for all possible allelel frequencies at each site. (Method: Fumagalli et al. 2013) 
```

```

Plot Time vs space Fst plots for each chromosome. I'm using the Ringlet chromosome names as in the genome MS. Each chromosome is shown in two plots - plot1=MODC vs MODE, plot2= MODC vs Mus. 

Modern graphs
```
p1 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761647.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("Mod Core vs Mod Expanding LR761647.1") + theme(axis.title.x = element_blank())

p2 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761648.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("Mod Core vs Mod Expanding LR761648.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())

p3 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761649.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761649.1") + xlim(0,18830000) + theme(axis.title.x = element_blank())

p4 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761650.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761650.1") + xlim(0,18830000) + theme(axis.title.x = element_blank())

p5 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761651.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761651.1") + xlim(0,18830000) + theme(axis.title.x = element_blank())


p6 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761652.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761652.1") + xlim(0,18830000) + theme(axis.title.x = element_blank())


p7 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761653.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761653.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())

p8 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761654.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761654.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())

p9 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761655.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761655.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())

p10 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761656.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761656.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())


p11 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761657.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761657.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())

p12 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761658.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761658.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())

p13 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761659.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761659.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())

p14 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761660.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761660.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())

p15 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761661.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761661.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())


p16 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761662.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761662.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())


p17 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761663.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761663.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())

p18 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761664.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761664.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())

p19 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761665.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761665.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())

p20 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761666.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761666.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())


p21 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761667.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761667.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())

p22 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761668.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761668.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())

p23 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761669.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761669.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())

p24 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761670.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761670.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())

p25 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761671.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761671.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())


p26 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761672.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761672.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())


p27 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761673.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761673.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())

p28 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761674.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761674.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())

p29 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761675.1"),], aes(x=midpoint, y=fst)) + geom_point() + ggtitle("LR761675.1") + xlim(0,18830000)+ theme(axis.title.x = element_blank())

```




###### Synteny between modern and museum samples. 

BLAST with FlyBase & Enrichment analysis using PANTHER  

See [this](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.15188?casa_token=X-WHMot7TDcAAAAA%3Abn7IwwiinA44JDoEU-yuVV3iLk4RkXwcCU1av3_hKRG1hgKDNaCzPHbrEGlRCBk5j8bMcIW6ynjT) example paper



#### 4d. LD analyses

ngsLD (https://github.com/fgvieira/ngsLD) and plot r2 estimates using fit_LDdecay.R

Can use GL as input

See pigeon paper 


#### 4e. 




