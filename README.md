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

../../software/rename-master/rename 's/^[0-9]+-//' *

## Pipeline

### 1. Demultiplex and Adapter trimming

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

#####*pipeline

01a_museum_cutadapt_filtering_trimming.sh

01a_modern_cutadapt_filtering_trimming.sh

01b_museum_trimmomatic_filtering_trimming.sh

01b_modern_trimmomatic_filtering_trimming.sh

#####*wrapper

01a_parallel_cutadapt_bluecp3.sh

Edit the generated script above to submit from your home directory:

```
1. Set all paths to your home directory if necessary. 

2. Adjust the number of threads (PBS -t 1-xx) to equal the number of individuals to be analysed. 

3. Check that any empty arguments have been removed from the cutadapt command

4. You might have to set the path to cutadapt to find your local version

```

To incorporate new data (e.g. resequencing of some individuals to increase mean depth), new fastq files need to be adapter trimmed. Fastq files are concatenated after this using the script [concat.fastq.sh]()

Reseq data are kept in the following folders:
```
00_raw_data_museum2

01a_museum2_cutadapt_reads

01a_mus.concat_cutadapt_reads  ## concatenated museum1 and museum2 + all samples that didn't have reseq data added. I'll point to this folder when mapping

02a_museum2_mapped  ##see below
```

### 2. Map to Reference Genome

### 3. MapDamage: correct for Cytosine deamination in museum data

#### 3a. MapDamage run on museum data

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





