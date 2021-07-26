# ATLAS

[ATLAS](https://bitbucket.org/wegmannlab/atlas/wiki/Home) is a pipeline developed in Daniel Wegmann's lab to process raw bam files from ancient DNA to obtain accurate variant calls and estimates of genetic diversity. 
The pipeline takes into account post-mortem damage (PMD) and low sequence coverage, and recovers variatns more accurately than the state-of-the-art method: MapDamage + GATK. 

ATLAS can be used to estimate genetic diversity. We can also generate input for ANGSD using ATLAS. 


## Test-species: 

E3 Aphantopus hyperantus

I'm testing the raw data pre-processing and ATLAS pipeline on the UCL server (CS) 


## Pre-processing: 

1. Concat raw museum samples that were sequenced twice, and move all samples to a folder called raw_museum_FINAL

2. Remove adapter sequence using Cutadapt

3. Repair and Merge overlapping sequences using BBmap

4. Map to genome using BWA mem

5. Add RG, remove duplicates, Local realignment using GATK

6. Validate bam file using PicardTools ValidateSamFile

7. ATLAS: splitMerge

8. ATLAS: PMD

9. ATLAS: recalibrate

10. ATLAS: output ANGSD inputs

11. ATLAS: estimate genetic diversity in windows


What is the state of the art with poolseq data? Should we use these approaches in windows for our data? 







