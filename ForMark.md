# Summary of ANGSD filtering

Details of the pre-processing and ANGSD steps and scripts can be found in the [README.md](https://github.com/alexjvr1/Velocity2020/edit/master/README.md) file. 

I've summarised them here and linked specifically to the pages where I've investigated the effect of depth on ANGSD estimates, as well as the estimates of genetic diversity. 

Finally I link to the scripts used to generate the data I've sent Mark

## Pre-processing 

1a. Trim adapter sequence using cutadapt

1b. Concatenate raw museum data for samples that have been sequenced twice.

1c. Repair problems in museum paired end (PE) data for data from 1b. (BBrepair)

1d. Merge overlapping PE reads in museum data (BBmerge)

2. Map and process
2a. Map museum and modern data to Sanger genome using BWA mem

2b. Correct museum data for possible deamination (MapDamage -> output = corrected bam file)


## ANGSD

