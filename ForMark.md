# Summary of ANGSD filtering

## Pre-processing 

Details of the pre-processing steps and scripts can be found in the [README.md](https://github.com/alexjvr1/Velocity2020/edit/master/README.md) file. 

1a. Trim adapter sequence using cutadapt

1b. Concatenate resequenced museum data (some individuals have been sequenced >1)

1c. Repair problems in museum PE data for data from 1.2. (BBrepair)

1d. Merge overlapping PE reads in museum data (BBmerge)

2. Map and process
2a. Map museum and modern data to Sanger genome

2b. Correct museum data for possible deamination (MapDamage -> output = corrected bam file)

## ANGSD

