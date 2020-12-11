# Summary of ANGSD filtering

Details of the pre-processing and ANGSD steps and scripts can be found in the [README.md](https://github.com/alexjvr1/Velocity2020/edit/master/README.md) file. 

I've summarised them here and linked specifically to the pages where I've investigated the effect of depth on ANGSD estimates, as well as the estimates of genetic diversity. 

Finally I link to the scripts used to generate the data I've sent Mark

## Note on file naming convention

A typical file will be named as follows: 

```
MUS.LR761675.1.minDP20.MinIND10.HumanReadable.10k.saf.gz
```

pop.chromosome.filters.any additional info.file type

Files start with the population (MUS/MODC/MODE), followed by the chromosome name (LR7616xx.1). We're initially working with LR761675.1 which is the shortest chromosome (6.2Mb). 




## 1. Pre-processing 

1a. Trim adapter sequence using cutadapt

1b. Concatenate raw museum data for samples that have been sequenced twice.

1c. Repair problems in museum paired end (PE) data for data from 1b. (BBrepair)

1d. Merge overlapping PE reads in museum data (BBmerge)

2. Map and process
2a. Map museum and modern data to Sanger genome using BWA mem

2b. Correct museum data for possible deamination (MapDamage -> output = corrected bam file)


## 2. ANGSD


## 3. 

