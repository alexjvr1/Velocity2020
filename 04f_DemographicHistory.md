# Demographic History

## Goal:

Estimate the change in population size over time in each of the populations. 

I'll use two approaches: 

1. Stairway plot [(Liu & Fu 2015)](https://www.nature.com/articles/ng.3254)

2. fastSimCoal


### 1. Stairway plot

The package readme and download can be found [here](https://sites.google.com/site/jpopgen/stairway-plot)

Requirements: java1.7 or higher

No installation needed - just unzip the download and load java

```
module load languages/java-jdk-11.0.3
```

Create the folded SFS for each population 

```
#On the server

/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04b_ANGSD_FINAL/SFS_and_Fst


~/bin/angsd/misc/realSFS MODC/MODC.LR761675.1.minDP20.MinIND10.saf.idx -fold 1 > MODC/MODC.LR75.fold.sfs
~/bin/angsd/misc/realSFS MUS/MUS.LR761675.1.minDP20.MinIND10.saf.idx -fold 1 > MUS/MUS.LR75.fold.sfs
```


Visualise the SFS in R using heatmaps
```

```


Modify the stairway plot blueprint file for this run
```
/newhome/aj18951/software/
```
