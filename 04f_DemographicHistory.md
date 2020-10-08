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
~/bin/angsd/misc/realSFS MODE/MODE.LR761675.1.minDP20.MinIND10.saf.idx -fold 1 > MODE/MODE.LR75.fold.sfs


##While this is running it will output information we need: 
-> nSites: 5204985   ##MODC
-> nSites: 2093813   ##MUS
-> nSites: 5261260   ##MODE

##This info can also be obtained directly from the SFS by summing all the numbers. 
##Cat 1 = invariant sites, so remove this from the SFS before starting the analysis. There should be n catagories, where n=number of indivs. 



##Sample sizes of our filtered pops
MODE - 33
MODC - 36
MUS - 24
```




Modify the stairway plot blueprint file for this run
```
/newhome/aj18951/software/stairway_plot_v2.1.1/E3

nseq = pop size x 2 (i.e. haploid sequences)
L = total number of observed nucleotides. 
mu: 2.9 x 10(-9) #Mutation rate as estimated in H.melpomone https://mallet.oeb.harvard.edu/publications/estimation-spontaneous-mutation-rate-heliconius-melpomene
year_per_generation: 1 #Ringlet is univoltine
```


Plot the xx.final.summary output
```
ex.final <- read.table("two-epoch_fold.final.summary", header=T)
summary(ex.final)   
 mutation_per_site    n_estimation theta_per_site_median theta_per_site_2.5.
 Min.   :1.000e-09   Min.   :200   Min.   :0.0003688     Min.   :0.0002509  
 1st Qu.:1.326e-05   1st Qu.:200   1st Qu.:0.0004841     1st Qu.:0.0003666  
 Median :3.944e-05   Median :200   Median :0.0012178     Median :0.0008005  
 Mean   :5.691e-05   Mean   :200   Mean   :0.0009660     Mean   :0.0006809  
 3rd Qu.:8.545e-05   3rd Qu.:200   3rd Qu.:0.0012214     3rd Qu.:0.0009808  
 Max.   :3.788e-04   Max.   :200   Max.   :0.0012282     Max.   :0.0010354  
 theta_per_site_97.5.      year          Ne_median        Ne_2.5.     
 Min.   :0.0004922    Min.   :     2   Min.   : 7683   Min.   : 5227  
 1st Qu.:0.0011219    1st Qu.: 26529   1st Qu.:10086   1st Qu.: 7636  
 Median :0.0019101    Median : 78880   Median :25370   Median :16677  
 Mean   :0.0017166    Mean   :113829   Mean   :20124   Mean   :14185  
 3rd Qu.:0.0022742    3rd Qu.:170893   3rd Qu.:25445   3rd Qu.:20433  
 Max.   :0.0026079    Max.   :757606   Max.   :25588   Max.   :21572  
    Ne_97.5.        Ne_12.5.        Ne_87.5.    
 Min.   :10255   Min.   : 7289   Min.   : 7890  
 1st Qu.:23372   1st Qu.: 7743   1st Qu.:20988  
 Median :39794   Median :22714   Median :27703  
 Mean   :35763   Mean   :17266   Mean   :25244  
 3rd Qu.:47378   3rd Qu.:23502   3rd Qu.:30951  
 Max.   :54332   Max.   :23697   Max.   :39149  



ggplot(ex.final, aes(x=log(year), y=Ne_median/1000))+geom_line()
```

