# Test filters and overlap in loci across datasets

There's a big dropout in the number of loci kept in the datasets after applying ANGSD filters. 

This is particularly true for the modern data where only 2% of the data are kept! (although initial counts include invariant sites)

In the museum data ~63% of the loci are retained using the initial filters (inlcuding keeping SNPs with p-value 0.05)

### Aim:

1) How many loci are filtered with each option?

2) How does this affect the intersect between datasets? 

3) How much time does it add to the analysis? 


### Method: 

I'll use only the largest contig for the analysis. 
I'll start with the basic filters then add one each time. 


## 1. Subset the data to use only largest contig

We can use the [region](http://www.popgen.dk/angsd/index.php/Input#BAM_files) filter in ANGSD (-r) to run the analysis only on one contig. I've chosen the smallest contig: 

LR761675.1: 6196582 bp (6.1Mb)


## 2. Depth per individual after basic filtering:


needs:

```
-doCounts 1 -doDepth 1
```

To calculate depth per individual and for all individuals jointly




### Basic filter set

Filters we have to include: 

-remove_bads 1 : remove reads with 255 flag (not primary, failure and duplicate reads) (1=default)

-uniqueOnly 1 : remove reads with multiple best hits

-minMapQ 20 : PHRED 20. This should already be in place during the mapping.

-minQ 20 : PHRED 20 for individual base score.

-only_proper_pairs 0 : NBNB THIS flag is changed to 0 because some of the reads in my final mus files are not properly paired! ***OLD FILTER include only properly paired reads (default) and should already have been applied to the museum reads prior to this.

-trim 0 : We're not trimming any data

-baq 1 : estimate base alignment quality using samtools method. BAQ 1 = more stringent, but might remove too many loci. BAQ 2 = extended BAQ which discovers more variants but include more false positives. See discussions [here](https://github.com/ANGSD/angsd/issues/97) and [here](https://github.com/ANGSD/angsd/issues/106).
For Ringlet data we don't lose very many loci when using BAQ 1, so I will opt for this more stringent approach (pending Het checks, as there has been at least one report of too many Het sites being removed with BAQ 1 in ancient DNA). 

###ALLELE FREQUENCY ESTIMATION

-doMajorMinor 4 : Force Major allele based on reference. The minor allele is then inferred using doMajorMinor 1. This option needs to be used when calculating SFS for multiple populations as ANGSD otherwise determines a minor allele within each population. I.e. this may not be the same across all the populations.

-ref [..fasta] : For doMM 4 above we need to specify a reference genome.



MUS:
```
~/bin/angsd/angsd -b MUS.poplist -checkBamHeaders 1 -minQ 20 -minMapQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -r LR761675.1:  -doCounts 1 -dumpCounts 2 -doDepth 1
	-> angsd version: 0.933-18-gfd1a21a (htslib: 1.10.2-61-g8859b09) build(May  6 2020 14:42:05)
	-> No '-out' argument given, output files will be called 'angsdput'
[bammer_main] 48 samples in 48 input files
	-> Parsing 48 number of samples 
	-> Region lookup 1/1
	-> Printing at chr: LR761675.1 pos:6174870 chunknumber 5900 contains 488 sitess
	-> Done reading data waiting for calculations to finish
	-> Done waiting for threads
	-> Output filenames:
		->"angsdput.arg"
		->"angsdput.pos.gz"
		->"angsdput.counts.gz"
		->"angsdput.depthSample"
		->"angsdput.depthGlobal"
	-> Tue May 26 15:14:56 2020
	-> Arguments and parameters for all analysis are located in .arg file
	-> Total number of sites analyzed: 4837675
	-> Number of sites retained after filtering: 4836300 
	[ALL done] cpu-time used =  78.55 sec
	[ALL done] walltime used =  189.00 sec

```




MOD.CORE
```
~/bin/angsd/angsd -b MODC.poplist -checkBamHeaders 1 -minQ 20 -minMapQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -r LR761675.1:  -doCounts 1 -dumpCounts 2 -doDepth 1 -out MODC

~/bin/angsd/angsd -b MODC.poplist -checkBamHeaders 1 -minQ 20 -minMapQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -r LR761675.1:  -doCounts 1 -dumpCounts 2 -doDepth 1 -out MODC
	-> angsd version: 0.933-18-gfd1a21a (htslib: 1.10.2-61-g8859b09) build(May  6 2020 14:42:05)
[bammer_main] 38 samples in 38 input files
	-> Parsing 38 number of samples 
	-> Region lookup 1/1

	-> Allocated ~ 10 million nodes to the nodepool, this is not an estimate of the memory usage
	-> Printing at chr: LR761675.1 pos:6018684 chunknumber 3000 contains 1624 sites
	-> Done reading data waiting for calculations to finish
	-> Done waiting for threads
	-> Output filenames:
		->"MODC.arg"
		->"MODC.pos.gz"
		->"MODC.counts.gz"
		->"MODC.depthSample"
		->"MODC.depthGlobal"
	-> Tue May 26 15:23:20 2020
	-> Arguments and parameters for all analysis are located in .arg file
	-> Total number of sites analyzed: 5721913
	-> Number of sites retained after filtering: 5719993 
	[ALL done] cpu-time used =  155.34 sec
	[ALL done] walltime used =  343.00 sec
```

MOD.EXP
```
~/bin/angsd/angsd -b MODE.poplist -checkBamHeaders 1 -minQ 20 -minMapQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -r LR761675.1:  -doCounts 1 -dumpCounts 2 -doDepth 1 -out MODE
	-> angsd version: 0.933-18-gfd1a21a (htslib: 1.10.2-61-g8859b09) build(May  6 2020 14:42:05)
[bammer_main] 40 samples in 40 input files
	-> Parsing 40 number of samples 
	-> Region lookup 1/1
	-> Printing at chr: LR761675.1 pos:1862919 chunknumber 2700 contains 496 sites
	-> Allocated ~ 10 million nodes to the nodepool, this is not an estimate of the memory usage

	-> Allocated ~ 20 million nodes to the nodepool, this is not an estimate of the memory usage

	-> Allocated ~ 30 million nodes to the nodepool, this is not an estimate of the memory usage
	-> Printing at chr: LR761675.1 pos:4345792 chunknumber 6400 contains 530 sitess


	-> Printing at chr: LR761675.1 pos:6128202 chunknumber 8900 contains 468 sitess
	-> Done reading data waiting for calculations to finish
	-> Done waiting for threads
	-> Output filenames:
		->"MODE.arg"
		->"MODE.pos.gz"
		->"MODE.counts.gz"
		->"MODE.depthSample"
		->"MODE.depthGlobal"
	-> Tue May 26 15:33:45 2020
	-> Arguments and parameters for all analysis are located in .arg file
	-> Total number of sites analyzed: 5701053
	-> Number of sites retained after filtering: 5697000 
	[ALL done] cpu-time used =  223.80 sec
	[ALL done] walltime used =  511.00 sec

```


#### Check overlap of these sites between datasets & assess depths 

copy everything to mac: 

```
/Users/alexjvr/2018.postdoc/Velocity2020/E3/Test.ANGDSstats/DEPTH

scp bluecp3:/newhome/aj18951/E3*/04*/DEPTH*/* .

##unzip the pos data

gunzip *pos.gz

#print the first two columns (i.e. not depth) to a second file for each dataset so that we can compare them

for i in $(ls *pos); do awk -F "\t" '{print $1,"\t",$2}' $i >> $i2

##how many SNPs in each dataset? 

>wc -l *pos2
 
 5719994 MODC.pos2
 5697001 MODE.pos2
 4836301 MUS.pos2

#find the overlap
#comm compares two sorted files. column 1 prints lines only in file1, line2 = file2, line3 = overlap. Suppress columns with -12

comm -12 MODC.pos2 MODE.pos2 |wc -l
4,721,101      ###~83% of each dataset

comm -12 MODC.pos2 MUS.pos2 |wc -l
4,027,018      ###70% MODC, 83% MUS

comm -12 MODE.pos2 MUS.pos2 |wc -l
4,025,715      ###70% MODC, 83% MUS
```


Missingness per individual
```
##How many 0 coverage sites per sample?

MODC.pos <- read.table("MODC.counts.gz", header=T)
MODC.res <- colSums(MODC.pos==0)/nrow(MODC.pos)*100   ##percentage of loci that will drop out because of low or no coverage. Assuming 2x limit
MODC.res  ##38 indivs named 0-37 here
 ind0TotDepth  ind1TotDepth  ind2TotDepth  ind3TotDepth  ind4TotDepth 
     32.23030      27.61680      30.20782      28.81853      33.22506 
 ind5TotDepth  ind6TotDepth  ind7TotDepth  ind8TotDepth  ind9TotDepth 
     34.02079      32.48880      27.74065      27.25474      48.22835 
ind10TotDepth ind11TotDepth ind12TotDepth ind13TotDepth ind14TotDepth 
     26.31533      25.62216      25.86388      26.18867      27.44773 
ind15TotDepth ind16TotDepth ind17TotDepth ind18TotDepth ind19TotDepth 
     32.30359      29.77352      27.31704      21.92863      42.87257 
ind20TotDepth ind21TotDepth ind22TotDepth ind23TotDepth ind24TotDepth 
     27.73068      32.29497      31.19343      30.51126      26.86881 
ind25TotDepth ind26TotDepth ind27TotDepth ind28TotDepth ind29TotDepth 
     29.41853      34.12212      26.58886      34.01071      35.02628 
ind30TotDepth ind31TotDepth ind32TotDepth ind33TotDepth ind34TotDepth 
     29.51528      35.99779      24.00923      24.84856      25.33979 
ind35TotDepth ind36TotDepth ind37TotDepth 
     22.66097      22.71679      22.70361 
     
MUS.pos <- read.table("MUS.counts.gz", header=T)
MODE.pos <- read.table("MODE.counts.gz", header=T)

MODE.res <- colSums(MODE.pos==0)/nrow(MODE.pos)*100 
MUS.res <- colSums(MUS.pos==0)/nrow(MUS.pos)*100 

library(reshape2)
MODE.res.melt <- melt(MODE.res)
MODE.res.melt$pop <- "MODE"
MODC.res.melt <- melt(MODC.res)
MODC.res.melt$pop <- "MODC"
MUS.res.melt <- melt(MUS.res)
MUS.res.melt$pop <- "MUS"
indivDepth <- rbind(MODE.res.melt, MODC.res.melt, MUS.res.melt)
indivDepth$colour <- indivDepth$pop
indivDepth$colour <- gsub("MODE", "MOD", indivDepth$colour)
indivDepth$colour <- gsub("MODC", "MOD", indivDepth$colour)

library(ggplot2)
pdf("Missingness.PerIndiv.pdf")
ggplot(indivDepth, aes(y=value, group=pop, colour=pop)) + geom_boxplot() + ggtitle("Ringlet Missingness per individual") + ylab("Missingness (%)") + xlab("Populations")                              
dev.off()

```

![alt.txt][missingness]

[missingness]:https://user-images.githubusercontent.com/12142475/83029353-82b21180-a02a-11ea-9322-ae4179e8c2e0.png


We'll exclude the 6 individuals from MOD.exp that have not sequenced well by removing them from the poplist: 

```
AH-02-2019-74
AH-02-2019-75
AH-02-2019-77
AH-02-2019-78
AH-02-2019-79
AH-02-2019-80
```


Distribution of depths (still with all 40 MODE indivs): 
```
##In R
##/Users/alexjvr/2018.postdoc/Velocity2020/E3/Test.ANGDSstats

MODC.pos <- read.table("MODC.pos.gz", header=T)
MODE.pos <- read.table("MODE.pos.gz", header=T)
MUS.pos <- read.table("MUS.pos.gz", header=T)
pdf("GlobalDepthHistograms.pdf")
par(mfrow=c(3,1))
hist((log10(MODC.pos$totDepth)), main="Mod.Core n=38", ylab="frequency", xlab="")
hist((log10(MODE.pos$totDepth)), main="Mod.Exp n=40", ylab="frequency", xlab="")
hist((log10(MUS.pos$totDepth)), main="Museum n=48", ylab="frequency", xlab="log10 Depth across all loci")
dev.off()
```

totalDepth Histogram 

![alt.txt][GDH]

[GDH]:https://user-images.githubusercontent.com/12142475/83019981-bdfb1300-a01f-11ea-99a7-e19a53ea493e.png


We'll set the maxDP for each population as mean + 2xSD. This is changed in the bash scripts. 
```
summary(MODC.pos$totDepth)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    1.0    87.0   117.0   113.1   137.0  6052.0 
sd(MODC.pos$totDepth)
[1] 107.3653

##MAXDEPTH: 113.1+ 2*(107.3)= 328X

summary(MODE.pos$totDepth)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    1.0   144.0   196.0   186.1   221.0 13702.0 
sd(MODE.pos$totDepth)
[1] 217.3255

##MAXDEPTH: 186.1+ 2*(217.3)= 621X


summary(MUS.pos$totDepth)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1.00   24.00   50.00   52.22   74.00 5672.00 
sd(MUS.pos$totDepth)
[1] 45.88092

##MAXDEPTH: 52.2+ 2*(45.9)= 144X
```


###### Depth filter


How many loci do we lose with a Different filters? 
```
/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04a_ANGSD_TESTS

##################
      MODE 
##################
-> Total number of sites analyzed: 5701053
-> Number of sites retained after filtering: 5667963 (MAXDP 621X)
-> Number of sites retained after filtering: 5496763 (MINDP 3x)
-> Number of sites retained after filtering: 5059904 (MinInd 10) 
-> Number of sites retained after filtering: 4655316 (MinInd 18)

-> Number of sites retained after filtering: 4831669 (MAXDP, MINDP, MinInd 10, -C50) 

-> Number of sites retained after filtering: 42814 (MAXDP, MINDP, rmTriallelic, baq, C,  P-val 0.001)
-> Number of sites retained after filtering: 45450 (MAXDP, MINDP, rmTriallelic, baq, C,  P-val 0.01)
-> Number of sites retained after filtering: 48340 (MAXDP, MINDP, rmTriallelic, baq, C,  P-val 0.05)
-> Number of sites retained after filtering: 58816 (noMaxDP and minDP=2)


-> Number of sites retained after filtering: 52197 (P-val 0.01, C, No baq)

###The -C filter (adjust for excessive mismatches) is removing 25% of the loci. Why is that? 
###The baq filter (base alignment quality) is also removing a lot of loci
###What if we remove the SNP_pval filter? 






110k for the OLD dataset vs ~50k for NEW dataset. What is the difference??

/newhome/aj18951/bin/angsd/angsd 
-b MODE.poplist 
-checkBamHeaders 1 
-ref ../RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna 
-minQ 20 
-minMapQ 20 
-uniqueOnly 1 
-remove_bads 1 
-only_proper_pairs 0 (1 for original)
-r LR761675.1: 
-doCounts 1 
-dumpCounts 2 
-doDepth 1 
-setMaxDepth 621 (new ONLY)
-setMinDepthInd 2 (old ONLY)
-setMinDepthInd 3 (new ONLY)
-C 50 (new ONLY)
-doMajorMinor 4 
-GL 1 
-doMaf 1 
-doSaf 1 
-anc ../RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna 
-rmTriallelic 1 
-SNP_pval 0.05 
-baq 1 (new ONLY) 
-out MODE.forSFS.PVAL0.05 
-doGeno 8 (old ONLY)
-doPost 1 (old ONLY) 



##################
      MODC 
##################
-> Total number of sites analyzed: 5721913
-> Number of sites retained after filtering: 5680716 (MAXDP 328x)
-> Number of sites retained after filtering: 5537574 (MINDP 3x)
-> Number of sites retained after filtering: 4708426 (MinInd 10)
-> Number of sites retained after filtering: 2890825 (MinInd 18)

-> Number of sites retained after filtering: 4507703 (MAXDP, MINDP, MinInd 10, -C50)

-> Number of sites retained after filtering: 49867 (MAXDP, MINDP, rmTriallelic, baq, C,  P-val 0.001) 
-> Number of sites retained after filtering: 54202 (MAXDP, MINDP, rmTriallelic, baq, C,  P-val 0.01)

-> Number of sites retained after filtering: 80191 (Pval 0.01, no minDP)



##################
       MUS
##################
-> Total number of sites analyzed: 4837675
-> Number of sites retained after filtering: 4807320 (MAXDP 166X)
-> Number of sites retained after filtering: 4253545 (MINDP 3X)
-> Number of sites retained after filtering: 1190194 (MinInd 10)
-> Number of sites retained after filtering: 89424 (MinInd 18)

-> Number of sites retained after filtering: 254283 (MAXDP, MINDP, MinInd 10, -C50)

-> Number of sites retained after filtering: 22575 (MAXDP, MINDP, rmTriallelic, baq, C,  P-val 0.001) 
-> Number of sites retained after filtering: 28705 (MAXDP, MINDP, rmTriallelic, baq, C,  P-val 0.01) 
-> Number of sites retained after filtering: 33273 (MAXDP, MINDP, rmTriallelic, baq, C,  P-val 0.05)
```


###### Overlap between datasets: 

```
##Prepare outputs: 

##unzip the pos data

gunzip *pos.gz

#print the first two columns (i.e. not depth) to a second file for each dataset so that we can compare them

for i in $(ls *pos); do awk -F "\t" '{print $1,"\t",$2}' $i >> $i.2; done
for i in $(ls *pos.2); do sort $i >> $i.3; done

#find the overlap
#comm compares two sorted files. column 1 prints lines only in file1, line2 = file2, line3 = overlap. Suppress columns with -12
#e.g. 

comm -12 MODC.pos2 MODE.pos2 |wc -l


#############
MODC vs MODE
#############

comm -12 MODC.MAXDP.MINDP.pos.2.3 MODE.MAXDP.MINDP.pos.2.3 | wc -l    
5,423,881										##MAXDP and MINDP (98% of each dataset retained)
comm -12 MODC.MAXDP.MINDP.MININD10.pos.2.3 MODE.MAXDP.MINDP.MININD10.pos.2.3 | wc -l
4,476,830										##MAXDP and MINDP & MinIND 10 (90-95% retained)
comm -12 MODC.MAXDP.MINDP.MININD18.pos.2.3 MODE.MAXDP.MINDP.MININD18.pos.2.3 | wc -l							 
2,852,432										##MAXDP and MINDP & MinIND 18 (62-98%)								

#############
MODC vs MUS
#############

comm -12 MODC.MAXDP.MINDP.pos.2.3 MUS.MAXDP.MINDP.pos.2.3 | wc -l    
4,231,209										##MAXDP and MINDP (78 - 87% retained)
comm -12 MODC.MAXDP.MINDP.MININD10.pos.2.3 MUS.MAXDP.MINDP.MININD10.pos.2.3 | wc -l
244,875											##MAXDP and MINDP & MinIND 10 (5-22%)
comm -12 MODC.MAXDP.MINDP.MININD18.pos.2.3 MUS.MAXDP.MINDP.MININD18.pos.2.3 | wc -l						 
63,938											##MAXDP and MINDP & MinIND 18 (1-6%)					
```


