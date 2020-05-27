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

-baq 1 : estimate base alignment quality using samtools method.

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

MODC.pos <- read.table("MODC.pos", header=T)
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
     






```



Distribution of depths: 
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

GlobalDepth Histogram 

![alt.txt][GDH]

[GDH]:(https://user-images.githubusercontent.com/12142475/83019981-bdfb1300-a01f-11ea-99a7-e19a53ea493e.png)



How many loci do we lose with a minDP filter? 
```


```

