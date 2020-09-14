# Effect of sequencing depth on various ANGSD estimates

We expect that sequencing depth is correlated with the final Genotype likelihood scores, and hence the Fst and other estimates calculated in ANGSD. 

Here I explore these relationships. 

### Estimating depth 

#### ANGSD

Sequencing depth is estimated in ANGSD as a [distribution of depths](http://www.popgen.dk/angsd/index.php/Depth) either per sample or overall (global)

This is useful to get an idea of the mean and variance of the sequencing depth. e.g. here I plot depth for LR761675.1 for each indiv in MUS, MODC, and MODE. 
```
library(ggplot2)
library(reshape2)

##MODC 38 indivs
MODC.sampleDepth <- read.table("MODC.LR761675.1.depthSample", header=F)
colnames(MODC.sampleDepth) <- paste(rep(0:(ncol(MODC.sampleDepth)-1)), "x", sep="") ##ncol (330) -1
MODC.sampleDepth$Indiv <- paste("Indiv", rep(1:38), sep="")
mdata <- melt(MODC.sampleDepth, id="Indiv")

pdf("MODC.depthSample.LR761675.pdf")
ggplot(mdata[77:570,], aes(x=variable, y=value, group=Indiv, colour=Indiv)) + geom_line() + ggtitle("MODC LR761675 Distribution of sample depth") + xlab("Depth") + ylab("Count")  ##exclude 0x and 1x because of filters
dev.off()


##MODE 40 indivs
MODE.sampleDepth <- read.table("MODE.LR761675.1.depthSample", header=F)
colnames(MODE.sampleDepth) <- paste(rep(0:(ncol(MODE.sampleDepth)-1)), "x", sep="")
MODE.sampleDepth$Indiv <- paste("Indiv", rep(1:40), sep="")
mdata <- melt(MODE.sampleDepth, id="Indiv")

pdf("MODE.depthSample.LR761675.pdf")
ggplot(mdata[81:800,], aes(x=variable, y=value, group=Indiv, colour=Indiv)) + geom_line() + ggtitle("MODE LR761675 Distribution of sample depth") + xlab("Depth") + ylab("Count")  ##exclude 0x and 1x because of filters
dev.off()


##MUS 48 indivs
## The museum sequences are short and all paired, so the coverage is almost exclusively in multiples of 2x
MUS.sampleDepth <- read.table("MUS.LR761675.1.depthSample", header=F)
colnames(MUS.sampleDepth) <- paste(rep(0:(ncol(MUS.sampleDepth)-1)), "x", sep="")
MUS.sampleDepth$Indiv <- paste("Indiv", rep(1:nrow(MUS.sampleDepth)), sep="")
mdata <- melt(MUS.sampleDepth, id="Indiv")
mdata$number <- gsub("x","",mdata$numeric)
mdata.even <- mdata[mdata$number %% 2==0,] ##keep only rows with even coverage

pdf("MUS.depthSample.LR761675.pdf")
ggplot(mdata.even[49:480,], aes(x=variable, y=value, group=Indiv, colour=Indiv)) + geom_line() + ggtitle("MUS LR761675 Distribution of sample depth") + xlab("Depth") + ylab("Count")  ##exclude 0x and 1x because of filters
dev.off()
```
The lack of 1x loci is due to the filter (minDP 2x), while the 0x loci is a reflection of the "gappiness" in the data - i.e. loci not sequenced in that indiv but that occur in at least one other indiv. 


![alt_txt][MODC.depth]

[MODC.depth]:https://user-images.githubusercontent.com/12142475/91590047-9be09c00-e952-11ea-9624-56776c72b231.png

![alt_txt][MODE.depth]

[MODE.depth]:https://user-images.githubusercontent.com/12142475/91590147-cc283a80-e952-11ea-815b-f13491446f17.png


![alt_txt][MUS.depth]

[MUS.depth]:https://user-images.githubusercontent.com/12142475/92608341-214a4180-f2ad-11ea-8470-3522dfea6b99.png



And for the global depth for each chromosome (LRxx) for each population. 

```
#Copy Global Depth for each population: 

pwd
/Users/alexjvr/2020.postdoc/Velocity/E3/ANGSD_FINAL/DepthEstimates

##MODC

scp aj18951@bluecrystalp3.acrc.bris.ac.uk:/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04a_ANGSD_FINAL/SFS_and_Fst/MODC/*LR*depthGlobal .

#Read into R
library(ggplot2)
library(reshape2)
#The multiplot function from here: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)

fnames <- list.files(pattern="MODC.*depthGlobal")
tables <- lapply(fnames, read.table, header=F)
globalDepth.list <- do.call(rbind, tables) #paste each table below each other, so here we'll have 29 rows. 
colnames(globalDepth.list) <- paste(0:328, "x", sep="") #coverage starts at 0x
t.globalDepth <- t(globalDepth.list) #transpose to calculate the proportion of reads at depth x
freqs <- scale(t.globalDepth, center=F, scale=colSums(t.globalDepth)) #calculate proportions
colnames(freqs) <- paste("LR",47:75, sep="")
freqs.melt <- melt(freqs[1:200,])
freqs.melt$prop <- freqs.melt$value*100

#function to add breaks in the x-axis scale
every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}
p1.MODC.GlobalDepth <- ggplot(freqs.melt, aes(x=Var1, y=value, colour=Var2)) + geom_point()+ ggtitle("MODC Distribution of depth across chromosomes") + xlab("Depth") + ylab("Proportion of reads")+scale_x_discrete(breaks=every_nth(n=25))

pdf("MODC.GlobalDepth.perChr.pdf")
p1.MODC.GlobalDepth
dev.off()



##MODE

fnames <- list.files(pattern="MODE.*depthGlobal")
tables <- lapply(fnames, read.table, header=F)
globalDepth.list <- do.call(rbind, tables) #paste each table below each other, so here we'll have 29 rows. 
colnames(globalDepth.list) <- paste(0:621, "x", sep="") #coverage starts at 0x
t.globalDepth <- t(globalDepth.list) #transpose to calculate the proportion of reads at depth x
freqs <- scale(t.globalDepth, center=F, scale=colSums(t.globalDepth)) #calculate proportions
colnames(freqs) <- paste("LR",47:75, sep="")
freqs.melt <- melt(freqs)
freqs.melt$prop <- freqs.melt$value*100

every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}

p2.MODE.GlobalDepth <- ggplot(freqs.melt, aes(x=Var1, y=value, colour=Var2)) + geom_point()+ ggtitle("MODE Distribution of depth across chromosomes") + xlab("Depth") + ylab("Proportion of reads")+scale_x_discrete(breaks=every_nth(n=50))

pdf("MODE.GlobalDepth.perChr.pdf")
p2.MODE.GlobalDepth
dev.off()


##MUS

fnames <- list.files(pattern="MUS.*depthGlobal")
tables <- lapply(fnames, read.table, header=F)
globalDepth.list <- do.call(rbind, tables) #paste each table below each other, so here we'll have 29 rows. 
colnames(globalDepth.list) <- paste(0:144, "x", sep="") #coverage starts at 0x
t.globalDepth <- t(globalDepth.list) #transpose to calculate the proportion of reads at depth x
freqs <- scale(t.globalDepth, center=F, scale=colSums(t.globalDepth)) #calculate proportions
colnames(freqs) <- paste("LR",47:75, sep="")
freqs.melt <- melt(freqs)
freqs.melt$prop <- freqs.melt$value*100

every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}

p3.MUS.GlobalDepth <- ggplot(freqs.melt, aes(x=Var1, y=value, colour=Var2)) + geom_point()+ ggtitle("MUS Distribution of depth across chromosomes") + xlab("Depth") + ylab("Proportion of reads")+scale_x_discrete(breaks=every_nth(n=25))

pdf("MUS.GlobalDepth.perChr.pdf")
p3.MUS.GlobalDepth
dev.off()

multiplot(p1.MODC.GlobalDepth,p2.MODE.GlobalDepth,p3.MUS.GlobalDepth)
```


![alt_txt][Global.depth]

[Global.depth]:https://user-images.githubusercontent.com/12142475/92648675-b0714c80-f2e1-11ea-8d19-0e1dbe884a8d.png



The data used here is from the GLs estimated for the SFS on which the Fst calculations are based. 

i.e. Samtools GL model; derived allele frequencies polarised by the reference genome so that we can compare the three populations (-doMajorMinor 4). 
```
cat MODC.LR761675.1.arg
	-> Command: 
/newhome/aj18951/bin/angsd/angsd -b MODC.poplist -checkBamHeaders 1 -minQ 20 -minMapQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -r LR761675.1 -GL 1 -doSaf 1 -anc /newhome/aj18951/E3_Aphantopus_hyperantus_2020/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna -ref /newhome/aj18951/E3_Aphantopus_hyperantus_2020/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna -doCounts 1 -setMinDepthInd 2 -setMaxDepth 328 -doMajorMinor 4 -out /newhome/aj18951/E3_Aphantopus_hyperantus_2020/04a_ANGSD_FINAL/SFS_and_Fst/MODC/MODC.LR761675.1 -C 50 -baq 1 -doDepth 1 -maxDepth 328 -dumpCounts 2 
```

data located on bluecrystal 
```
/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04a_ANGSD_FINAL/SFS_and_Fst
```

And plots drawn on mac
```
/Users/alexjvr/2020.postdoc/Velocity/E3/ANGSD_FINAL/DepthEstimates
```


#### Samtools v.1.8

As it is simpler to plot individual mean depth vs Fst/GL/etc, I will also use depth estimates obtained directly from the bam files. 


I'm initially using only one chromosome: LR761675.1

Depth was estimated for LR..75 using [this](https://github.com/alexjvr1/Velocity2020/blob/master/LR75.depth.sh) script.

This calculates a depth for each individual at each locus. These files are very big as they contain every sequenced locus for an individual. 

Rename files using mv. 
```
##There are various versions of rename and mv on linux, so first check what is on your system
#rename --version
#mv --version
mv (GNU coreutils) 8.4
Copyright (C) 2010 Free Software Foundation, Inc.
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.

Written by Mike Parker, David MacKenzie, and Jim Meyering.


#I'm using mv to rename my files
#e.g. 

for i in $(ls *depth); do mv "$i" "${i/depth/depth.LR75}"; done 
#where ../xxx/yyy  is the text in the name (xxx), and the replacement text (yyy)
```


1. Combine these files together in the right order (i.e. ensure the bp positions line up and missing data is filled in as 0 or NA).

The size of the files means the easiest way to do this is with a database manager e.g. SQL. But if that isn't on the server the job can be split into batches of indivs and combined using awk or R. 

2. Extract data only from the loci for which we have GLs. 

3. Plot against Fst mid-point (for window based analysis), and for nucleotide diversity



### Estimates from ANGSD to plot against depth

#### 1. Fst




#### 2. nucleotide diversity




#### 2. GL



### "ESS" vs Fst

We're attempting to find a way to model the expected Fst so that we can identify possible outlier loci. Mark has written a script to estimate a score (0-1) for each locus based on the GLs which should be correlated with ESS (although it isn't ESS). 

This R script is [here]()

