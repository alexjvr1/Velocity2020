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
MUS.sampleDepth <- read.table("MUS.LR761675.1.depthSample", header=F)
colnames(MUS.sampleDepth) <- paste(rep(0:(ncol(MUS.sampleDepth)-1)), "x", sep="")
MUS.sampleDepth$Indiv <- paste("Indiv", rep(1:nrow(MUS.sampleDepth)), sep="")
mdata <- melt(MUS.sampleDepth, id="Indiv")

pdf("MUS.depthSample.LR761675.pdf")
ggplot(mdata[97:480,], aes(x=variable, y=value, group=Indiv, colour=Indiv)) + geom_line() + ggtitle("MUS LR761675 Distribution of sample depth") + xlab("Depth") + ylab("Count")  ##exclude 0x and 1x because of filters
dev.off()
```
The lack of 1x loci is due to the filter (minDP 2x), while the 0x loci is a reflection of the "gappiness" in the data - i.e. loci not sequenced in that indiv but that occur in at least one other indiv. 


![alt_txt][MODC.depth]

[MODC.depth]:https://user-images.githubusercontent.com/12142475/91590047-9be09c00-e952-11ea-9624-56776c72b231.png

![alt_txt][MODE.depth]

[MODE.depth]:https://user-images.githubusercontent.com/12142475/91590147-cc283a80-e952-11ea-815b-f13491446f17.png


![alt_txt][MUS.depth]

[MUS.depth]:


And for the global depth for each chromosome (LRxx) for each population. 

```


```

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







### Estimates from ANGSD to plot against depth

#### 1. Fst




#### 2. nucleotide diversity




#### 2. GL



