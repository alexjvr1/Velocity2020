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




### Estimates from ANGSD to plot against depth


#### Depth from ANGSD

ANGSD outputs the global depth at a site with the -doCounts and -dumpCounts options:

Here I'm using contig CADXM010000001.1
```
/newhome/aj18951/bin/angsd/angsd -b MODE.poplist -checkBamHeaders 1 -minQ 20 -minMapQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -r CADCXM010000001.1 -GL 1 -doSaf 1 -anc /newhome/aj18951/E3_Aphantopus_hyperantus_2020/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna -ref /newhome/aj18951/E3_Aphantopus_hyperantus_2020/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna -doCounts 1 -setMinDepthInd 2 -setMaxDepth 621 -doMajorMinor 4 -out /newhome/aj18951/E3_Aphantopus_hyperantus_2020/04a_ANGSD_FINAL/SFS_and_Fst/MODE/MODE.CADCXM010000001.1 -C 50 -baq 1 -dumpCounts 2 -doDepth 1 -maxDepth 621 

```

Note that I filter for maxDP based on the population meanDP.

I've set a minDepth per individual at a locus (2x), but no minInd filter for this particular run. 



#### 1. Fst

ANGSD calculates Fst in windows. To correlate to depth per position I'll make the windows really small. 

```
#GLs and SAF have previously been estimated for the whole genome. 
#Data on bluecrystal:/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04a_ANGSD_FINAL/SFS_and_Fst/

##Load necessary modules to run in terminal
module load languages/gcc-6.1

#call realSFS from ANGSD to estimate 2D-SFS
~/bin/angsd/misc/realSFS MODC/MODC.CADCXM010000001.1.saf.idx MODE/MODE.CADCXM010000001.1.saf.idx > MODC.MODE.CDX01.ml

#prepare for window-based analysis by indexing
~/bin/angsd/misc/realSFS fst index MODC/MODC.CADCXM010000001.1.saf.idx MODE/MODE.CADCXM010000001.1.saf.idx -sfs MODC.MODE.CDX01.ml -fstout CDX01.here

#The global Fst is similar to that found for LR75 (unw:0.020322; w:0.126046)
~/bin/angsd/misc/realSFS fst stats CDX01.here.fst.idx 
	-> Assuming idxname:CDX01.here.fst.idx
	-> Assuming .fst.gz file: CDX01.here.fst.gz
	-> FST.Unweight[nObs:36019]:0.040754 Fst.Weight:0.106320
0.040754	0.106320


#create two window-based outputs to plot: 10bp and 1bp "windows", non-overlapping
~/bin/angsd/misc/realSFS fst stats2 CDX01.here.fst.idx -win 10 -step 10 > CDX01.window10step10.fst
~/bin/angsd/misc/realSFS fst stats2 CDX01.here.fst.idx -win 1 -step 1 > CDX01.window1step1.fst

```


Download depth estimated by ANGSD and 3pop Fst file to mac

```
/Users/alexjvr/2020.postdoc/Velocity/E3/ANGSD_FINAL/DepthEstimates/IndivDepth

scp bluecp3:/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04a_ANGSD_FINAL/SFS_and_Fst/MODE/MODE.CADCXM010000001.1.pos.gz .
scp bluecp3:/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04a_ANGSD_FINAL/SFS_and_Fst/MODC/MODC.CADCXM010000001.1.pos.gz .
scp bluecp3:/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04a_ANGSD_FINAL/SFS_and_Fst/*step1.fst .
```

Remember to add an extra header (fst) in to the .fst file. 

Draw plots in R
```
library(dplyr)
library(ggplot2)

CDX.fst <- read.table("CDX01.window10step10.fst", header=T)
colnames(CDX.fst)<- c("region", "chr", "pos", "Nsites", "fst") #change headers so that "pos" header is in all files
MODE <- read.table(gzfile("MODE.CADCXM010000001.1.pos.gz"), header=T)
MODC <- read.table(gzfile("MODC.CADCXM010000001.1.pos.gz"), header=T)

MODE$pop <- "MODE"
MODC$pop <- "MODC"
MODE.MODC <- bind_rows(MODE, MODC) ##use dplyr to join the two tables together
MODE.MODC2 <- left_join(MODE, MODC, by="pos", suffix=c(".E", ".C"))

MODC.MODE2.fst <- left_join(MODC.MODE2, CDX.fst, by="pos") ##join fst data. Remember to change the header n the CDX.fst file to match "pos"
MODC.MODE2.fst <- MODC.MODE2.fst[complete.cases(MODC.MODE2.fst),]


#Plot depth vs fst
ggplot(MODE.MODC.fst[complete.cases(MODE.MODC.fst),], aes(y=fst, x=totDepth, colour=pop))+geom_point()

#plot fst and depth across the contig
#First I'm plotting fst vs pos for each pop, and colouring by depth. 

ggplot(MODC.MODE2.fst, aes(x=pos, y=fst, col=totDepth.E))+geom_point()
par(new=TRUE)
ggplot(MODC.MODE2.fst, aes(x=pos, y=fst, col=totDepth.C))+geom_point()

```

![alt_txt][depth.fst.MODE]

[depth.fst.MODE]:https://user-images.githubusercontent.com/12142475/93209477-19c7e400-f756-11ea-8a75-e3f03c28df40.png


![alt_txt][fstvsdepth.MODC]

[fstvsdepth.MODC]:https://user-images.githubusercontent.com/12142475/93209474-17fe2080-f756-11ea-9e5a-1feb119f71c9.png


We can see that the really high Fst (>0.6) is associated with low coverage. Let's look at these loci: 
```
MODC.MODE2.fst[which(MODC.MODE2.fst$fst>0.6),c("pos", "fst", "totDepth.E", "totDepth.C")]
        pos      fst totDepth.E totDepth.C
12823 78265 0.772495         18          9
13370 78875 0.780562         14         19
13380 78885 0.783590         13         18
```

And let's see if there's an overall correlation between depth and Fst: 
```
##plot depth vs fst coloured by pop
MODE.MODC.fst$col <- MODE.MODC.fst$pop
MODE.MODC.fst$col <- gsub("MODE", "blue",MODE.MODC.fst$col)
MODE.MODC.fst$col <- gsub("MODC", "green", MODE.MODC.fst$col)
ggplot(MODE.MODC.fst, aes(x=totDepth, y=fst, col=pop))+geom_point()

##And zoom in on the 0-150x
ggplot(MODE.MODC.fst[which(MODE.MODC.fst$totDepth<150),], aes(x=totDepth, y=fst, col=pop))+geom_point()
```

![alt_txt][DP.Fst]

[DP.Fst]:https://user-images.githubusercontent.com/12142475/93212532-9f4d9300-f75a-11ea-93e5-8de0879b15c7.png


![alt_txt][DP150]

[DP150]:https://user-images.githubusercontent.com/12142475/93212518-9c52a280-f75a-11ea-9241-610e1d6283e4.png



#### 2. nucleotide diversity




#### 2. GL



### "ESS" vs Fst

We're attempting to find a way to model the expected Fst so that we can identify possible outlier loci. Mark has written a script to estimate a score (0-1) for each locus based on the GLs which should be correlated with ESS (although it isn't ESS). 

This R script is [here]()

