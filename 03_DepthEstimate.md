# Effect of sequencing depth on various ANGSD estimates

We expect that sequencing depth is correlated with the final Genotype likelihood scores, and hence the Fst and other estimates calculated in ANGSD. 

Here I explore these relationships. 

### Estimating depth 

#### BedTools

Depth per bam file can be estimated with bedtools: 
```
module load apps/bedtools2
for i in $(ls *bam); do bedtools genomecov -ibam $i -g ../RefGenome/*fna > $i.depth

#And calculate the mean depth per chromosome for each individual
for i in $(ls *depth); do ls $i && awk -F "\t" '$1 ~ LR761675.1 {mult+=($2*$3); sum+=$3} END {print mult/sum}' $i; done
```



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

** Why does LR50 (Z) have higher coverage than the other chromosomes? There is a bit of variance between chromosomes, but the Z seems to have much higher mean depth than the other chromosomes. Could this have something to do with over-clustering due to mis-assembly of the reference Z-chr? 

LR50 is quite large, but isn't the largest chrs. Do we expect gene duplications here which could lead to mis-assembly of the ref? Or mis-mapping?


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

For when to use -fold see [here](https://github.com/ANGSD/angsd/issues/259). 

We need to use -fold when 1) calculating the 2Dsfs AND 2) the fst (fst index pop1.saf.idx pop2.saf.idx -sfs 2Dsfs.idx -fold 1 -fstout pop1.pop2.fold.fst)

Leaving -fold 1 out in the second part changes the final Fst. 

I'm using Bhatia's Fst ([Bhatia et al. 2013](https://pubmed.ncbi.nlm.nih.gov/23861382/); -whichFst 1) as suggested by the [ANGSD developers](https://github.com/ANGSD/angsd/issues/239), as this Fst is more robust to differences in sample size. 

```
#GLs and SAF have previously been estimated for the whole genome. 
#Data on bluecrystal:/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04a_ANGSD_FINAL/SFS_and_Fst/

##Load necessary modules to run in terminal
module load languages/gcc-6.1

#call realSFS from ANGSD to estimate 2D-SFS
#we need to calculate the FOLDED sfs and then specify folded SFS for the fst index. 

~/bin/angsd/misc/realSFS MODC/MODC.CADCXM010000001.1.saf.idx MODE/MODE.CADCXM010000001.1.saf.idx -fold 1 > MODC.MODE.CDX01.ml
~/bin/angsd/misc/realSFS MODC/MODC.CADCXM010000001.1.saf.idx MUS/MUS.CADCXM010000001.1.saf.idx -fold 1 > MODC.MUS.CDX01.ml
~/bin/angsd/misc/realSFS MUS/MUS.CADCXM010000001.1.saf.idx MODE/MODE.CADCXM010000001.1.saf.idx -fold 1 > MUS.MODE.CDX01.ml


#prepare for window-based analysis by indexing
~/bin/angsd/misc/realSFS fst index MODC/MODC.CADCXM010000001.1.saf.idx MODE/MODE.CADCXM010000001.1.saf.idx -sfs MODC.MODE.CDX01.ml -whichFst 1 -fold 1 -fstout CDX01.Bhatia.here


#Or for all three populations simultaneously
~/bin/angsd/misc/realSFS fst index MODC/MODC.CADCXM010000001.1.saf.idx MODE/MODE.CADCXM010000001.1.saf.idx MUS/MUS.CADCXM010000001.1.saf.idx -sfs MODC.MODE.CDX01.ml -sfs MODC.MUS.CDX01.ml -sfs MUS.MODE.CDX01.ml -whichFst 1 -fold 1 -fstout 3pop.CDX01.Bhatia.here

~/bin/angsd/misc/realSFS fst stats2 3pop.CDX01.Bhatia.here.fst.idx -win 50000 -step 10000 > 3pop.win50k.step10k.fst


#The global Fst is similar to that found for LR75 (unw:0.020322; w:0.126046)
~/bin/angsd/misc/realSFS fst stats CDX01.here.fst.idx 
	-> Assuming idxname:CDX01.Bhatia.here.fst.idx
	-> Assuming .fst.gz file: CDX01.Bhatia.here.fst.gz
	-> FST.Unweight[nObs:35958]:0.037287 Fst.Weight:0.108665
0.037287	0.108665


##Estimate Fst in windows: 
#create two window-based outputs to plot: 10bp and 1bp "windows", non-overlapping
~/bin/angsd/misc/realSFS fst stats2 CDX01.Bhatia.here.fst.idx -win 10 -step 10 > CDX01.window10step10.Bhatia.fst

##If you need to estimate Fst for e.g. a specific gene or part of the chromosome, print the per site fst
#then calc the weighted fst for the regions by calculating the ratio between the sum of As and the sum of B:
#Colnames are: Chr, pos, A = alpha from the reynolds 1983 (or Bhatia), B = alpha + beta

~/bin/angsd/misc/realSFS fst print CDX01.Bhatia.here.fst.idx > CDX01.Bhatia.here.fst.idx.HR
	-> Assuming idxname:CDX01.Bhatia.here.fst.idx
	-> Assuming .fst.gz file: CDX01.Bhatia.here.fst.gz
	-> Information from index file: nSites:35958 nChrs:1
	-> Population[0]: MODC/MODC.CADCXM010000001.1.saf.idx
	-> Population[1]: MODE/MODE.CADCXM010000001.1.saf.idx
		0	CADCXM010000001.1	35958	8
choose:1 
pars->stop:-1 ppos:146558 first:0 last:35958
CADCXM010000001.1	5917	0.000095	0.001535
CADCXM010000001.1	5918	0.000095	0.001536
CADCXM010000001.1	5919	0.000104	0.001603
CADCXM010000001.1	5920	0.000095	0.001536
CADCXM010000001.1	5921	0.000095	0.001535
CADCXM010000001.1	5922	0.000095	0.001535
CADCXM010000001.1	5923	0.000095	0.001536
CADCXM010000001.1	5924	0.000095	0.001535
CADCXM010000001.1	5925	0.000095	0.001535
CADCXM010000001.1	5926	0.000095	0.001535



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

CDX.fst <- read.table("CDX01.window10step10.Bhatia.fst", header=T)
colnames(CDX.fst)<- c("region", "chr", "pos", "Nsites", "fst") #change headers so that "pos" header is in all files
MODE <- read.table(gzfile("MODE.CADCXM010000001.1.pos.gz"), header=T)
MODC <- read.table(gzfile("MODC.CADCXM010000001.1.pos.gz"), header=T)

MODE$pop <- "MODE"
MODC$pop <- "MODC"
MODE.MODC <- bind_rows(MODE, MODC) ##use dplyr to join the two tables together
MODE.MODC2 <- inner_join(MODE, MODC, by="pos", suffix=c(".E", ".C"))

MODC.MODE2.fst <- inner_join(MODE.MODC2, CDX.fst, by="pos") ##join fst data. Remember to change the header n the CDX.fst file to match "pos"


#Plot depth vs fst
ggplot(MODC.MODE2.fst, aes(y=fst, x=totDepth, colour=pop))+geom_point()

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
1282 78265 0.718374         18          9
1337 78875 0.732677         14         19
1338 78885 0.739345         13         18
1746 84105 0.636810         11          8
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

##Correlation test
cor.test(MODC.MODE2.fst[which(MODC.MODE2.fst$totDepth.E>50),11], MODC.MODE2.fst[which(MODC.MODE2.fst$totDepth.E>50),3])

	Pearson's product-moment correlation

data:  MODC.MODE2.fst[which(MODC.MODE2.fst$totDepth.E > 50), 11] and MODC.MODE2.fst[which(MODC.MODE2.fst$totDepth.E > 50), 3]
t = -6.3473, df = 1528, p-value = 2.882e-10
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.2087162 -0.1110554
sample estimates:
      cor 
-0.160278 


cor.test(MODC.MODE2.fst[which(MODC.MODE2.fst$totDepth.C>50),11], MODC.MODE2.fst[which(MODC.MODE2.fst$totDepth.C>50),6])

	Pearson's product-moment correlation

data:  MODC.MODE2.fst[which(MODC.MODE2.fst$totDepth.C > 50), 11] and MODC.MODE2.fst[which(MODC.MODE2.fst$totDepth.C > 50), 6]
t = -4.6814, df = 849, p-value = 3.315e-06
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 -0.22345309 -0.09241244
sample estimates:
       cor 
-0.1586313 

##I find highly significant correlations between depth and Fst until ~150x within each population.
##Depth is negatively correlated with Fst, so we expect higher Fsts for low depth sites. 
##Check Fst peaks for depth distribution
#We expect that this will be smoothed out by the window-based approach

```

![alt_txt][DP.Fst]

[DP.Fst]:https://user-images.githubusercontent.com/12142475/93212532-9f4d9300-f75a-11ea-93e5-8de0879b15c7.png


![alt_txt][DP150]

[DP150]:https://user-images.githubusercontent.com/12142475/93212518-9c52a280-f75a-11ea-9241-610e1d6283e4.png



It looks like <50x produces Fst > 0.4. I need to rerun SFS estimates with: 

1) MinInd 18, and minIndDepth 2x set 

2) Problematic indivs removed


##### Fst vs Depth for filtered dataset: LR75

```
#On server: calculate fst per bp

~/bin/angsd/misc/realSFS fst stats2 MODC.MODE.LR75.minDP20.minIND10.fstout.fst.idx  -win 1 -step 2 > MODC.MODE.LR75.minDP20.minIND10.win1.step2.fst 
~/bin/angsd/misc/realSFS fst stats2 MODC.MUS.LR75.minDP20.minIND10.fstout.fst.idx -win 1 -step 2 > MODC.MUS.LR75.minDP20.minIND10.win1.step2.fst 
~/bin/angsd/misc/realSFS fst stats2 MODE.MUS.LR75.minDP20.minIND10.fstout.fst.idx -win 1 -step 2 > MODE.MUS.LR75.minDP20.minIND10.win1.step2.fst



##On Mac
/Users/alexjvr/2020.postdoc/Velocity/E3/ANGSD_FINAL/DepthEstimates/FstVsDepth
scp bluecp3:/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04b_ANGSD_FINAL/SFS_and_Fst/*step2.fst .


#In R
library(dplyr)
library(reshape2)
library(ggplot2)

#Read in all data
MODC.MODE.fst <- read.table("MODC.MODE.LR75.minDP20.minIND10.win1.step2.fst", header=T)
colnames(MODC.MODE.fst)[2]<- "pos" #change headers so that "pos" header is in all files
colnames(MODC.MODE.fst)[4]<- "fst"

MODC.MUS.fst <- read.table("MODC.MUS.LR75.minDP20.minIND10.win1.step2.fst", header=T)
colnames(MODC.MUS.fst)[2]<- "pos" #change headers so that "pos" header is in all files
colnames(MODC.MUS.fst)[4]<- "fst"

MODE.MUS.fst <- read.table("MODE.MUS.LR75.minDP20.minIND10.win1.step2.fst", header=T)
colnames(MODE.MUS.fst)[2]<- "pos" #change headers so that "pos" header is in all files
colnames(MODE.MUS.fst)[4]<- "fst"

MUS <- read.table(gzfile("../DiversityVsDepth/MUS.LR761675.1.minDP20.MinIND10.pos.gz"), header=T)
MODE <- read.table(gzfile("../DiversityVsDepth/MODE.LR761675.1.minDP20.MinIND10.pos.gz"), header=T)
MODC <- read.table(gzfile("../DiversityVsDepth/MODC.LR761675.1.minDP20.MinIND10.pos.gz"), header=T)

MODE$pop <- "MODE"
MODC$pop <- "MODC"
MUS$pop <- "MUS"

##MODC.MODE
#df1
MODE.MODC <- bind_rows(MODE, MODC) ##use dplyr to join the two tables together
MODE.MODC1.fst <- left_join(MODE.MODC, MODC.MODE.fst, by="pos")
MODE.MODC1.fst <- MODE.MODC1.fst[complete.cases(MODE.MODC1.fst),]

#df2
MODC.MODE2 <- left_join(MODC, MODE, by="pos", suffix=c(".C", ".E")) ##use dplyr to join the two tables together
MODC.MODE2.fst <- left_join(MODC.MODE2, MODC.MODE.fst, by="pos") ##join fst data. Remember to change the header n the CDX.fst file to match "pos"



#Plot depth vs fst
ggplot(MODE.MODC1.fst[which(MODE.MODC1.fst$pos<100000&MODE.MODC1.fst$totDepth<300),], aes(x=totDepth, y=fst, color=pop))+geom_point()

#Plot fst and depth across the contig
#First I'm plotting fst vs pos for each pop, and colouring by depth. 

ggplot(MODC.MODE2.fst[1:100000,], aes(x=pos, y=fst, col=totDepth.E))+geom_point()
par(new=TRUE)
ggplot(MODC.MODE2.fst[1:100000,], aes(x=pos, y=fst, col=totDepth.C))+geom_point()


##MODC.MUS
#df1
MODC.MUS <- bind_rows(MODC, MUS) ##use dplyr to join the two tables together
MODC.MUS1.fst <- left_join(MODC.MUS, MODC.MUS.fst, by="pos")
MODC.MUS1.fst <- MODC.MUS1.fst[complete.cases(MODC.MUS1.fst),]

#df2
MODC.MUS2 <- left_join(MODC, MUS, by="pos", suffix=c(".C", ".M")) ##use dplyr to join the two tables together
MODC.MUS2.fst <- left_join(MODC.MUS2, MODC.MUS.fst, by="pos") ##join fst data. Remember to change the header n the CDX.fst file to match "pos"



#Plot depth vs fst
ggplot(MODC.MUS1.fst[which(MODC.MUS1.fst$pos<100000&MODC.MUS1.fst$totDepth<100),], aes(x=totDepth, y=fst, color=pop))+geom_point()

#Plot fst and depth across the contig
#First I'm plotting fst vs pos for each pop, and colouring by depth. 

ggplot(MODC.MUS2.fst[1:100000,], aes(x=pos, y=fst, col=totDepth.C))+geom_point()
par(new=TRUE)
ggplot(MODC.MUS2.fst[1:100000,], aes(x=pos, y=fst, col=totDepth.M))+geom_point()

##MODE.MUS
#df1
MODE.MUS <- bind_rows(MODE, MUS) ##use dplyr to join the two tables together
MODE.MUS1.fst <- left_join(MODE.MUS, MODE.MUS.fst, by="pos")
MODE.MUS1.fst <- MODE.MUS1.fst[complete.cases(MODE.MUS1.fst),]

#df2
MODE.MUS2 <- left_join(MODE, MUS, by="pos", suffix=c(".E", ".M")) ##use dplyr to join the two tables together
MODE.MUS2.fst <- left_join(MODE.MUS2, MODE.MUS.fst, by="pos") ##join fst data. Remember to change the header n the CDX.fst file to match "pos"



#Plot depth vs fst
ggplot(MODE.MUS1.fst[which(MODE.MUS1.fst$pos<100000&MODE.MUS1.fst$totDepth<150),], aes(x=totDepth, y=fst, color=pop))+geom_point()

#Plot fst and depth across the contig
#First I'm plotting fst vs pos for each pop, and colouring by depth. 

ggplot(MODE.MUS2.fst[1:100000,], aes(x=pos, y=fst, col=totDepth.E))+geom_point()
par(new=TRUE)
ggplot(MODE.MUS2.fst[1:100000,], aes(x=pos, y=fst, col=totDepth.M))+geom_point()


#Adjust dataset to see if lower depth affects Fst
MODE.MUS2.fst.sub <- MODE.MUS2.fst[which(MODE.MUS2.fst$totDepth.E<150),]
ggplot(MODE.MUS2.fst.sub[1:100000,], aes(x=pos, y=fst, col=totDepth.E))+geom_point()

MODE.MUS2.fst.sub <- MODE.MUS2.fst[which(MODE.MUS2.fst$totDepth.M<50),]
ggplot(MODE.MUS2.fst.sub[1:100000,], aes(x=pos, y=fst, col=totDepth.M))+geom_point()

```

![alt_txt][fst.depth.plot1]

[fst.depth.plot1]:https://user-images.githubusercontent.com/12142475/94037206-08698200-fdbd-11ea-82eb-9ea3a9942abc.png

![alt_txt][fst.depth.plot2]

[fst.depth.plot2]:https://user-images.githubusercontent.com/12142475/94038070-f2a88c80-fdbd-11ea-823a-ac954e9ccf44.png

![alt_txt][fst.depth.plot3]

[fst.depth.plot3]:https://user-images.githubusercontent.com/12142475/94038089-f76d4080-fdbd-11ea-9f69-b79fb9e4ad32.png


![alt_txt][fst.depth.plot4]

[fst.depth.plot4]:https://user-images.githubusercontent.com/12142475/94038194-12d84b80-fdbe-11ea-98e9-9ece4d1f3fc9.png

![alt_txt][fst.depth.plot5]

[fst.depth.plot5]:https://user-images.githubusercontent.com/12142475/94038421-52069c80-fdbe-11ea-816a-32a2ee7a2200.png

![alt_txt][fst.depth.plot6]

[fst.depth.plot6]:https://user-images.githubusercontent.com/12142475/94038430-5468f680-fdbe-11ea-8627-d89661e5172e.png


![alt_txt][fst.depth.plot7]

[fst.depth.plot7]:https://user-images.githubusercontent.com/12142475/94038945-03a5cd80-fdbf-11ea-8cf3-9bd45dbca6e8.png

![alt_txt][fst.depth.plot8]

[fst.depth.plot8]:https://user-images.githubusercontent.com/12142475/94038956-056f9100-fdbf-11ea-8124-863e39c02193.png


##### Lower Depth plots

![alt_txt][fst.depth.plot9]

[fst.depth.plot9]:https://user-images.githubusercontent.com/12142475/94039898-34d2cd80-fdc0-11ea-8e0d-60b361bd52b3.png


#### 2. nucleotide diversity

Estimate thetas for each population (Still using CADX..01)

```
/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04a_ANGSD_OLD/SFS_and_Fst

module load languages/gcc-6.1

#create folded SFS
~/bin/angsd/misc/realSFS MODE/MODE.CADCXM010000001.1.saf.idx -fold 1 > MODE/MODE.CADCXM010000001.1.sfs
~/bin/angsd/misc/realSFS saf2theta MODE/MODE.CADCXM010000001.1.saf.idx -sfs MODE/MODE.CADCXM010000001.1.sfs -outname MODE.CADX01

#calculate thetas with 1bp non-overlapping windows
~/bin/angsd/misc/thetaStat do_stat MODE.CADX01.thetas.idx -win 1 -step 1

#Information from index file:
		0	CADCXM010000001.1	43168	8	80

##same for MODC
~/bin/angsd/misc/realSFS MODC/MODC.CADCXM010000001.1.saf.idx -fold 1 > MODC/MODC.CADCXM010000001.1.sfs
~/bin/angsd/misc/realSFS saf2theta MODC/MODC.CADCXM010000001.1.saf.idx -sfs MODC/MODC.CADCXM010000001.1.sfs -outname MODC.CADX01

#calculate thetas with 1bp non-overlapping windows
~/bin/angsd/misc/thetaStat do_stat MODC.CADX01.thetas.idx -win 1 -step 1

#Information from index file:
		0	CADCXM010000001.1	49178	8	76

#Copy to Mac

pwd
/Users/alexjvr/2020.postdoc/Velocity/E3/ANGSD_FINAL/SFS

scp bluecp3:/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04a_ANGSD_OLD/SFS_and_Fst/*pestPG .

##In R (see previous section for Fst tables)

MODC.theta <- read.table("../../SFS/MODC.CADX01.thetas.idx.pestPG", header=F)
colnames(MODC.theta) <- c("region", "chr", "pos", "tW","tP","tF","tH","tL", "Tajima", "fuf", "fud", "fayh", "zeng", "nSites")

MODE.theta <- read.table("../../SFS/MODE.CADX01.thetas.idx.pestPG", header=F)
colnames(MODE.theta) <- c("region", "chr", "pos", "tW","tP","tF","tH","tL", "Tajima", "fuf", "fud", "fayh", "zeng", "nSites")

MODE.MODC.theta <- left_join(MODE.theta, MODC.theta, by="pos", suffix=c(".E", ".C"))
MODE.MODC.theta <- MODE.MODC.theta[complete.cases(MODE.MODC.theta),]
dim(MODE.MODC.theta)
[1] 36017    27

MODE.MODC.theta.depth <- left_join(MODE.MODC.theta, MODC.MODE2, by="pos")
MODE.MODC.theta.depth <- MODE.MODC.theta.depth[complete.cases(MODE.MODC.theta.depth),]

dim(MODE.MODC.theta.depth)
[1] 36017    33

#Check if theta is correlated with depth
ggplot(MODE.MODC.theta.depth, aes(x=totDepth.C, y=tW.C))+ geom_point()
ggplot(MODE.MODC.theta.depth, aes(x=totDepth.E, y=tW.E))+ geom_point()

#This shows quite a strange pattern of 0 and 0.2 Wattersons theta. Why would it be so binary? 
#There is some indication that low Depth results in more noise in the estimates.

#Check if the high Wtheta is concentrated around the Fst peak: 
ggplot(MODE.MODC.theta.depth, aes(x=pos, y=tW.E))+ geom_point()
ggplot(MODE.MODC.theta.depth, aes(x=pos, y=tW.C))+ geom_point()
```

Left = MODC, Right = MODE


![alt_txt][depth.vsTheta]

[depth.vsTheta]:https://user-images.githubusercontent.com/12142475/93225095-c19add00-f769-11ea-9a6e-3208e27aea3f.png

Fig1: Wattersons theta vs depth


![alt_txt][pos.Theta]

[pos.Theta]:https://user-images.githubusercontent.com/12142475/93225087-bf388300-f769-11ea-8257-20ac24ff9fdf.png

Fig2. Wattersons theta across the contig. This shows the biggest variance in theta around the biggest Fst peak, which is also where we have a big variance in sequencing depth. 



##### Nucleotide Diversity with LR75


Copy pestPG for win1step1 (i.e. calculated for every bp) for LR75 from bluecrystal to mac

```
pwd
/Users/alexjvr/2020.postdoc/Velocity/E3/ANGSD_FINAL/DiversityEstimates

scp bluecp3:/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04b_ANGSD_FINAL/SFS_and_Fst/M*/*LR761675.1.minDP20.MinIND10.win1.step1.thetasWindow.gz.pestPG .

scp bluecp3:/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04b_ANGSD_FINAL/SFS_and_Fst/M*/*LR761675.1.minDP20.MinIND10.pos.gz .

##In R (see previous section for Fst tables)

library(dplyr)
library(reshape2)
library(ggplot2)

##MODC
MODC.theta <- read.table("MODC.LR761675.1.minDP20.MinIND10.win1.step1.thetasWindow.gz.pestPG", header=F)
colnames(MODC.theta) <- c("region", "chr", "pos", "tW","tP","tF","tH","tL", "Tajima", "fuf", "fud", "fayh", "zeng", "nSites")

MODC <- read.table(gzfile("MODC.LR761675.1.minDP20.MinIND10.pos.gz"), header=T)
MODC$pop <- "MODC"

MODC1 <- left_join(MODC.theta, MODC, by="pos")

##MODE
MODE.theta <- read.table("MODE.LR761675.1.minDP20.MinIND10.win1.step1.thetasWindow.gz.pestPG", header=F)
colnames(MODE.theta) <- c("region", "chr", "pos", "tW","tP","tF","tH","tL", "Tajima", "fuf", "fud", "fayh", "zeng", "nSites")

MODE <- read.table(gzfile("MODE.LR761675.1.minDP20.MinIND10.pos.gz"), header=T)
MODE$pop <- "MODE"

MODE1 <- left_join(MODE.theta, MODE, by="pos")


##MUS
MUS.theta <- read.table("MUS.LR761675.1.minDP20.MinIND10.win1.step1.thetasWindow.gz.pestPG", header=F)
colnames(MUS.theta) <- c("region", "chr", "pos", "tW","tP","tF","tH","tL", "Tajima", "fuf", "fud", "fayh", "zeng", "nSites")

MUS <- read.table(gzfile("MUS.LR761675.1.minDP20.MinIND10.pos.gz"), header=T)
MUS$pop <- "MUS"

MUS1 <- left_join(MUS.theta, MUS, by="pos")

##plots
#Check if theta is correlated with depth
ggplot(MODC[10000,], aes(x=totDepth, y=tW))+ geom_point() + ggtitle("MODC, LR75")
ggplot(MODE[10000,], aes(x=totDepth, y=tW))+ geom_point()+ gg_title("MODE, LR75")
ggplot(MUS[10000,], aes(x=totDepth, y=tW))+ geom_point()+ gg_title("MUS, LR75")

#This shows quite a strange pattern of 0 and 0.2 Wattersons theta. Why would it be so binary? 
#There is some indication that low Depth results in more noise in the estimates.

#Check if the high Wtheta is concentrated around the Fst peak: 
ggplot(MODE.MODC.theta.depth, aes(x=pos, y=tW.E))+ geom_point()
ggplot(MODE.MODC.theta.depth, aes(x=pos, y=tW.C))+ geom_point()

```



#### Nucleotide Diversity

For LR75

The log scaled diversity estimates can be obtained using ANGSD

```
/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04b_ANGSD_FINAL/SFS_and_Fst

module load languages/gcc-6.1

~/bin/angsd/misc/thetaStat print MODC/MODC.LR761675.1.minDP20.MinIND10.thetas.idx > MODC/MODC.LR761675.1.minDP20.minIND10.thetas.HumanReadable

~/bin/angsd/misc/thetaStat print MODE/MODE.LR761675.1.minDP20.MinIND10.thetas.idx > MODE/MODE.LR761675.1.minDP20.minIND10.thetas.HumanReadable

~/bin/angsd/misc/thetaStat print MUS/MUS.LR761675.1.minDP20.MinIND10.thetas.idx > MUS/MUS.LR761675.1.minDP20.minIND10.thetas.HumanReadable



```

Draw the same plots for Nucleotide diversity as for Watterson's theta. See scripts [here](https://github.com/alexjvr1/Velocity2020/blob/master/NucelotideDiversityPlot.md) 





## 2. GL

### "ESS" vs Fst

We're attempting to find a way to model the expected Fst so that we can identify possible outlier loci. Mark has written a script to estimate a score (0-1) for each locus based on the GLs which should be correlated with ESS (although it isn't ESS). 

This R script is [here](https://github.com/alexjvr1/Velocity2020/blob/master/ESSforGLs.R)

First we have to write the GL file for CADX..01 for the same filter options used throughout this doc. 

Rerun the script for CADX..01 adding the beagle GL (-doGlf 2) output option: 

```
qsub 04a_MODE.SAF_withGL.sh 

#!/bin/bash
#PBS -N E3_MODE.SAF  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=20:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)
##PBS -t 1-87 #array job

#Set filters
MININD=""
MINMAF=""
#PVAL="0.05"
MINQ="20"
minMAPQ="20"
minDP="2"
maxDP="621"
POP="MODE"
C="50"
POP="MODE"
POPLIST="MODE.poplist"
SPECIESDIR="/newhome/aj18951/E3_Aphantopus_hyperantus_2020"
PP=1 #use all reads. Flag 1 uses only proper pairs, but	MUS has	merged reads. NB to filter for proper pair reads in the bamfiles using samtools before this point

#OUTNAME="MODE"

#run job in working directory
cd $PBS_O_WORKDIR 

#load modules
module load languages/gcc-6.1
angsd=~/bin/angsd/angsd

#Define variables
#REGION=$(sed "${PBS_ARRAYID}q;d" regions)
REGION="CADCXM010000001.1"

#estimate SAF for modern expanding population using ANGSD

time $angsd -b $POP.poplist -checkBamHeaders 1 -minQ 20 -minMapQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs $PP -r $REGION \
-GL 1 -doSaf 1 -anc $SPECIESDIR/RefGenome/*fna -ref $SPECIESDIR/RefGenome/*fna -doCounts 1 -setMinDepthInd 2 -setMaxDepth $maxDP -doMajorMinor 4\
 -out $SPECIESDIR/04a_ANGSD_FINAL/SFS_and_Fst/$POP/$POP.$REGION.new -C $C -baq 1 -dumpCounts 2 -doDepth 1 -maxDepth $maxDP -doGlf 2
```


And plot in R on the mac: 
```
/Users/alexjvr/2020.postdoc/Velocity/E3/ANGSD_FINAL/SFS/Score.Mark

Rscript script2.R MODE.CADCXM010000001.1.minDP2.old.beagle.gz MODE.minDP2.old.scorevec
[1] "number of individuals is  40"
[1] "maximum 'N_eff' is  34.6517613256096"
[1] "minimum 'N_eff' is  0.599662013119725"


##R
library(dplyr)
library(ggplot2)

MODE.GL <- read.table(gzfile("../../SFS/MODE.CADCXM010000001.1.minDP2.old.beagle.gz"), header=T)

head(MODE.GL)
                  marker allele1 allele2     Ind0   Ind0.1   Ind0.2     Ind1
1 CADCXM010000001.1_5740       3       0 0.333333 0.333333 0.333333 0.333333
2 CADCXM010000001.1_5741       1       0 0.333333 0.333333 0.333333 0.333333
3 CADCXM010000001.1_5742       3       0 0.333333 0.333333 0.333333 0.333333
4 CADCXM010000001.1_5743       0       1 0.333333 0.333333 0.333333 0.333333
5 CADCXM010000001.1_5744       2       0 0.333333 0.333333 0.333333 0.333333
6 CADCXM010000001.1_5745       2       0 0.333333 0.333333 0.333333 0.333333

dim(MODE.GL)
[1] 43168   123


MODE.GL.names <- MODE.GL$marker
MODE.GL.names <- as.data.frame(MODE.GL.names)
MODE.GL.names$pos <- MODE.GL$marker
MODE.GL.names$pos <- gsub("CADCXM010000001.1_", "", MODE.GL.names$pos)
MODE.GL.names$pos <- as.numeric(MODE.GL.names$pos)


CDX.fst <- read.table("MODE.MODC.CDX01.window1step1.minDP2.old.fst", header=T) #read in Fst. Remember to add "fst" as header in the last column using nano
colnames(CDX.fst) <- c("index", "marker", "pos", "N", "fst")

MODE.score <- read.table("MODE.minDP2.old.scorevec", header=T)
MODE.score$pos <- MODE.GL.names$pos  
colnames(MODE.score) <- c("score", "pos")

MODE.GL.sub <- left_join(CDX.fst, MODE.score, by="pos")

ggplot(MODE.GL.sub, aes(y=fst, x=score))+ geom_point()
```

![alt_txt][score.fig]

[score.fig]:https://user-images.githubusercontent.com/12142475/93860526-21880b00-fcb7-11ea-8792-7d0ea5bafbfd.png



And plot with the minDP10 filter
```
/Users/alexjvr/2020.postdoc/Velocity/E3/ANGSD_FINAL/SFS/Score.Mark

Rscript script2.R MODE.CADCXM010000001.1.minDP20.MinIND10.beagle.gz MODE.minDP20.scorevec
[1] "number of individuals is  33"
[1] "maximum 'N_eff' is  32.9684065049473"
[1] "minimum 'N_eff' is  7.74582069797937"

##R
library(dplyr)
library(ggplot2)

MODE.GL <- read.table(gzfile("MODE.CADCXM010000001.1.minDP20.MinIND10.beagle.gz"), header=T)

dim(MODE.GL)
[1] 15556   102

#get a list of marker names
MODE.GL.names <- MODE.GL$marker
MODE.GL.names <- as.data.frame(MODE.GL.names)
MODE.GL.names$pos <- MODE.GL$marker
MODE.GL.names$pos <- gsub("CADCXM010000001.1_", "", MODE.GL.names$pos)
MODE.GL.names$pos <- as.numeric(MODE.GL.names$pos)

CDX.fst <- read.table("MODC.MODE.CDX.minDP20minInd10.win1bp.fst", header=T) #read in Fst. Remember to add "fst" as final column name in nano
colnames(CDX.fst) <- c("index", "marker", "pos", "N", "fst")

MODE.score <- read.table("MODE.minDP20.scorevec", header=T)  ##calculated score
MODE.score$pos <- MODE.GL.names$pos  
colnames(MODE.score) <- c("score", "pos")

MODE.GL.sub <- left_join(CDX.fst, MODE.score, by="pos")   ##join right hand side df to left and keep all rows that overlap with left df (using pos)

ggplot(MODE.GL.sub, aes(y=fst, x=score))+ geom_point()
```

Based on 11519 bp

![alt_txt][score2]

[score2]:https://user-images.githubusercontent.com/12142475/93649253-769bf680-fa03-11ea-95bb-3853c3c7fc1d.png



## Filter for missingness and check all depth plots

There is a high level of missingness in the data, and I've identified [individuals that sequenced poorly](https://github.com/alexjvr1/Velocity2020/blob/master/Missingness_Plots.md). I'll re-estimate GL scores after removing these indivs. I am also changing the minDP and minInd filters for all of the ANGSD datasets: 




I'm using contig CADCXM010000001.1 (152182bp) for the intial tests as these run really quickly

#### How many loci do we recover with minDP and minInd filters compared with minimal filters

I didn't realise that the minDP filters worked on global depth: 

-setMinDepth xx # discard site if total sequencing depth (all indivs together) is below xx

-setMaxDepth xx  # discard site if total depth (all indivs added together) is above xx

-setMinDepthInd xx # Discard individual if sequencing depth is below xx. This filter is only applied to analysis which uses allele count (-doCounts). Thus anything that uses GLs (e.g. SAF estimation) does not use this filter. 

To get around this issue I will use the global depth and minInd to filter: 

-minInd 18

-setMinDepth 36 ##i.e. 2X per individual for the min 18 genotyped individuals. 


```
#!/bin/bash
#PBS -N E3_MODC.SAF  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=20:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)
##PBS -t 1-87 #array job

#Set filters
N="36"
MININD="18"
MINMAF=""
MINQ="20"
minMAPQ="20"
minDP="36"
maxDP="621"
POP="MODC"
C="50"
POP="MODC"
POPLIST="MODC.36.poplist"
SPECIESDIR="/newhome/aj18951/E3_Aphantopus_hyperantus_2020"
PP=0 #use all reads. Flag 1 uses only proper pairs, but	MUS has	merged reads. NB to filter for proper pair reads in the bamfiles using samtools before this point

#OUTNAME="MODC"

#run job in working directory
cd $PBS_O_WORKDIR 

#load modules
module load languages/gcc-6.1
angsd=~/bin/angsd/angsd

#Define variables
#REGION=$(sed "${PBS_ARRAYID}q;d" regions)
REGION="LR761675.1"

#estimate SAF for modern expanding population using ANGSD

time $angsd -b $POP.$N.poplist -checkBamHeaders 1 -minQ 20 -minMapQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs $PP -r $REGION \
-GL 1 -doSaf 1 -anc $SPECIESDIR/RefGenome/*fna -ref $SPECIESDIR/RefGenome/*fna -doCounts 1 -setMinDepthInd $minDP -setMaxDepth $maxDP -doMajorMinor 4\
 -out $SPECIESDIR/04b_ANGSD_FINAL/SFS_and_Fst/$POP/$POP.$REGION -C $C -baq 1 -dumpCounts 2 -doDepth 1 -doGlf 2
```


Outputs:
```
CADCXM010000001.1 (len=152 kb)
		Total sites analyzed	minDP2(no minInd)	minDP36.minInd18	minDP20.minInd10	GATK
MODE (N=33/40) 	54k (36%)		43k			10k			15.5k			15.5k 

MODC (N=36/38)	57k (36%)		49k			8.9k			13.5			13.5k

MUS (N=24/48)	11k (7%)					351			691			691			


LR761675 (len=6.2 Mb)
		Total sites analyzed	setminDPInd2	minDP36.minInd18	minDP20.minInd10		GATK
MODE (N=33/40) 	5.7M (93%)					5.0M			5.3M			5.3M (87%)

MODC (N=36/38)	5.7M (93%)					4.7M			5.2M			5.2M (85%)

MUS (N=24/48)	4.7M (77%)		4.7M			0.1M			2.1M			2.1M (34%)

```

#### What is the distribution of missingness in these new datasets? 

Use the methods described [here](https://github.com/alexjvr1/Velocity2020/blob/master/Missingness_Plots.md) to estimate missingness



#### GL models GATK vs Samtools

I recovered exactly the same loci for all three datasets (MODC, MODE, MUS) using these strict filters. 



#### Run the same scripts for LR761675.1

1. Diff between GATK and Samtools? None - I recover the same loci


2. Proportion of loci retained after filters? 34% (MUS) and >80% (MOD)


#### Create SFS and hence fst and diversity estimates for these new files

On server 

MODC vs MODE
```
module load languages/gcc-6.1

#Use the unfolded SAF for each population to create a folded SFS
#realSFS pop1.unfolded.saf.idx pop2.unfolded.saf.idx -fold 1 >folded.sfs

~/bin/angsd/misc/realSFS MODC/MODC.CADCXM010000001.1.minDP20.MinIND10.saf.idx MODE/MODE.CADCXM010000001.1.minDP20.MinIND10.saf.idx -fold 1 > MODC.MODE.CDX.minDP20.minInd10.sfs

#Use the folded SFS and unfolded SAFs to prepare inputs for Fst calculation
#realSFS fst index pop1.unfolded.saf.idx pop2.unfolded.saf.idx -sfs folded.sfs -fstout persite

~/bin/angsd/misc/realSFS fst index MODC/MODC.CADCXM010000001.1.minDP20.MinIND10.saf.idx MODE/MODE.CADCXM010000001.1.minDP20.MinIND10.saf.idx -sfs MODC.MODE.CDX.minDP20.minInd10.sfs -fold 1 -fstout MODC.MODE.CDX.minDP20.minInd10.persite

 ~/bin/angsd/misc/realSFS fst stats  MODC.MODE.CDX.minDP20.minInd10.persite.fst.idx 
	-> Assuming idxname:MODC.MODE.CDX.minDP20.minInd10.persite.fst.idx
	-> Assuming .fst.gz file: MODC.MODE.CDX.minDP20.minInd10.persite.fst.gz
	-> FST.Unweight[nObs:11521]:0.016994 Fst.Weight:0.111027
0.016994	0.111027


##Estimate Fst every 1bp (non-overlapping windows)
#realSFS fst stat2 persite.fst.idx -win XXXX -step XXXX >window.fst

~/bin/angsd/misc/realSFS fst stats2 MODC.MODE.CDX.minDP20.minInd10.persite.fst.idx -win 1 -step 2 > MODC.MODE.minDP20.minInd10.win1bp.fst

win:1 step:1
nSites:11521 (7%)
```

##### same for LR75 dataset


And same for the MODC:MODE LRxx75 dataset
```
~/bin/angsd/misc/realSFS MODC/MODC.LR761675.1.minDP20.MinIND10.saf.idx MODE/MODE.LR761675.1.minDP20.MinIND10.saf.idx -fold 1 > MODC.MODE.LR75.minDP20.minInd10.sfs

~/bin/angsd/misc/realSFS MODC/MODC.LR761675.1.minDP20.MinIND10.saf.idx MODE/MODE.LR761675.1.minDP20.MinIND10.saf.idx -sfs MODC.MODE.LR761675.1.minDP20.MinIND10.fold.sfs -fold 1 -fstout MODC.MODE.LR75.minDP20.minIND10.fstout


~/bin/angsd/misc/realSFS fst stats  MODC.MODE.LR75.minDP20.minIND10.fstout.fst.idx
	-> Assuming idxname:MODC.MODE.LR75.minDP20.minIND10.fstout.fst.idx
	-> Assuming .fst.gz file: MODC.MODE.LR75.minDP20.minIND10.fstout.fst.gz
	-> FST.Unweight[nObs:5092573]:0.015304 Fst.Weight:0.131845
0.015304	0.131845


#using 5,092,713 sites from each pop (of 6196582bp = 82% coverage)

##create fst for each bp 

~/bin/angsd/misc/realSFS fst stats2 MODC.MODE.LR75.minDP20.minIND10.fstout.fst.idx -win 1 -step 2 > MODC.MODE.LR75.minDP20.minInd10.win1bp.fst


###Create fst output in windows for comparison with previous dataset

~/bin/angsd/misc/realSFS fst stats2 MODC.MODE.LR75.minDP20.minIND10.fstout.fst.idx -win 50000 -step 10000 > MODC.MODE.LR75.minDP20.minInd10.win50k.step10k.fst2

#using 5,092,713 sites from each pop (of 6196582bp = 82% coverage)
```



MODC:MUS LRxx75 dataset
```
~/bin/angsd/misc/realSFS MODC/MODC.LR761675.1.minDP20.MinIND10.saf.idx MUS/MUS.LR761675.1.minDP20.MinIND10.saf.idx -fold 1 > MODC.MUS.LR75.minDP20.minInd10.fold.sfs

#index for fst calculation
~/bin/angsd/misc/realSFS fst index MODC/MODC.LR761675.1.minDP20.MinIND10.saf.idx MUS/MUS.LR761675.1.minDP20.MinIND10.saf.idx -sfs MODC.MUS.LR75.minDP20.minInd10.fold.sfs -fold 1 -fstout MODC.MUS.LR75.minDP20.minIND10.fstout

#global fst
~/bin/angsd/misc/realSFS fst stats  MODC.MUS.LR75.minDP20.minIND10.fstout.fst.idx
	-> Assuming idxname:MODC.MUS.LR75.minDP20.minIND10.fstout.fst.idx
	-> Assuming .fst.gz file: MODC.MUS.LR75.minDP20.minIND10.fstout.fst.gz
	-> FST.Unweight[nObs:2084592]:0.018856 Fst.Weight:0.067054
0.018856	0.067054

#window based fst
~/bin/angsd/misc/realSFS fst stats2 MODC.MUS.LR75.minDP20.minIND10.fstout.fst.idx -win 50000 -step 10000 > MODC.MUS.LR75.minDP20.minInd10.win50k.step10k.fst

#using 2,084,592 sites from each pop (of 6196582bp = 33% coverage)

```

MODE:MUS LRxx75
```
#create the folded sfs
~/bin/angsd/misc/realSFS MODE/MODE.LR761675.1.minDP20.MinIND10.saf.idx MUS/MUS.LR761675.1.minDP20.MinIND10.saf.idx -fold 1 > MODE.MUS.LR75.minDP20.minInd10.fold.sfs

#index the sfs for window-based Fst estimation
~/bin/angsd/misc/realSFS fst index MODE/MODE.LR761675.1.minDP20.MinIND10.saf.idx MUS/MUS.LR761675.1.minDP20.MinIND10.saf.idx -sfs MODE.MUS.LR75.minDP20.minInd10.fold.sfs -fold 1 -fstout MODE.MUS.LR75.minDP20.minIND10.fstout

#Global fst estimate
~/bin/angsd/misc/realSFS fst stats  MODE.MUS.LR75.minDP20.minIND10.fstout.fst.idx
	-> Assuming idxname:MODE.MUS.LR75.minDP20.minIND10.fstout.fst.idx
	-> Assuming .fst.gz file: MODE.MUS.LR75.minDP20.minIND10.fstout.fst.gz
	-> FST.Unweight[nObs:2079057]:0.041143 Fst.Weight:0.168676
0.041143	0.168676

#window based fst
~/bin/angsd/misc/realSFS fst stats2 MODE.MUS.LR75.minDP20.minIND10.fstout.fst.idx -win 50000 -step 10000 > MODE.MUS.LR75.minDP20.minInd10.win50k.step10k.fst

win:50000 step:10000
nSites:2079057 (~33% coverage of genome)

```


### Thetas

Estimate diversity measures for each of the cleaned pops and LR75
```
#Estimate the folded sfs for population 
~/bin/angsd/misc/realSFS MODE/MODE.LR761675.1.minDP20.MinIND10.saf.idx -fold 1 > MODE/MODE.LR761675.1.minDP20.MinIND10.sfs
~/bin/angsd/misc/realSFS MODC/MODC.LR761675.1.minDP20.MinIND10.saf.idx -fold 1 > MODC/MODC.LR761675.1.minDP20.MinIND10.sfs
~/bin/angsd/misc/realSFS MUS/MUS.LR761675.1.minDP20.MinIND10.saf.idx -fold 1 > MUS/MUS.LR761675.1.minDP20.MinIND10.sfs


#prep the saf files
~/bin/angsd/misc/realSFS saf2theta MODE/MODE.LR761675.1.minDP20.MinIND10.saf.idx -sfs MODE/MODE.LR761675.1.minDP20.MinIND10.sfs -outname MODE/MODE.LR761675.1.minDP20.MinIND10
~/bin/angsd/misc/realSFS saf2theta MODC/MODC.LR761675.1.minDP20.MinIND10.saf.idx -sfs MODC/MODC.LR761675.1.minDP20.MinIND10.sfs -outname MODC/MODC.LR761675.1.minDP20.MinIND10
~/bin/angsd/misc/realSFS saf2theta MUS/MUS.LR761675.1.minDP20.MinIND10.saf.idx -sfs MUS/MUS.LR761675.1.minDP20.MinIND10.sfs -outname MUS/MUS.LR761675.1.minDP20.MinIND10


#Calculate diversity stats

~/bin/angsd/misc/thetaStat do_stat MODE/MODE.LR761675.1.minDP20.MinIND10.thetas.idx -win 50000 -step 10000  -outnames MODE/MODE.LR761675.1.minDP20.MinIND10.thetasWindow.gz
~/bin/angsd/misc/thetaStat do_stat MODC/MODC.LR761675.1.minDP20.MinIND10.thetas.idx -win 50000 -step 10000  -outnames MODC/MODC.LR761675.1.minDP20.MinIND10.thetasWindow.gz
~/bin/angsd/misc/thetaStat do_stat MUS/MUS.LR761675.1.minDP20.MinIND10.thetas.idx -win 50000 -step 10000  -outnames MUS/MUS.LR761675.1.minDP20.MinIND10.thetasWindow.gz

#And to plot against depth: 
~/bin/angsd/misc/thetaStat do_stat MODE/MODE.LR761675.1.minDP20.MinIND10.thetas.idx -win 1 -step 1  -outnames MODE/MODE.LR761675.1.minDP20.MinIND10.win1.step1.thetasWindow.gz
~/bin/angsd/misc/thetaStat do_stat MODC/MODC.LR761675.1.minDP20.MinIND10.thetas.idx -win 1 -step 1  -outnames MODC/MODC.LR761675.1.minDP20.MinIND10.win1.step1.thetasWindow.gz
~/bin/angsd/misc/thetaStat do_stat MUS/MUS.LR761675.1.minDP20.MinIND10.thetas.idx -win 1 -step 1  -outnames MUS/MUS.LR761675.1.minDP20.MinIND10.win1.step1.thetasWindow.gz


```



Plot against depth 
```
scp bluecp3:/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04b_ANGSD_FINAL/SFS_and_Fst/M*/M*.LR761675.1.minDP20.MinIND10.pos.gz  .
scp bluecp3:/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04b_ANGSD_FINAL/SFS_and_Fst/M*/M*win1.step1.thetasWindow.gz.pestPG  .
##Delete after using as they're really big files


##R
library(reshape2)
library(dplyr)
library(ggplot2)

MUS.count <- read.table(gzfile("MUS.LR761675.1.minDP20.MinIND10.pos.gz"), header=T)
MODE.count <- read.table(gzfile("MODE.LR761675.1.minDP20.MinIND10.pos.gz"), header=T)
MODC.count <- read.table(gzfile("MODC.LR761675.1.minDP20.MinIND10.pos.gz"), header=T)

#remember to remove the #from the header in nano to read the headers into  R
#If the file is very big it'll be better to use sed:
#sed -i 's/^#//' *pestPG
MUS.div <- read.table("MUS*pestPG", header=T)
MODC.div <- read.table("MODC.LR761675.1.minDP20.MinIND10.win1.step1.thetasWindow.gz.pestPG", header=T)
MODE.div <- read.table("MODE.LR761675.1.minDP20.MinIND10.win1.step1.thetasWindow.gz.pestPG", header=T)

MUS <- left_join(MUS.div, MUS.count, by="pos")
MODC <- left_join(MODC.div, MODC.count, by="pos")
MODE <- left_join(MODE.div, MODE.count, by="pos")

pdf("E3.tWvsDepth.pdf")
ggplot(MUS[1:100000,], aes(x=totDepth, y=tW))+geom_point()+ggtitle("tW vs Depth, 1M loci: MUS")
ggplot(MODC[1:100000,], aes(x=totDepth, y=tW))+geom_point()+ggtitle("tW vs Depth, 1M loci: MODC")
ggplot(MODE[1:100000,], aes(x=totDepth, y=tW))+geom_point()+ggtitle("tW vs Depth, 1M loci: MODE")
dev.off()
```


![alt_txt][tWvsDepth.MUS]

[tWvsDepth.MUS]:https://user-images.githubusercontent.com/12142475/94028045-6a24ee80-fdb3-11ea-817a-6f0f822bf4c5.png

Low depth does not seem to be a problem. But I need to adjust my maxDP filter. 



## MaxDepth

I need a filter for maxDP as well.

The maxDP cut-off I usually use is meanDP + (2xSD of depth). 

Heng Li has found a correlation between depth and site quality [preprint here](http://arxiv.org/pdf/1404.0929.pdf), and a full blog post on Blue Collar Bioinformatics [here](http://bcb.io/2014/05/12/wgs-trio-variant-evaluation/). 

"In short, with whole genome samples, it was found that high coverage can lead to inflated locus quality scores. Heng proposed that for read depths greater than the mean depth plus 2-3 times the square root of mean depth that the quality score will be twice as large as the depth in real variants and below that value for false variants."
quoted from the [ddocent site by Jonathan Puritz](http://bcb.io/2014/05/12/wgs-trio-variant-evaluation/)

To check that my filters are appropriate I'll plot depth vs sequence quality: 

```


```



## Depth vs Number of usable variants

Problem: Finding comparative diversity estimates is proving really difficult because the Museum dataset is not comparable to the Modern datasets. 

I've tried various things to match the data up: 1) filtering to the same mean depth, 2) subsampling to the same number of reads, 3) using only variants we're certain of (e.g. with genotype calls). 

The museum data seems to be of worse quality than the modern data. e.g. if we compare the mean depth of the filtered bam files to the number of loci for which we have data usable by ANGSD

```
# use the *depthSample outputs from ANGSD
# This gives the number of sites usable by ANGSD per individual 
MODC <- read.table("MODC.LR761675.1.NoPCRDup.OCT19.depthSample")
MODE <- read.table("MODE.LR761675.1.noPCRDUP.OCT19.depthSample")
MUS <- read.table("MUS.LR761675.1.NoPCRDup.OCT19.depthSample")
MODC.t <- t(MODC)
MODE.t <- t(MODE)
MUS.t <- t(MUS)
MODC.t <- MODC.t[2:nrow(MODC.t),]
MODE.t <- MODE.t[2:nrow(MODE.t),]
MUS.t <- MUS.t[2:nrow(MUS.t),]
colSums(MODC.t)
colSums(MUS.t)
colSums(MODE.t)


#Mean depth is estimted with samtools from the filtered bam files. 
#But we can also estimate depth from the depthSample frequency files

##MODE
MODE <- read.table("MODE.LR761675.1.noPCRDUP.OCT19.depthSample")
MODE.t <- t(MODE)
#colnames(MODE.t) <- paste("Indiv", 1:ncol(MODE.t), sep="")
COV <- 0:(nrow(MODE.t)-1) #create the vector with coverage for each row (1x, 2x, etc)
MODE.t.TIMESCOV <- MODE.t*COV #multiply read count by coverage class
colSums(MODE.t.TIMESCOV)/colSums(MODE.t) #mean of a frequency table = sum(frequencies x class)/sum(frequencies). This gives the mean coverage

##The above includes class 0x (i.e. loci for which we don't have data for that individual). We can exclude these: 
#This gives the mean coverage per indiv for all loci with data
colSums(MODE.t.TIMESCOV[2:nrow(MODE.t.TIMESCOV),],)/colSums(MODE.t[2:nrow(MODE.t),],) #mean of a frequency table = sum(frequencies x class)/sum(frequencies). 

##MODC
MODC <- read.table("MODC.LR761675.1.NoPCRDup.OCT19.depthSample")
MODC.t <- t(MODC)
#colnames(MODC.t) <- paste("Indiv", 1:ncol(MODC.t), sep="")
COV <- 0:(nrow(MODC.t)-1) #create the vector with coverage for each row (1x, 2x, etc)
MODC.t.TIMESCOV <- MODC.t*COV #multiply read count by coverage class
colSums(MODC.t.TIMESCOV)/colSums(MODC.t) #mean of a frequency table = sum(frequencies x class)/sum(frequencies). This gives the mean coverage

##The above includes class 0x (i.e. loci for which we don't have data for that individual). We can exclude these: 
#This gives the mean coverage per indiv for all loci with data
colSums(MODC.t.TIMESCOV[2:nrow(MODC.t.TIMESCOV),],)/colSums(MODC.t[2:nrow(MODC.t),],) #mean of a frequency table = sum(frequencies x class)/sum(frequencies). 

##MUS
MUS <- read.table("MUS.LR761675.1.NoPCRDup.OCT19.depthSample")
MUS.t <- t(MUS)
#colnames(MUS.t) <- paste("Indiv", 1:ncol(MUS.t), sep="")
COV <- 0:(nrow(MUS.t)-1) #create the vector with coverage for each row (1x, 2x, etc)
MUS.t.TIMESCOV <- MUS.t*COV #multiply read count by coverage class
colSums(MUS.t.TIMESCOV)/colSums(MUS.t) #mean of a frequency table = sum(frequencies x class)/sum(frequencies). This gives the mean coverage

##The above includes class 0x (i.e. loci for which we don't have data for that individual). We can exclude these: 
#This gives the mean coverage per indiv for all loci with data
colSums(MUS.t.TIMESCOV[2:nrow(MUS.t.TIMESCOV),],)/colSums(MUS.t[2:nrow(MUS.t),],) #mean of a frequency table = sum(frequencies x class)/sum(frequencies). 

```

# Obs het vs depth


We can also estimate the per sample observed heterozygosity by creating saf/sfs for each individual: 
```
#!/bin/bash
#PBS -N E3_LR.MUS.SAF  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=20:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)
#PBS -t 1-24 #array job

#Set filters
N="24"
MININD="10"
MINQ="20"
minMAPQ="20"
minDP="3"
maxDP="621"
POP="MUS"
C="50"
POPLIST="MUS.24.poplist.dedup"
SPECIESDIR="/newhome/aj18951/E3_Aphantopus_hyperantus_2020"
PP=0 #use all reads. Flag 1 uses only proper pairs, but	MUS has	merged reads. NB to filter for proper pair reads in the bamfiles using samtools before this point

#run job in working directory
cd $PBS_O_WORKDIR 

#load modules
module load languages/gcc-6.1
angsd=~/bin/angsd/angsd

#Define variables
#REGION=$(sed "${PBS_ARRAYID}q;d" regions)
REGION="LR761675.1"
INDIV=$(sed "${PBS_ARRAYID}q;d" MUS.24.poplist.dedup)
NAME=`echo ${INDIV} | awk -F "bams/" '{print $2}' | awk -F ".merged" '{print $1}'`

#estimate SAF for modern expanding population using ANGSD

time $angsd -i ${INDIV} -checkBamHeaders 1 -minQ 20 -minMapQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs $PP -r $REGION \
-GL 1 -doSaf 1 -anc $SPECIESDIR/RefGenome/*fna -ref $SPECIESDIR/RefGenome/*fna -doCounts 1 -setMinDepthInd $minDP -setMaxDepth $maxDP -doMajorMinor 4\
 -out $SPECIESDIR/04b_ANGSD_FINAL/SFS_and_Fst/OCT19.DeDup.minDP3/$POP.$NAME.$REGION.NoPCRDup.OCT19 -C $C -baq 1 -dumpCounts 2 -doDepth 1 -doGlf 2 -noTrans 1
```

Estimate the per sample SFS to obtain the individual observed het. The second value is the number of observed heterozygotes
```
for i in $(ls *AH*saf.idx); do ~/bin/angsd/misc/realSFS $i > $i.ml; done
```

And find the proportion for each sample
```
for i in $(ls *ml); do awk '{sum+=$1+$2+$3} END {print $2/sum}' $i; done

##In R
library(ggplot2)
depth <- read.table(LociVsDepth, header=T)
head(depth)
          Indiv    Loci Depth.sam pop NoDups.Loci  ObsHet Depth.ANGSD
1 AH-01-1900-01 1775320       5.6 MUS       59398 0.00077         8.3
2 AH-01-1900-04  932681       2.3 MUS       27656 0.00073         2.2
3 AH-01-1900-06 1021676       2.7 MUS       31678 0.00073         2.6
4 AH-01-1900-08 1018349       2.7 MUS       32190 0.00074         2.7
5 AH-01-1900-09 1398015       4.0 MUS       46278 0.00081         4.2
6 AH-01-1900-10  952989       2.2 MUS       26540 0.00065         2.4
  Depth.ANGSD_minusmissing
1                      9.8
2                      5.6
3                      5.7
4                      6.0
5                      6.4
6                      6.3

ggplot(depth, aes(x=Depth.ANGSD_minusmissing, y=ObsHet, colour=pop))+geom_point()+ylab("Proportion of observed Heterozygote sites from ANGSD (>3X Depth)")+xlab("Mean Depth of filtered ANGSD file")

```

![alt_txt][depthVShet]

[depthVShet]:https://user-images.githubusercontent.com/12142475/96591161-dd505080-12de-11eb-892c-fb778a181e99.png


Depth vs number of loci recovered
```
ggplot(depth, aes(x=Depth.ANGSD_minusmissing, y=NoDups.Loci, colour=pop))+geom_point()+ylab("Number of loci with data in ANGSD (>3X Depth)")+xlab("Mean Depth of filtered ANGSD file")

ggplot(depth, aes(x=Depth.sam, y=NoDups.Loci, colour=pop))+geom_point()+ylab("Number of loci with data in ANGSD (>3X Depth)")+xlab("Mean Depth of filtered BAM file")

```

![alt_txt][ANGSDdepthvsANGSDloci]

[ANGSDdepthvsANGSDloci]:https://user-images.githubusercontent.com/12142475/96593583-91eb7180-12e1-11eb-8b97-e0cd1a62ad3d.png


![alt_txt][BAMdepthvsANGSDloci]

[BAMdepthvsANGSDloci]:https://user-images.githubusercontent.com/12142475/96594037-16d68b00-12e2-11eb-9b3e-fa717638e033.png
