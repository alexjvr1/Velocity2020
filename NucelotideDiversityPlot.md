# Plot Nucleotide diversity across the genome

See Fig 3 in [Feng et al. 2019](https://reader.elsevier.com/reader/sd/pii/S0960982218316099?token=65EE80B75A634FD527312349DBECB9E11DFF80648448692DC1E7DA07DD7E7A6552DFB5DD8B92B7420E4E045F526A4075#go_to_%22%C3%BE%C3%BF\u0000f\u0000l\u0000i\u0000n\u0000k\u00004%22)


Thetas calculated in ANGSD in windows (-win 50kb -step 10kb)

## Plots in R

Copy diversity estimates to mac
```
pwd
/Users/alexjvr/2020.postdoc/Velocity/E3/ANGSD_FINAL/DiversityEstimates

scp bluecp3:/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04b_ANGSD_FINAL/SFS_and_Fst/M*/*pestPG .

#remove # in header of each file to access the colnames
```


R
```
library(reshape2)
library(ggplot2)

MODC.div <- read.table("MODC.LR761675.1.minDP20.MinIND10.thetasWindow.gz.pestPG", header=T)
MODE.div <- read.table("MODE.LR761675.1.minDP20.MinIND10.thetasWindow.gz.pestPG", header=T)
MUS.div <- read.table("MUS.LR761675.1.minDP20.MinIND10.thetasWindow.gz.pestPG", header=T)

head(MUS.div)
  X.indexStart.indexStop..firstPos_withData.lastPos_withData..WinStart.WinStop.
1                                        (2006,15523)(10505,60004)(10000,60000)
2                                        (5037,20226)(20000,70000)(20000,70000)
3                                        (6658,24122)(30000,80088)(30000,80000)
4                                       (10554,28003)(40305,90391)(40000,90000)
5                                     (12427,30833)(50000,100704)(50000,100000)
6                                     (15523,33806)(60004,111849)(60000,110000)
         Chr WinCenter        tW       tP       tF       tH       tL    Tajima
1 LR761675.1     35000  93.87998 49.42940 286.0133 18.69926 34.06433 -1.731692
2 LR761675.1     45000 112.42744 64.86303 325.6139 25.48888 45.17596 -1.549229
3 LR761675.1     55000 127.18538 73.59931 368.2641 29.05320 51.32625 -1.543967
4 LR761675.1     65000 125.75956 72.24172 365.4616 28.68044 50.46108 -1.559386
5 LR761675.1     75000 132.53109 75.67231 388.6090 30.33121 53.00176 -1.572537
6 LR761675.1     85000 133.71968 76.20060 390.0751 30.06521 53.13291 -1.576734
        fuf       fud     fayh      zeng nSites
1 -4.695567 -5.347124 0.249719 -0.481952  13517
2 -4.335533 -4.976488 0.267276 -0.452760  15189
3 -4.339249 -4.987789 0.267355 -0.451619  17464
4 -4.366201 -5.014383 0.264403 -0.453349  17449
5 -4.425061 -5.088568 0.261167 -0.454424  18406
6 -4.399422 -5.049676 0.263384 -0.456384  18283


##Headers explained
# tW: Watterson's theta
# tP:
# tF:
# tH:
# tL:
# Tajima: Tajima's D
# fuf: 
# fud:
# fayh: Fay's H
# Zeng: 
# nsites: number of sites in the window 

#These outputs are the sum of estimates across the window. So to get a mean estimate for the window we need to divide by the number of sites

MODC.div$tW.rat <- MODC.div$tW/MODC.div$nSites 
MODE.div$tW.rat <- MODE.div$tW/MODE.div$nSites
MUS.div$tW.rat <- MUS.div$tW/MUS.div$nSites

##combine into one dataset and melt to plot with ggplot
library(dplyr)

Pop3.NucDiv <- left_join(MODE.div, MODC.div, by="WinCenter", suffix=c(".e", ".c"))
Pop3.NucDiv <- left_join(Pop3.NucDiv, MUS.div, by="WinCenter", suffix=c(".mod", ".mus"))
NucDiv.melt <- melt(Pop3.NucDiv, id.vars=c("WinCenter", "tW.rat.c", "tW.rat.e", "tW.rat"))
colnames(NucDiv.melt) <- c("WinCenter", "MODC", "MODE", "MUS", "variable", "value")
NucDiv.melt2 <- (melt(NucDiv.melt[,2:4]))

pdf("E3.tW.minDP20.minInd10.pdf")
ggplot(NucDiv.melt2, aes(x=value, fill=variable)) + geom_histogram(alpha=0.6, bins=30) + scale_fill_manual(values=c("#2E8B57", "#46CC7C", "#DAA520"))+ggtitle("Watterson's theta for Ringlet: LR75")+xlab("Watterson's theta")+ylab("frequency")
dev.off()

```

![alt_txt][tW]

[tW]:https://user-images.githubusercontent.com/12142475/94013602-368d9880-fda2-11ea-8b04-70e9d4c29a9a.png


```
##And plot tW across the chromosome for each population

NucDiv.melt2 <- melt(NucDiv.melt[,1:4], id.vars="WinCenter") 
head(NucDiv.melt2)
  WinCenter variable       value
1     35000     MODC 0.004617283
2     45000     MODC 0.005054221
3     55000     MODC 0.005237145
4     65000     MODC 0.005299250
5     75000     MODC 0.005200878
6     85000     MODC 0.005375766

pdf("E3.tW_linegraph.pdf")
ggplot(NucDiv.melt2, aes(x=WinCenter, y=value, color=variable))+geom_line()+scale_color_manual(values=c("#2E8B57", "#46CC7C", "#DAA520"))+ggtitle("Watterson's theta for Ringlet: LR75")+xlab("LR75 position (bp)")+ylab("Watterson's theta")
dev.off()
```

![alt_txt][tW2]

[tW2]:https://user-images.githubusercontent.com/12142475/94015633-0b587880-fda5-11ea-8581-b4fa275dffbd.png

### Use Adobe Illustrator to turn these into a pretty plot

![alt_txt][E3.tWFig]

[E3.tWFig]:https://user-images.githubusercontent.com/12142475/94018753-d77f5200-fda8-11ea-85e6-ce7b1c51f24d.png


And similarly for nucleotide diversity

![alt_txt][nucDiv]

[nucDiv]:https://user-images.githubusercontent.com/12142475/95187939-adc02500-07c3-11eb-9001-5256a6b121d1.png

On Dropbox:
```
~/Dropbox/NERC\ leps\ 2015/Velocity\ WP1\ analysis/1_RingletMS_2020/E3_Fig2_DiversityStats.pdf
```


A comparison of NucDiv from ANGSD and from filtered genotype calls from samtools/bcftools call pipeline

![alt_txt][nucDiv.comp]

[nucDiv.comp]:https://user-images.githubusercontent.com/12142475/95516128-aca31980-09b6-11eb-9805-cac045abbde8.png


Tajima's D

![alt_txt][TajD]

[TajD]:https://user-images.githubusercontent.com/12142475/95207692-0998a700-07e0-11eb-820c-3961f79ecac3.png




## Per site vs window-based estimates

I've calculated Watterson's theta from ANGSD as the cumulative window based estimate (tW) divided by the number of sites with data in that window (nsites). 

To see if this tallys with ANGSDs per site estimates, I will manually calculate Watterson's theta from the per site estimates. 

Per site estimates can be obtained from [ANGSD](http://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests)
```
pwd
/Users/alexjvr/2020.postdoc/Velocity/E3/ANGSD_FINAL/DiversityEstimates/OCT26

##R
library(dplyr)
library(reshape2)
library(ggplot2)

MODC.HR <- read.table("../MODC.head1Mil.HumanReadable", header=T)
MODE.HR <- read.table("../MODE.head1Mil.HumanReadable", header=T)
MUS.HR <- read.table("../MUS.head1Mil.HumanReadable", header=T)

pop12 <- left_join(MODE.HR, MODC.HR, by="Pos")
pop123 <- left_join(pop12, MUS.HR, by="Pos")
pop123 <- pop123[complete.cases(pop123),]
pop123.m <- melt(pop123, id=c("Pos", "Watterson.x", "Watterson.y", "Watterson"))
colnames(pop123.m) <- c("Pos", "MODC", "MODE", "MUS", "variable", "value")

#log values are reported. Here I change them to real values
pop123.m$MODC.exp <- exp(pop123.m$MODC)
pop123.m$MODE.exp <- exp(pop123.m$MODE)
pop123.m$MUS.exp <- exp(pop123.m$MUS)

#To calculate window-based estimates: 
##This isn't exactly right because the windows should be from pos x+50k, not the first 50k datapoints. BUT it gives us an idea of the trend in the data. 
require(zoo)
test <- rollapply(pop123.m$MODE.exp, width=50000, by=10000, FUN=mean)
test.t <- as.data.frame(test)
colnames(test.t)[1] <- "MODE"
test.t$MUS <- rollapply(pop123.m$MUS.exp, width=50000, by=10000, FUN=mean)
test.t$MODC <- rollapply(pop123.m$MODC.exp, width=50000, by=10000, FUN=mean)
test.t.melt <- melt(test.t, id="POS")
head(test.t.melt)
ggplot(test.t.melt, aes(x=POS, y=value, by=variable))+geom_line()


#Here we want to look at the distribution of values comparing all three populations
pop123.mm <- melt(pop123.m[1:4],id="Pos")
pop123.mm$exp <- exp(pop123.mm$value) #log values are reported. Here I change them to real values
#With invariant sites
ggplot(pop123.mm, aes(x=exp, fill=variable))+geom_histogram()+ scale_fill_manual(values=c("#2E8B57", "#46CC7C", "#DAA520"))

#And without the invariant sites
ggplot(pop123.mm[which(pop123.mm$exp>0.01),], aes(x=exp, fill=variable))+geom_histogram()+ scale_fill_manual(values=c("#2E8B57", "#46CC7C", "#DAA520"))

#Log scores
ggplot(pop123.mm, aes(x=value, fill=variable))+geom_histogram()+ scale_fill_manual(values=c("#2E8B57", "#46CC7C", "#DAA520"))

```

Window-based estimates from single site tW

![alt_txt][SS_window_tW]

[SS_window_tW]:https://user-images.githubusercontent.com/12142475/97296461-9c10f100-1848-11eb-8981-bda851e9c4e6.png


tW per site for all sites
![alt_txt][SS_all]

[SS_all]:https://user-images.githubusercontent.com/12142475/97296023-01181700-1848-11eb-9a64-fcb22ea6cba4.png


tW per site excluding invariant sites
![alt_txt][SS_variantsonly]

[SS_variantsonly]:https://user-images.githubusercontent.com/12142475/97296073-13925080-1848-11eb-846b-f12eb9159c18.png


histogram of log values
![alt_txt][Hist.log]

[Hist.log]:https://user-images.githubusercontent.com/12142475/97295879-cc0bc480-1847-11eb-8719-c119c60f7013.png



## Same sample sizes

It seems that ANGSD doesn't correct for n (sample size) at each site. So I will run a crude correction based on the sample size, although I know that the genotyping rate is different between the different populations. 
```
library(dplyr)

MUS.div$tW.rat <- MUS.div$tW/MUS.div$nSites
MODC.div$tW.rat <- MODC.div$tW/MODC.div$nSites
MODE.div$tW.rat <- MODE.div$tW/MODE.div$nSites

#Divided by sample size
MUS.div$tW.rat.corr <- MUS.div$tW.rat/24
MODC.div$tW.rat.corr <- MODC.div$tW.rat/36
MODE.div$tW.rat.corr <- MODE.div$tW.rat/33

#Corrected for genotyping rate (GR). These numbers are based on what I expect the average GR is - 10/24 indivs for MUS, and around 80-90% for the MOD data. 
MUS.div$tW.rat.corr <- MUS.div$tW.rat*(0.4*24)
MODC.div$tW.rat.corr <- MODC.div$tW.rat*(0.8*36)
MODE.div$tW.rat.corr <- MODE.div$tW.rat*(0.8*33)


Pop3.NucDiv <- left_join(MODE.div, MODC.div, by="WinCenter", suffix=c(".e", ".c"))
Pop3.NucDiv <- left_join(Pop3.NucDiv, MUS.div, by="WinCenter", suffix=c(".mod", ".mus"))
NucDiv.melt3 <- melt(Pop3.NucDiv, id.vars=c("WinCenter", "tW.rat.corr.c", "tW.rat.corr.e", "tW.rat.corr"))
colnames(NucDiv.melt3) <- c("WinCenter", "MODC", "MODE", "MUS", "variable", "value")
NucDiv.melt4 <- (melt(NucDiv.melt3[,2:4]))

pdf("E3.tW.corr.pdf")
ggplot(NucDiv.melt4, aes(x=value, fill=variable)) + geom_histogram(alpha=0.6, bins=30) + scale_fill_manual(values=c("#2E8B57", "#46CC7C", "#DAA520"))+ggtitle("Watterson's theta for Ringlet: LR75")+xlab("Watterson's theta")+ylab("frequency")
dev.off()
```

![alt_txt][tW.sizeCorr]

[tW.sizeCorr]:


I'll reduce the modern datasets to the same sample size as MUS, and reduce missingness in the dataset by accepting min 75% genotyping rate. (i.e. 18/24 individuals)



## Pixy

[Pixy](https://pixy.readthedocs.io/en/latest/installation.html) has been developed specifically to work with patchy data.

See the preprint paper [here](https://www.biorxiv.org/content/10.1101/2020.06.27.175091v1.full.pdf). 

```
#Installed pixy on the server
#bzzjrb/software/pixy

```

We need vcf files with all sites. See below


Load modules 
```
#Callum installed pixy in my home directory from source. This should now run once I've loaded python3
module load languages/python-anaconda3-5.2.0

pixy --version
```


## VCF files AllSites

Change variant calling script
```
#Change -v to 0 in variant calling submission script to output all sites. 

#The raw bcf files need to be processed to "see" missing data
for i in $(ls *bcf); do bcftools filter -S . -O u -e 'FMT/DP=0' $i |bcftools view -O b -o $i.withmissing.bcf; done
```

Use vcftools to find the number of variants within each individual by finding the proportion of missingness for each individual. 
```
module load apps/vcftools-0.1.12b

vcftools --bcf file.bcf --missing-indv

```

Estimate depth from the bam files (I calculated this before using samtools depth)

Use these data to plot number of variants vs depth. 
```
pwd
/Users/alexjvr/2020.postdoc/Velocity/E3/ANGSD_FINAL/DiversityEstimates/OCT26

scp bzzjrbxx/newhome/bzzjrb/E3/03_variants/*imiss .

head vcf.totalsites
Indiv	Loci	Depth.sam	pop	NoDups.Loci	ObsHet	Depth.ANGSD	Depth.ANGSD_minusmissing	Variants	TotalSites
AH-01-1900-01	1775320	5.6	MUS	59398	0.00077	8.3	9.8	148326	3599207
AH-01-1900-04	932681	2.3	MUS	27656	0.00073	2.2	5.6	94013	2232132
AH-01-1900-06	1021676	2.7	MUS	31678	0.00073	2.6	5.7	104071	2484876
AH-01-1900-08	1018349	2.7	MUS	32190	0.00074	2.7	6.0	101351	2419390
AH-01-1900-09	1398015	4.0	MUS	46278	0.00081	4.2	6.4	132229	3190913
AH-01-1900-10	952989	2.2	MUS	26540	0.00065	2.4	6.3	94058	2234742
AH-01-1900-11	1181275	3.2	MUS	39958	0.00079	3.4	6.0	115635	2756403
AH-01-1900-13	1030830	2.6	MUS	28073	0.00071	2.4	6.0	103704	2477715
AH-01-1900-14	967191	2.0	MUS	30263	0.00064	2.4	5.6	89469	2105816

R
library(ggplot2)

data <- read.table("vcf.totalsites", header=T)
pdf("E3.rawvcf_plots.26Oct2020.pdf")
ggplot(data, aes(x=Depth.sam, y=TotalSites, group=pop))+geom_point(aes(colour=pop))+ggtitle("Total sites in raw vcf file")
ggplot(data, aes(x=Depth.sam, y=Variants, group=pop))+geom_point(aes(colour=pop))+ggtitle("Variants in raw vcf file")
ggplot(data, aes(x=Depth.sam, y=propVars, group=pop))+geom_point(aes(colour=pop))+ggtitle("Proportion of sites that are SNPs")
dev.off()
```

![alt_txt][fig1]

[fig1]:https://user-images.githubusercontent.com/12142475/97458835-02723e00-1933-11eb-81d5-cd498dfcf077.png


![alt_txt][fig2]

[fig2]:https://user-images.githubusercontent.com/12142475/97458934-1d44b280-1933-11eb-9955-d108d5dd34bc.png


![alt_txt][fig3]

[fig3]:https://user-images.githubusercontent.com/12142475/97459070-42d1bc00-1933-11eb-873c-60e472aeaa8a.png



### Verify with more chromosomes

To verify if this pattern holds true across chromosomes, I will estimate prop vars for 3 more chromosomes: 

```
LR761672-74

#These numbers are for the full dataset; i.e. with all indivs included

Number of sites

        LR72      LR73        LR74       LR75
MODE   7592213    6794090     6228061    5698793

MODC   7649355    6846041     6255121    5720701

MUS    7013370    5935074     5728727    5071269 

```


For each dataset we'll estimate the mean depth per individual with vcftools from the raw dataset which includes all sites: 
```
#e.g
vcftools --bcf job0001.sites.raw.bcf.withmissing.bcf --depth

#use the idepth file to generate depth per indiv for each chromosome
```

And the missingness for the sites and variants bcf files
```
for i in $(ls *withmissing.bcf); do vcftools --bcf $i --missing-indv; done

```

Copy all the data into a table on my computer: 
```
/Users/alexjvr/2020.postdoc/Velocity/E3/ANGSD_FINAL/DiversityEstimates/OCT26/MultiChrs_vcfstats

R

library(ggplot2)
data <- read.table("MultiChrs_vcfstats", header=T)

pdf("E3.4Chrs_vcfout.plots.29Oct2020.pdf")
ggplot(data, aes(x=MeanDP, y=Sites, group=POP))+geom_point(aes(colour=POP, pch=CHR))+ggtitle("Total sites in raw vcf file")
ggplot(data, aes(x=MeanDP, y=Sites, color=CHR))+geom_point(aes(pch=POP))+ ggtitle("Total sites in raw vcf file")

ggplot(data, aes(x=MeanDP, y=Variants, group=POP))+geom_point(aes(colour=POP, pch=CHR))+ggtitle("Variants in raw vcf file")
ggplot(data, aes(x=MeanDP, y=Variants, group=CHR))+geom_point(aes(colour=CHR, pch=POP))+ggtitle("Variants in raw vcf file")

ggplot(data, aes(x=MeanDP, y=PropVars, group=POP))+geom_point(aes(colour=POP, pch=CHR))+ggtitle("Proportion of sites that are SNPs")
ggplot(data, aes(x=MeanDP, y=PropVars, group=CHR))+geom_point(aes(colour=CHR, pch=POP))+ggtitle("Proportion of sites that are SNPs")

dev.off()


```


![alt_txt][MultiChr1]

[MultiChr1]:


# Pixy estimates

Installed on server: 
```
module load languages/python-anaconda3-5.2.0

/newhome/bzzjrb/

#e.g. calculate for full MUS dataset
/cm/shared/languages/Anaconda3-5.2.0/bin/pixy --stats pi --vcf /newhome/bzzjrb/E3/03_variants/LR75/MUS.LR75.raw.withmissing.vcf --zarr_path /newhome/bzzjrb/E3/pixy/ --reuse_zarr yes --variant_filter_expression 'DP>=10' --invariant_filter_expression 'DP>=10' --outfile_prefix ./MUS.48 --populations /newhome/bzzjrb/MUS.48.names --window_size 10000

/cm/shared/languages/Anaconda3-5.2.0/bin/pixy --stats pi --vcf /newhome/bzzjrb/E3/03_variants/LR75/MODE.LR75.raw.withmissing.33.recode.vcf --zarr_path /newhome/bzzjrb/E3/pixy/ --reuse_zarr yes --variant_filter_expression 'DP>=10' --invariant_filter_expression 'DP>=10' --outfile_prefix ./MODE.test --populations /newhome/bzzjrb/MODE.names --window_size 10000

/cm/shared/languages/Anaconda3-5.2.0/bin/pixy --stats pi --vcf /newhome/bzzjrb/E3/03_variants/LR75/MODC.LR75.raw.withmissing.36.recode.vcf --zarr_path /newhome/bzzjrb/E3/pixy/ --reuse_zarr yes --variant_filter_expression 'DP>=10' --invariant_filter_expression 'DP>=10' --outfile_prefix ./MODC.test --populations /newhome/bzzjrb/MODC.names --window_size 10000
```

copy to mac and plot
```
#
pwd


#R
library(ggplot2)
library(dplyr)
library(reshape2)

MODE <- read.table("MODE.test_pi.txt", head=T)
MODC <- read.table("MODC.test_pi.txt", head=T)

#Plot pixy pi for modern populations
MOD <- bind_rows(MODC, MODE)
ggplot(MOD, aes(x=window_pos_1, y=avg_pi, color=pop))+geom_line()+scale_color_manual(values=c("#2E8B57", "#46CC7C"))+ggtitle("Pixy pi: Modern, unfiltered vcf (meanDP10x)")

#Plot difference in pi for modern populations
MOD2 <- left_join(MODC, MODE, by="window_pos_1", suffix=c(".c", ".e"))
MOD2$diff <- MOD2$avg_pi.c-MOD2$avg_pi.e
ggplot(MOD2, aes(x=window_pos_1, y=diff))+geom_line()+ggtitle("Pixy pi: Difference in modern pi (meanDP10x)")


MUS <- read.table("MUS.48_pi.txt", header=T)
POP3 <- bind_rows(MOD, MUS)

#There's a weird MUS point with pi of 1 and some with 0 that probably aren't true. Second plot removes these
ggplot(POP3, aes(x=window_pos_1, y=avg_pi, color=pop))+geom_line()+scale_color_manual(values=c("#2E8B57", "#46CC7C", "#DAA520"))+ggtitle("Pixy pi: meanDP10x only")
ggplot(POP3[which(POP3$avg_pi<1&POP3$avg_pi>0),], aes(x=window_pos_1, y=avg_pi, color=pop))+geom_line()+scale_color_manual(values=c("#2E8B57", "#46CC7C", "#DAA520"))+ggtitle("Pixy pi: meanDP10x only")

```

![alt_txt][pixy1]

[pixy1]:https://user-images.githubusercontent.com/12142475/99673951-49081380-2a6d-11eb-9b5d-39058610dba1.png


Filter vcf files and recalculate
```
#If I filter for 50% missingness I lose most of the loci in the MUS dataset - we end up with 10% of the marker density compared to modern populations. 
#Using the datasets with poorly sequenced indivs removed, filter for quality, depth and missingness. 
#Pixy calculations

vcftools --vcf MODE.LR75.raw.withmissing.33.recode.vcf --minQ 25 --min-meanDP 6 --max-missing 0.3 --recode --recode-INFO-all --out FILTERED/MODE.LR75.33.filtered.withmissing

After filtering, kept 33 out of 33 Individuals
Outputting VCF file...
After filtering, kept 985263 out of a possible 5698793 Sites
Run Time = 180.00 seconds

vcftools --vcf MODC.LR75.raw.withmissing.36.recode.vcf --minQ 25 --min-meanDP 6 --max-missing 0.3 --recode --recode-INFO-all --out FILTERED/MODC.LR75.36.filtered.withmissing

After filtering, kept 36 out of 36 Individuals
Outputting VCF file...
After filtering, kept 45749 out of a possible 5720701 Sites

vcftools --vcf MUS.LR75.raw.withmissing.24.recode.vcf --minQ 25 --min-meanDP 6 --max-missing 0.42 --recode --recode-INFO-all --out FILTERED/MUS.LR75.24.filtered.withmissing

After filtering, kept 24 out of 24 Individuals
Outputting VCF file...
After filtering, kept 5778 out of a possible 5071269 Sites


vcftools --vcf MUS.LR75.raw.withmissing.vcf --minQ 25 --min-meanDP 6 --max-missing 0.2 --recode --recode-INFO-all --out FILTERED/MUS.LR75.48.filtered.withmissing


After filtering, kept 48 out of 48 Individuals
Outputting VCF file...
After filtering, kept 3527 out of a possible 5071269 Sites

```

Run pixy
```

```


Mean diversity estimates for MODE and MODC are similar. But 10 fold less for MUS. 
```
#MODE
#How many windows with no data: 

awk -F "\t" '{sum+=$5} END {print sum/(621-0)}' MODE.filtered_pi.txt 
0.00468793


#MODC
#Windows with no data
awk -F "\t" '{print $5}' MODC.filtered_pi.txt |grep NA |wc -l
546

awk -F "\t" '{sum+=$5} END {print sum/(617-546)}' MODC.filtered_pi.txt 
0.0071652

#MUS
#Windows with no data
awk -F "\t" '{print $5}' MUS.filtered_pi.txt |grep NA |wc -l
570

#10kb windows = 614, 24 indivs
awk -F "\t" '{sum+=$5} END {print sum/(614-570)}' MUS.filtered_pi.txt 


#50kb windows = 124 (wc -l MUS.filtered_pi.txt)
awk -F "\t" '{sum+=$5} END {print sum/(124-89)}' MUS.filtered_pi.txt 
0.0117212

```


Copy to mac and plot
```
cd /Users/alexjvr/2020.postdoc/Velocity/E3/ANGSD_FINAL/DiversityEstimates/PIXY
scp bzzjrb@bluecrystalp3.acrc.bris.ac.uk:/newhome/bzzjrb/*txt .

#####
## R
#####

library(ggplot2)
library(dplyr)
library(reshape2)

MODE <- read.table("MODE.filtered_pi.txt", head=T)
MODC <- read.table("MODC.filtered_pi.txt", head=T)

#Plot pixy pi for modern populations
MOD <- bind_rows(MODC, MODE)
ggplot(MOD, aes(x=window_pos_1, y=avg_pi, color=pop))+geom_line()+scale_color_manual(values=c("#2E8B57", "#46CC7C"))+ggtitle("Pixy pi: Modern, unfiltered vcf (meanDP10x)")

#Plot difference in pi for modern populations
MOD2 <- left_join(MODC, MODE, by="window_pos_1", suffix=c(".c", ".e"))
MOD2$diff <- MOD2$avg_pi.c-MOD2$avg_pi.e
ggplot(MOD2, aes(x=window_pos_1, y=diff))+geom_line()+ggtitle("Pixy pi: Difference in modern pi (meanDP10x)")


MUS <- read.table("MUS.filtered_pi.txt", header=T)
POP3 <- bind_rows(MOD, MUS)

#There's a weird MUS point with pi of 1 and some with 0 that probably aren't true. Second plot removes these
ggplot(POP3, aes(x=window_pos_1, y=avg_pi, color=pop))+geom_line()+scale_color_manual(values=c("#2E8B57", "#46CC7C", "#DAA520"))+ggtitle("Pixy pi: meanDP10x only")
ggplot(POP3[which(POP3$avg_pi<1&POP3$avg_pi>0),], aes(x=window_pos_1, y=avg_pi, color=pop))+geom_line()+scale_color_manual(values=c("#2E8B57", "#46CC7C", "#DAA520"))+ggtitle("Pixy pi: meanDP10x only")

```

![alt_txt][pixy2]

[pixy2]:https://user-images.githubusercontent.com/12142475/99698697-fdfbf980-2a88-11eb-92d1-c0dc9adc1355.png
