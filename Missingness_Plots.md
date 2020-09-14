# Missingness in GL files

Missing data in ANGSD GL files are coded as equal likelihoods (0.3/0.3/0.3). 

i.e. for individuals and loci for which there is no data, instead of missing data (e.g. -9 or 0 in vcf files) the data is encoded as equal likelihoods. 

Missingness will affect things like estimates of LD, estimating the SAF and SFS, and the PCA, so it's important to know how patchy your data is. 


1. Estimate GLs (e.g. -GL 2) and write to beagle format (-doGlf 2). See [ANGSD documentation](http://www.popgen.dk/angsd/index.php/Genotype_Likelihoods#Output)

2. Read into R and draw heatmaps for each population

```
R
#install.packages("varhandle")
#install.packages("gplots")
#install.packages("RColorBrewer")
library(varhandle)
library(gplots)
library(RColorBrewer)

###MODE
################
MODE.GL <- read.table(gzfile("MODE.LR761675.1.1.GLF2.PVAL.MIMNAF0.05.MININD18.beagle.gz"), header=T)
bloc.MODE <- MODE.GL[,4:ncol(MODE.GL)] #The first three cols are 1) Chr name and position, 2) allele 1, 3) allele 2
bloc3.MODE <- as.matrix(unfactor(bloc.MODE))
#bloc4 <- bloc3[1:5,] #test subset of data on plot
bloc4.MODE <- bloc3.MODE[1:500,] # We can look at just a subset of the loci to get an idea of the patchiness of the data

#Change all the extra Ind labels to NA
labCol.MODE <- colnames(bloc3.MODE)
labCol.MODE <- gsub("Ind\\d{1}\\.1", "NA", labCol.MODE)
labCol.MODE <- gsub("Ind\\d{2}\\.1", "NA", labCol.MODE)
labCol.MODE <- gsub("Ind\\d{2}\\.2", "NA", labCol.MODE)
labCol.MODE <- gsub("Ind\\d{1}\\.2", "NA", labCol.MODE)
labCol.MODE[labCol.MODE=="NA"] <- NA

p1.MODE <- heatmap.2(bloc4.MODE, trace="none", col=brewer.pal(11,"RdBu"),scale="none", Colv=F, Rowv=T, dendrogram="none", labCol=labCol.MODE, adjCol=c(0.5,1))
#scale - the GLs are already rescaled to 1
#Colv=F to keep the column order
#Colv=T to check for loci with more missing data. We've used minInd 18, so this should look okay. 
#labCol and adjCol: labels for the x-axis. We have 3 cols per indiv, so I've used a single label, and moved it to be in the middle of the three cols. 

####MODC
#############
MODC.GL <- read.table(gzfile("SAM/MODC.LR761675.1.1.GLF2.PVAL0.05.MINMAF0.05.MININD18.DOMAF1.beagle.gz"), header=T)
bloc.MODC <- MODC.GL[,4:ncol(MODC.GL)] #The first three cols are 1) Chr name and position, 2) allele 1, 3) allele 2
bloc3.MODC <- as.matrix(unfactor(bloc.MODC))
#bloc4 <- bloc3[1:5,] #test subset of data on plot
bloc4.MODC <- bloc3.MODC[1:500,] # We can look at just a subset of the loci to get an idea of the patchiness of the data

#Change all the extra Ind labels to NA
labCol.MODC <- colnames(bloc3.MODC)
labCol.MODC <- gsub("Ind\\d{1}\\.1", "NA", labCol.MODC)
labCol.MODC <- gsub("Ind\\d{2}\\.1", "NA", labCol.MODC)
labCol.MODC <- gsub("Ind\\d{2}\\.2", "NA", labCol.MODC)
labCol.MODC <- gsub("Ind\\d{1}\\.2", "NA", labCol.MODC)
labCol.MODC[labCol.MODC=="NA"] <- NA

p2.MODC <- heatmap.2(bloc4.MODC, trace="none", col=brewer.pal(11,"RdBu"),scale="none", Colv=F, Rowv=T, dendrogram="none", labCol=labCol.MODC, adjCol=c(0.5,1))

####MUS
###########
MUS.GL <- read.table(gzfile("MUS.LR761675.1.1.GLF2.PVAL.MIMNAF0.05.MININD18.beagle.gz"), header=T)
bloc.MUS <- MUS.GL[,4:ncol(MUS.GL)] #The first three cols are 1) Chr name and position, 2) allele 1, 3) allele 2
bloc3.MUS <- as.matrix(unfactor(bloc.MUS))
#bloc4 <- bloc3[1:5,] #test subset of data on plot
bloc4.MUS <- bloc3.MUS[1:500,] # We can look at just a subset of the loci to get an idea of the patchiness of the data

#Change all the extra Ind labels to NA
labCol.MUS <- colnames(bloc3.MUS)
labCol.MUS <- gsub("Ind\\d{1}\\.1", "NA", labCol.MUS)
labCol.MUS <- gsub("Ind\\d{2}\\.1", "NA", labCol.MUS)
labCol.MUS <- gsub("Ind\\d{2}\\.2", "NA", labCol.MUS)
labCol.MUS <- gsub("Ind\\d{1}\\.2", "NA", labCol.MUS)
labCol.MUS[labCol.MUS=="NA"] <- NA

p3.MUS <- heatmap.2(bloc4.MUS, trace="none", col=brewer.pal(11,"RdBu"),scale="none", Colv=F, Rowv=T, dendrogram="none", labCol=labCol.MUS, adjCol=c(0.5,1))



##Figures

pdf("MODE.GLheatmap.pdf")
heatmap.2(bloc4.MODE, trace="none", col=brewer.pal(11,"RdBu"),scale="none", Colv=F, Rowv=T, dendrogram="none", labCol=labCol.MODE, adjCol=c(0.5,1))
dev.off()

pdf("MODC.GLheatmap.pdf")
p2.MODC
dev.off()

pdf("MUS.GLheatmap.pdf")
p3.MUS
dev.off()

```


![alt_txt][heatmap]

[heatmap]:https://user-images.githubusercontent.com/12142475/93075083-347a5a00-f67d-11ea-8614-0e04e1f48835.png




We can see that, in order, 1) MODE has a few individuals with very little coverage that we need to remove, but overall very dense coverage. 

2) MODC the patchiness looks even between loci and individuals (apart from a couple indivs), but overall less dense genotyping than in MODE. 

3) MUS has quite sparse data, and clearly quite a few indivs that have genotyped very poorly. 

It's evident that we need to filter on missingness within individuals. We can also see that 0.3 GL is synonymous with missing data. So we can quantify the proportion of missing data within an individual by calculating the proportion of 0.3 GL within one col of each individual. 

```
library(reshape2)

####MODE
###########

imiss_MODE.GL <- MODE.GL
imiss_MODE.GL[imiss_MODE.GL=="0.333333"] <- NA

mean(is.na(imiss_MODE.GL)) ##calculates the overall proportion of missingness
#[1] 0.2755835

MODE.propNA <- colMeans(is.na(imiss_MODE.GL))*100  ##calculates the proportion of NA in each column, i.e. for each individual

##Then we can visualize the missingness
melt.MODE.propNA <- melt(MODE.propNA) #plotting is easier with long data
melt.MODE.propNA$Ind <- rownames(melt.MODE.propNA)
df.new.MODE <- melt.MODE.propNA[seq(4, nrow(melt.MODE.propNA),3),]  ##get rid of the first 3 lines, and keep only one row per indiv.
#df.new.MODE <- df.new.MODE[order(df.new.MODE$value)]
p4.MODE <- hist(df.new.MODE$value, breaks=seq(0,100,10))  ##we can see that only a few individuals have more than 50% missing data
df.new.MODE[df.new.MODE$value>50,]  ##Identify these individuals
         value   Ind
Ind18 56.74097 Ind18
Ind32 99.34785 Ind32
Ind34 99.51919 Ind34
Ind35 99.65731 Ind35
Ind31 99.66780 Ind31
Ind36 99.68354 Ind36
Ind33 99.79544 Ind33


####MODC
###########

imiss_MODC.GL <- MODC.GL
imiss_MODC.GL[imiss_MODC.GL=="0.333333"] <- NA

mean(is.na(imiss_MODC.GL)) ##calculates the overall proportion of missingness
#[1] 0.3642511

MODC.propNA <- colMeans(is.na(imiss_MODC.GL))*100  ##calculates the proportion of NA in each column, i.e. for each individual

##Then we can visualize the missingness
melt.MODC.propNA <- melt(MODC.propNA) #plotting is easier with long data
melt.MODC.propNA$Ind <- rownames(melt.MODC.propNA)
df.new.MODC <- melt.MODC.propNA[seq(4, nrow(melt.MODC.propNA),3),]  ##get rid of the first 3 lines, and keep only one row per indiv.
#df.new.MODC <- df.new.MODC[order(df.new.MODC$value)]
p5.MODC <- hist(df.new.MODC$value, breaks=seq(0,100,10))  ##we can see that only a few individuals have more than 50% missing data
df.new.MODC[df.new.MODC$value>50,]  ##Identify these individuals
         value   Ind
Ind9  63.03054  Ind9
Ind19 55.15271 Ind19

####MUS
###########

imiss_MUS.GL <- MUS.GL
imiss_MUS.GL[imiss_MUS.GL=="0.333333"] <- NA

mean(is.na(imiss_MUS.GL)) ##calculates the overall proportion of missingness
#[1] 0.2755835

MUS.propNA <- colMeans(is.na(imiss_MUS.GL))*100  ##calculates the proportion of NA in each column, i.e. for each individual

##Then we can visualize the missingness
melt.MUS.propNA <- melt(MUS.propNA) #plotting is easier with long data
melt.MUS.propNA$Ind <- rownames(melt.MUS.propNA)
df.new.MUS <- melt.MUS.propNA[seq(4, nrow(melt.MUS.propNA),3),]  ##get rid of the first 3 lines, and keep only one row per indiv.
#df.new.MUS <- df.new.MUS[order(df.new.MUS$value)]
p6.MUS <- hist(df.new.MUS$value, breaks=seq(0,100,10))  ##we can see that only a few individuals have more than 50% missing data
df.new.MUS[df.new.MUS$value>50,]  ##Identify these individuals
         value   Ind
Ind1  81.49237  Ind1
Ind2  78.29443  Ind2
Ind4  55.15530  Ind4
Ind6  86.43632  Ind6
Ind11 75.02297 Ind11
Ind16 82.66863 Ind16
Ind17 72.06396 Ind17
Ind18 72.72560 Ind18
Ind25 96.52637 Ind25
Ind28 63.90369 Ind28
Ind29 93.21816 Ind29
Ind30 84.59842 Ind30
Ind31 64.65723 Ind31
Ind35 64.95130 Ind35
Ind36 51.62654 Ind36
Ind37 56.69914 Ind37
Ind40 57.72836 Ind40
Ind41 73.14832 Ind41
Ind42 84.81897 Ind42
Ind43 86.43632 Ind43
Ind44 60.59548 Ind44
Ind45 65.09833 Ind45
Ind46 66.91785 Ind46
Ind47 79.48906 Ind47


par(mfrow=c(1,3))
plot(p4.MODE, xlim=c(0,100), ylim=c(0,20))
plot(p5.MODC,xlim=c(0,100), ylim=c(0,20))
plot(p6.MUS,xlim=c(0,100), ylim=c(0,20))


```

![alt_txt][Missing]

[Missing]:https://user-images.githubusercontent.com/12142475/92750803-558e3280-f37f-11ea-8d0f-4f4c4e2d73af.png


We also want to see all the sites where heterozygotes are over- or under-represented. 

