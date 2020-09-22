# PCAngsd

Population structure using PCAngsd

PCAngsd can deal with low coverage data (2x), but is affected by missingness. 

Input = Beagle GLs

Instead of calling a combined GL for all populations I can join the individual GLs by matching the first three columns of the beagle file (marker, allele1, allele2). I.e. private alleles will be excluded here. 

I can also filter for missingness in this way. 


I've created three GL datasets with minInd 10 and minDepth 20X (i.e. loci are kept if they have at least 2x coverage in 10 individuals within a population). 
```
#on mac
/Users/alexjvr/2020.postdoc/Velocity/E3/ANGSD_FINAL/SFS

R
MODE.GL.minDP20 <- read.table(gzfile("MODE.LR761675.1.minDP20.MinIND10.beagle.gz"), header=T) 
dim(MODE.GL.minDP20)
[1] 5261260     102

MODC.GL.minDP20 <- read.table(gzfile("MODC.LR761675.1.minDP20.MinIND10.beagle.gz"), header=T) 
dim(MODC.GL.minDP20)
[1] 5204985     111

MUS.GL.minDP20 <- read.table(gzfile("MUS.LR761675.1.minDP20.MinIND10.beagle.gz"), header=T) 
dim(MUS.GL.minDP20)
[1] 2093813      75

#######################################################################################
#Identify loci within pop that have <80% genotyping rate and remove them
#######################################################################################

library(reshape2)

#See below for each indiv
```

### MODE
```
imiss_MODE.GL.minDP20 <- MODE.GL.minDP20
imiss_MODE.GL.minDP20[imiss_MODE.GL.minDP20=="0.333333"] <- NA

mean(is.na(imiss_MODE.GL.minDP20)) ##calculates the overall proportion of missingness
#[1] 0.09138739

MODE.loci.propNA <- rowMeans(is.na(imiss_MODE.GL.minDP20))*100  ##calculates the proportion of NA in each row, i.e. for each locus

#Remove these loci
melt.MODE.loci.propNA <- melt(MODE.loci.propNA)
melt.MODE.loci.propNA$marker <- MODE.GL.minDP20$marker 
dim(melt.MODE.loci.propNA[which(melt.MODE.loci.propNA$value>20),])  #loci with >20% missingness
[1] 815553      2

MODE.markerstoremove <- (melt.MODE.loci.propNA[which(melt.MODE.loci.propNA$value>20),])
MODE.GL.minDP20.clean <- MODE.GL.minDP20[which(!MODE.GL.minDP20$marker %in% MODE.markerstoremove$marker),]  #remove problematic loci
dim(MODE.GL.minDP20.clean)
[1] 4445707     102


#Check which indivs have >20% missingness in the clean dataset
imiss_MODE.GL.minDP20 <- MODE.GL.minDP20.clean
imiss_MODE.GL.minDP20[imiss_MODE.GL.minDP20=="0.333333"] <- NA
MODE.indiv.propNA <- colMeans(is.na(imiss_MODE.GL.minDP20))*100  ##calculates the proportion of NA in each column, i.e. for each individual
melt.MODE.indiv.propNA <- melt(MODE.indiv.propNA)                ##plotting is easier with long data
melt.MODE.indiv.propNA$Indiv <- rownames(melt.MODE.indiv.propNA)
melt.MODE.indiv.propNA.new <- melt.MODE.indiv.propNA[seq(4, nrow(melt.MODE.indiv.propNA),3),]  ##get rid of the first 3 lines, and keep only one row per indiv.
melt.MODE.indiv.propNA.new[which(melt.MODE.indiv.propNA.new$value>20),]

[1] value Indiv
<0 rows> (or 0-length row.names)

melt.MODE.indiv.propNA.new
           value Indiv
Ind0   1.4480037  Ind0
Ind1   1.4096296  Ind1
Ind2   1.7483158  Ind2
Ind3   2.3698143  Ind3
Ind4   3.4244947  Ind4
Ind5   2.5063280  Ind5
Ind6   1.5945045  Ind6
Ind7   1.8299002  Ind7
Ind8   1.3003781  Ind8
Ind9   4.0188883  Ind9
Ind10  7.1385046 Ind10
Ind11  3.5150090 Ind11
Ind12  1.9820694 Ind12
Ind13  3.2800632 Ind13
Ind14  3.4366862 Ind14
Ind15  3.0831092 Ind15
Ind16  1.6373324 Ind16
Ind17  1.5841125 Ind17
Ind18  4.9873732 Ind18
Ind19 11.8568768 Ind19
Ind20 13.2753013 Ind20
Ind21  1.6282450 Ind21
Ind22  1.6046042 Ind22
Ind23  0.7728804 Ind23
Ind24  3.5716254 Ind24
Ind25  3.3341828 Ind25
Ind26 17.5604915 Ind26
Ind27  6.6139087 Ind27
Ind28  3.6319982 Ind28
Ind29  1.3915897 Ind29
Ind30  4.7006697 Ind30
Ind31  2.0929180 Ind31
Ind32  3.3313262 Ind32
```

### MODC
```
imiss_MODC.GL.minDP20 <- MODC.GL.minDP20
imiss_MODC.GL.minDP20[imiss_MODC.GL.minDP20=="0.333333"] <- NA

mean(is.na(imiss_MODC.GL.minDP20)) ##calculates the overall proportion of missingness
[1] 0.2608227

MODC.loci.propNA <- rowMeans(is.na(imiss_MODC.GL.minDP20))*100  ##calculates the proportion of NA in each row, i.e. for each locus

#Remove these loci
melt.MODC.loci.propNA <- melt(MODC.loci.propNA)
melt.MODC.loci.propNA$marker <- MODC.GL.minDP20$marker 
dim(melt.MODC.loci.propNA[which(melt.MODC.loci.propNA$value>20),])  #loci with >20% missingness
[1] 2982476       2   ##about 50% of the loci

MODC.markerstoremove <- (melt.MODC.loci.propNA[which(melt.MODC.loci.propNA$value>20),])
MODC.GL.minDP20.clean <- MODC.GL.minDP20[which(!MODC.GL.minDP20$marker %in% MODC.markerstoremove$marker),]  #remove problematic loci
dim(MODC.GL.minDP20.clean)
[1] 2222509     111


#Check which indivs have >20% missingness in the clean dataset
imiss_MODC.GL.minDP20 <- MODC.GL.minDP20.clean
imiss_MODC.GL.minDP20[imiss_MODC.GL.minDP20=="0.333333"] <- NA
MODC.indiv.propNA <- colMeans(is.na(imiss_MODC.GL.minDP20))*100  ##calculates the proportion of NA in each column, i.e. for each individual
melt.MODC.indiv.propNA <- melt(MODC.indiv.propNA)                ##plotting is easier with long data
melt.MODC.indiv.propNA$Indiv <- rownames(melt.MODC.indiv.propNA)
melt.MODC.indiv.propNA.new <- melt.MODC.indiv.propNA[seq(4, nrow(melt.MODC.indiv.propNA),3),]  ##get rid of the first 3 lines, and keep only one row per indiv.
melt.MODC.indiv.propNA.new[which(melt.MODC.indiv.propNA.new$value>20),]

         value Indiv   ##two indivs with just over 20% missing data. We'll keep these in for now
Ind27 20.49279 Ind27
Ind29 20.97013 Ind29
```

### MUS
```
imiss_MUS.GL.minDP20 <- MUS.GL.minDP20
imiss_MUS.GL.minDP20[imiss_MUS.GL.minDP20=="0.333333"] <- NA

mean(is.na(imiss_MUS.GL.minDP20)) ##calculates the overall proportion of missingness
[1] 0.4453347

MUS.loci.propNA <- rowMeans(is.na(imiss_MUS.GL.minDP20))*100  ##calculates the proportion of NA in each row, i.e. for each locus

#Remove these loci
melt.MUS.loci.propNA <- melt(MUS.loci.propNA)
melt.MUS.loci.propNA$marker <- MUS.GL.minDP20$marker 
dim(melt.MUS.loci.propNA[which(melt.MUS.loci.propNA$value>20),])  #loci with >20% missingness
[1] 2042187       2    ##This is basically all the loci! But we're just trying to draw a PCA, so 50k SNPs should be enough. We'll have to see how many overlap with MODC and MODE in the final dataset. What will be our min cut-off here? This is just one chr so we could include markers from all chrs? But this will be a computational nightmare when the dataset gets bigger (and is arguably unnecessary for a PCA...)

MUS.markerstoremove <- (melt.MUS.loci.propNA[which(melt.MUS.loci.propNA$value>20),])
MUS.GL.minDP20.clean <- MUS.GL.minDP20[which(!MUS.GL.minDP20$marker %in% MUS.markerstoremove$marker),]  #remove problematic loci
dim(MUS.GL.minDP20.clean)
[1] 51626    75


#Check which indivs have >20% missingness in the clean dataset
imiss_MUS.GL.minDP20 <- MUS.GL.minDP20.clean
imiss_MUS.GL.minDP20[imiss_MUS.GL.minDP20=="0.333333"] <- NA
MUS.indiv.propNA <- colMeans(is.na(imiss_MUS.GL.minDP20))*100  ##calculates the proportion of NA in each column, i.e. for each individual
melt.MUS.indiv.propNA <- melt(MUS.indiv.propNA)                ##plotting is easier with long data
melt.MUS.indiv.propNA$Indiv <- rownames(melt.MUS.indiv.propNA)
melt.MUS.indiv.propNA.new <- melt.MUS.indiv.propNA[seq(4, nrow(melt.MUS.indiv.propNA),3),]  ##get rid of the first 3 lines, and keep only one row per indiv.
melt.MUS.indiv.propNA.new[which(melt.MUS.indiv.propNA.new$value>20),]

[1] value Indiv        ###7 indivs with around 20% missing data. 
Ind1  22.80440  Ind1
Ind5  22.65138  Ind5
Ind7  21.08627  Ind7
Ind11 21.16763 Ind11
Ind13 23.55015 Ind13
Ind22 23.30221 Ind22
Ind23 22.38794 Ind23
```

#### Combine cleaned data together

Create a new beagle file by combining all the cleaned data together. Join by site and allele (ie. if different alleles are present in the different populations we're ignoring these data for now). 

```
library(dplyr)
#library(compare) #to test if dfs have intersected correctly

#Find markers that occur in all three cleaned datasets
MODE.marker <- as.data.frame(MODE.GL.minDP20.clean$marker)
colnames(MODE.marker) <- "marker"

MODC.marker <- as.data.frame(MODC.GL.minDP20.clean$marker)
colnames(MODC.marker) <- "marker"

MUS.marker <- as.data.frame(MUS.GL.minDP20.clean$marker)
colnames(MUS.marker) <- "marker"

POP2clean.markers <- intersect(MUS.marker, MODE.marker)  ##dplyr can only intersect 2 data frames at a time
POP3clean.markers <- intersect(POP2clean.markers, MODC.marker)  

dim(POP3clean.markers)
[1] 29651     1

#POP3clean.markers.test <- intersect(MUS.marker, MODE.marker, MODC.marker)  #not the same dimension as POPO3clean.markers.test, so this doesn't work correctly. 
#compare(POP3clean.markers.test, POP3clean.markers) #allowAll=T to ignore col and row order
#FALSE [FALSE]  

##Subset each dataset to keep only the overlapping markers
MODC.clean.sub  <- MODC.GL.minDP20.clean[which(MODC.GL.minDP20.clean$marker %in% POP3clean.markers$marker),]  
MODE.clean.sub  <- MODE.GL.minDP20.clean[which(MODE.GL.minDP20.clean$marker %in% POP3clean.markers$marker),]
MUS.clean.sub  <- MUS.GL.minDP20.clean[which(MUS.GL.minDP20.clean$marker %in% POP3clean.markers$marker),]

##Join datasets together
pops2.clean <- left_join(MODC.clean.sub, MODE.clean.sub, by="marker", suffix=c(".c", ".e"))
pops3.clean <- left_join(pops2.clean, MUS.clean.sub, by="marker", suffix=c(".mod", ".mus"))

##Rename individuals, but rememeber the population oder (MODC, MODE, MUS)
cols.toremove <- as.data.frame(grep("allele", colnames(pops3.clean)))[-(1:2),] #find all the columns with "allele". We need to remove all but the first two
[1]   112 113 213 214
pops3.clean2  <- pops3.clean[-cols.toremove]

ncol(pops3.clean2) #count the number of columns
[1] 282

#create a new vector with all the colnames. First 3 cols = "marker", "allele1", "allele2". Rest= indivs. The total=(ncol-3)/3
pops3.names <- paste("Ind", rep((0:((ncol(pops3.clean2)-6)/3)), each=3), sep="")  #total indivs= [ncol(df)-6 (first 3 cols+ -3 because Ind name starts at 0)] % 3 (3 GLs for each indiv)
appendix.pop3 <- paste(".",rep((as.character(c("0","1","2"))), times=93), sep="")  #we need to append .1 and .2 to the second and third col for each indiv
pop3.names2 <- paste(pops3.names, appendix.pop3, sep="")  #append 
pop3.names2 <- gsub(".0$", "", pop3.names2) #remove .0 from the first col for each indiv

colnames(pops3.clean2) <- c("marker", "allele1", "allele2", pop3.names2) #rename colnames of the cleaned dataset

colnames(pops3.clean2)  #check that this has worked
  [1] "marker"  "allele1" "allele2" "Ind0"    "Ind0.1"  "Ind0.2"  "Ind1"   
  [8] "Ind1.1"  "Ind1.2"  "Ind2"    "Ind2.1"  "Ind2.2"  "Ind3"    "Ind3.1" 
 [15] "Ind3.2"  "Ind4"    "Ind4.1"  "Ind4.2"  "Ind5"    "Ind5.1"  "Ind5.2" 

dim(pops3.clean2)
[1] 29651   282

write.table(pops3.clean2, "POPS3.clean.beagle", quote=F, sep="\t", row.names=F)

```


# Run PCAdapt

```
/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04b_ANGSD_FINAL/SFS_and_Fst

module load languages/python-anaconda3-5.2.0

python ~/software/pcangsd/pcangsd.py -h

##gzip beagle file
gzip *beagle


##calculate the covariance matrix
python ~/software/pcangsd/pcangsd.py -beagle *beagle.gz -o test1 -threads 10

#download to mac and plot in R

pwd /Users/alexjvr/2020.postdoc/Velocity/E3/ANGSD_FINAL/PCAngsd/PCA.minInd10.minDP20

scp bluecp:/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04b_ANGSD_FINAL/SFS_and_Fst/test1.cov .

#open R
library(reshape2)
library(ggplot2)

C <- as.matrix(read.table("test1.cov"))
e <- eigen(C)

forPCA <- as.data.frame(e$vectors)
forPCA$pop <- paste(rep(c("MODC", "MODE", "MUS"), times=c(36, 33,24)))  #popnames - get numbers and order from beagle files and order in which they were combined above
forPCAmelt <- melt(forPCA, id.vars=c("V1","V2","pop")) #long data for ggplot

e$values #find proportion explained by PC1 and PC2

e$values[1]
[1] 10.456
e$values[2]
[1] 4.618

pdf(file="E3.PCA_0.8genotypingrate.pdf")
ggplot(forPCAmelt, aes(x=V1, y=V2, color=pop)) + geom_point()+ ggtitle("PCA:E3 (80% genotyping rate, 25k loci)") + xlab(paste("PC1: 10.46%"))+ ylab(paste("PC2: 4.62%"))
dev.off()


```
