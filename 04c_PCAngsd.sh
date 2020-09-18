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

###MODE

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


###MODC

###MUS
imiss_MUS.GL.minDP20 <- MUS.GL.minDP20
imiss_MUS.GL.minDP20[imiss_MUS.GL.minDP20=="0.333333"] <- NA

mean(is.na(imiss_MUS.GL.minDP20)) ##calculates the overall proportion of missingness
[1] 0.4453347

MUS.loci.propNA <- rowMeans(is.na(imiss_MUS.GL.minDP20))*100  ##calculates the proportion of NA in each row, i.e. for each locus

#Remove these loci
melt.MUS.loci.propNA <- melt(MUS.loci.propNA)
melt.MUS.loci.propNA$marker <- MUS.GL.minDP20$marker 
dim(melt.MUS.loci.propNA[which(melt.MUS.loci.propNA$value>20),])  #loci with >20% missingness
[1] 2042187       2    ##This is basically all the loci! Let's plot the missingness to see if some indivs are the problem.

MUS.markerstoremove <- (melt.MUS.loci.propNA[which(melt.MUS.loci.propNA$value>20),])
MUS.GL.minDP20.clean <- MUS.GL.minDP20[which(!MUS.GL.minDP20$marker %in% MUS.markerstoremove$marker),]  #remove problematic loci
dim(MUS.GL.minDP20.clean)



#Check which indivs have >20% missingness in the clean dataset
imiss_MUS.GL.minDP20 <- MUS.GL.minDP20.clean
imiss_MUS.GL.minDP20[imiss_MUS.GL.minDP20=="0.333333"] <- NA
MUS.indiv.propNA <- colMeans(is.na(imiss_MUS.GL.minDP20))*100  ##calculates the proportion of NA in each column, i.e. for each individual
melt.MUS.indiv.propNA <- melt(MUS.indiv.propNA)                ##plotting is easier with long data
melt.MUS.indiv.propNA$Indiv <- rownames(melt.MUS.indiv.propNA)
melt.MUS.indiv.propNA.new <- melt.MUS.indiv.propNA[seq(4, nrow(melt.MUS.indiv.propNA),3),]  ##get rid of the first 3 lines, and keep only one row per indiv.
melt.MUS.indiv.propNA.new[which(melt.MUS.indiv.propNA.new$value>20),]

[1] value Indiv
<0 rows> (or 0-length row.names)
```
