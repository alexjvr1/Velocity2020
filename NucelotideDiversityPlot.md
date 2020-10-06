# Plot Nucleotide diversity across the genome

See Fig 3 in [Feng et al. 2019](https://reader.elsevier.com/reader/sd/pii/S0960982218316099?token=65EE80B75A634FD527312349DBECB9E11DFF80648448692DC1E7DA07DD7E7A6552DFB5DD8B92B7420E4E045F526A4075#go_to_%22%C3%BE%C3%BF\u0000f\u0000l\u0000i\u0000n\u0000k\u00004%22)


Thetas calculated in ANGSD in windows (-win 50kb -step 10kb)

## Plots in R

Copy diversity estiates to mac
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


Tajima's D

![alt_txt][TajD]

[TajD]:https://user-images.githubusercontent.com/12142475/95207692-0998a700-07e0-11eb-820c-3961f79ecac3.png

