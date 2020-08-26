# Figure 1

A map of all the sampling locations and a picture of Ringlet. Plus a PCA of the genomic structure. 

pwd /Users/alexjvr/2018.postdoc/BrownArgus_2018/BrownArgusDataMap 

```
##R

library(sp)  # classes for spatial data
library(raster)  # grids, rasters
library(rasterVis)  # raster visualisation
library(maptools)
library(rgeos)


##get UK map

require(spatial.tools)
elevation<-getData("alt", country = "GB")
x <- terrain(elevation, opt = c("slope", "aspect"), unit = "degrees")
plot(x)
slope <- terrain(elevation, opt = "slope")
aspect <- terrain(elevation, opt = "aspect")
hill <- hillShade(slope, aspect, 40, 270)
plot(hill, col = grey(0:100/100), legend = FALSE, main = "UK")
#plot(elevation, col = rainbow(25, alpha = 0.35), add = TRUE)

##BA data points
BA.SiteInfo <- read.table("BA_SiteInfo", header=T)  ##read table of info into R
UK.coords <- BA.SiteInfo ##change name of file so that it can be manipulated
coordinates(UK.coords) <- c("Long", "Lat")  # set spatial coordinates. Remember to check that these make sense
#plot(CH_coords)
#plot(CH_coords, pch = 20, cex = 1, col = "black")
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
proj4string(UK.coords) <- crs.geo  # define projection system of our data
summary(UK.coords)

##add columns for colour (host plant) and shape (new vs old sites). 

UK.coords$pch <- c(17, 15, 17,15, 15, 15, 17, 15, 17)
UK.coords$Colours <- c("gold1", "gold1", "darkorchid3", "gold1", "darkorchid3", "gold1", "darkorchid3", "gold1", "darkorchid3")

##plot map
pdf("BA.sampleMap_20180619.pdf")
plot(hill, col = grey(0:100/100), legend = FALSE, main = "UK")
plot(UK.coords, pch=UK.coords$pch, cex=2, col=UK.coords$Colours, add=T)
dev.off()

```
