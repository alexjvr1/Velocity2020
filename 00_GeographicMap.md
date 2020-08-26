# Figure 1

A map of all the sampling locations and a picture of Ringlet. Plus a PCA of the genomic structure. 

Velocity samples and coordinates are curated on Dropbox: 20Species_SampleSites_AJvR.xlsx

Ringlet data is curated in the Ringlet folder on Dropbox: /Velocity WP1 analysis/RingletMS_2020/Ringlet_Tables.xlsx

On Mac

pwd /Users/alexjvr/2020.postdoc/Velocity/E3/GeographicMap 

Copy all the sample details from the [Library Prep](https://docs.google.com/spreadsheets/d/1G9r50W0VV_ANZ19rIvqZpXWFemy2MW76_iXuyBuCQGA/edit#gid=1076895649) shared doc: 

Batch convert all the 100m grid references to Lat,Long [here](https://gridreferencefinder.com/batchConvert/batchConvert.php) 


Final site info input: E3.Ringlet_SiteInfo
```
Population	Site	Gridref_100m	Lat	Long	Pop	Pop2	pch	Colour
MODC	Pitton	SU205311	51.07886	-1.7087396	Modern	Mod.C	22	seagreen4
MODC	CaverelCopse	SU200310	51.077978	-1.7158825	Modern	Mod.C	22	seagreen4
MODC	CaverelCopse	SU199310	51.077982	-1.71731	Modern	Mod.C	22	seagreen4
MODC	CaverelCopse	SU195308	51.076197	-1.7230305	Modern	Mod.C	22	seagreen4
MODC	CaverelCopse	SU191310	51.078009	-1.7287296	Modern	Mod.C	22	seagreen4
MODC	MartinDownNNR	SU036201	50.980297	-1.9500924	Modern	Mod.C	22	seagreen4
MODC	MartinDownNNR	SU036203	50.982096	-1.9500905	Modern	Mod.C	22	seagreen4
MODE	Grantown	NJ033263	57.316975	-3.6071842	Modern	Mod.E	22	seagreen4
MODE	Dunphail	NJ015481	57.512331	-3.6458164	Modern	Mod.E	22	seagreen4
MUS	Museum	SU2602308620	50.87649918	-1.631500006	Museum	Museum	22	goldenrod
MUS	Museum	SU2987408208	50.87260056	-1.576799989	Museum	Museum	22	goldenrod
MUS	Museum	SU2998602236	50.81890106	-1.575700045	Museum	Museum	22	goldenrod
MUS	Museum	SU3887021934	50.99549866	-1.44749999	Museum	Museum	22	goldenrod
MUS	Museum	SU4197212417	50.90969849	-1.404399991	Museum	Museum	22	goldenrod

```


### Draw map in R
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

##E3 data points
E3.SiteInfo <- read.table("E3.Ringlet_SiteInfo", header=T)  ##read table of info into R
UK.coords <- E3.SiteInfo ##change name of file so that it can be manipulated
coordinates(UK.coords) <- c("Long", "Lat")  # set spatial coordinates. Remember to check that these make sense
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
proj4string(UK.coords) <- crs.geo  # define projection system of our data
UK.coords$Colour <- as.character(UK.coords$Colour) ##to plot in colour
summary(UK.coords)

##add columns for colour and shape (mus, mod.e, mod.c). I've already done this in the excel sheet

#UK.coords$pch <- c(21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21)

##plot map
pdf("E3.sampleMap_20200826.pdf")
plot(hill, col = grey(0:100/100), legend = FALSE, main = "UK")
plot(UK.coords, pch=UK.coords$pch, cex=2, col="black", bg=UK.coords$Colour, add=T)
dev.off()

```


### Supp Fig: Map of Core samples


Zoom in on the core sample sites to show that they're not all collected in exactly the same place. 

```
library(maps)
library(mapdata)
UK.coords.Core <- UK.coords[which(UK.coords$Pop!="Mod.E"),] #subset the data to exclude expanding pop

pdf("E3.ZoomedInCoreMap.pdf")
map("worldHires","UK", xlim=c(-2,-1), ylim=c(50.5,51.5), col="gray90", fill=TRUE)
points(UK.coords.Core$Long, UK.coords.Core$Lat, pch=UK.coords.Core$pch, col="black", bg=UK.coords.Core$Colour, cex=2)
dev.off()
```
