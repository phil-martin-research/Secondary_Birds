#analysis of biases 


#load packages
library(ggplot2)
library("ggmap")
library(maptools)
library(maps)
library(raster)
library(rgdal)

#load in data on tree cover
Trees_1<-raster("Data/Remote_sensing/Concensus_cover/consensus_full_class_1.tif")
Trees_2<-raster("Data/Remote_sensing/Concensus_cover/consensus_full_class_2.tif")
Trees_3<-raster("Data/Remote_sensing/Concensus_cover/consensus_full_class_3.tif")
Trees_4<-raster("Data/Remote_sensing/Concensus_cover/consensus_full_class_4.tif")
Tree_cover<-sum(Trees_1,Trees_2,Trees_3,Trees_4)
plot(Tree_cover)

#crop data so that it represents tropics only
Tree_cover_tropics1<-crop(Trees_1,extent(-180,180,-30,30))
Tree_cover_tropics2<-crop(Trees_2,extent(-180,180,-30,30))
Tree_cover_tropics3<-crop(Trees_3,extent(-180,180,-30,30))
Tree_cover_tropics4<-crop(Trees_4,extent(-180,180,-30,30))

#stick rasters together to get an estimate of total tree cover
Tree_cover_tropics<-sum(Tree_cover_tropics1,Tree_cover_tropics2,Tree_cover_tropics3,Tree_cover_tropics4)

#save raster
writeRaster(x=Tree_cover_tropics, filename="Data/Remote_sensing/Concensus_cover/Tree_cover_tropics.tif", format="GTiff", overwrite=TRUE)


#load data on cultivated land
Cult<-raster("Data/Remote_sensing/Concensus_cover/consensus_full_class_7.tif")
#crop data so that it represents tropics only
Cult_cover_tropics<-crop(Cult,extent(-180,180,-30,30))
writeRaster(x=Cult_cover_tropics, filename="Data/Remote_sensing/Concensus_cover/Cult_cover_tropics.tif", format="GTiff", overwrite=TRUE)

#load data on urban land
Urb<-raster("Data/Remote_sensing/Concensus_cover/consensus_full_class_9.tif")
#crop data so that it represents tropics only
Urb_cover_tropics<-crop(Urb,extent(-180,180,-30,30))
writeRaster(x=Urb_cover_tropics, filename="Data/Remote_sensing/Concensus_cover/Urb_cover_tropics.tif", format="GTiff", overwrite=TRUE)
