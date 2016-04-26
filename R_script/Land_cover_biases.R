#calculate mean cover of different land-uses in buffer around site points
library(raster)
library(ggplot2)
library(dismo)
library(plyr)
library(reshape)

Studies<-read.csv("Data/Study_location2.csv")
Studies2<-as.matrix(cbind(Studies$Longitude,Studies$Latitude))


#standard error function
stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

Tree<-raster("Data/Remote_sensing/Concensus_cover/Tree_cover_tropics.tif")
Cult<-raster("Data/Remote_sensing/Concensus_cover/Cult_cover_tropics.tif")
Urb<-raster("Data/Remote_sensing/Concensus_cover/Urb_cover_tropics.tif")


#extract information about tree cover
Studies$Tree<-extract(Tree,Studies2,buffer=2500,fun=mean,na.rm=T)
Studies$Cult<-extract(Cult,Studies2,buffer=2500,fun=mean,na.rm=T)
Studies$Urb<-extract(Urb,Studies2,buffer=2500,fun=mean,na.rm=T)

ggplot(Studies,aes(x=Tree))+geom_histogram()
ggplot(Studies,aes(x=Urb))+geom_histogram()
ggplot(Studies,aes(x=Cult))+geom_histogram()

#create 44 random points and extract data on tree cover, cultivation and urban areas
#first create a forest mask using 40% tree cover as a definition of forest
Forest<-calc(Tree,fun=function(x){ifelse(x>40,1,NA)})

#run the loop 10 times and calculate the median mean value for the cover around each of the points
Cov_rand_summary<-NULL
for (i in 1:1000){
  Rand_points<-randomPoints(Forest,n=44,tryf = 15)
  Rand_mat<-as.matrix(cbind(Rand_points[,1],Rand_points[,2]))
  Tree_rand<-extract(Tree,Rand_mat,buffer=2500,fun=mean,na.rm=T)
  Cult_rand<-extract(Cult,Rand_mat,buffer=2500,fun=mean,na.rm=T)
  Urb_rand<-extract(Urb,Rand_mat,buffer=2500,fun=mean,na.rm=T)
  Cov_rand<-data.frame(run=i,M_Tree=mean(Tree_rand,na.rm=T),M_Cult=mean(Cult_rand,na.rm=T),M_other=100-(mean(Cult_rand,na.rm=T)+mean(Tree_rand,na.rm=T)))
  Cov_rand_summary<-rbind(Cov_rand,Cov_rand_summary)
  print(i)
}

Cov_melt<-melt(Cov_rand_summary,id.vars="run")

Cover_quantile<-ddply(Cov_melt, .(variable), function(x) quantile(x$value,c(.025, .5, .975)))
names(Cover_quantile)<-c("variable","Lower","Median","Upper")
write.csv(Cover_quantile,"Output/Cover_quantiles.csv")

