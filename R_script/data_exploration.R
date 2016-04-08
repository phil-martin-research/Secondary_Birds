#this script is for analysis of differences in species richness and functional diversity metrics
#between primary and secondary forest using metrics accounting for species abundance

library(FD)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(GGally)
library(reshape)
library(lme4)
library(MuMIn)
library(ncf)
library(nlme)

#load in data
rm(list = ls())
FD_comp<-read.csv("Data/FD_abun_summary_comp.csv")
Location<-read.csv("Data/Study_location.csv")
FD_comp<-merge(FD_comp,Location,by="SiteID")

#do data exploration following Zuur et al. 2010

################################################
#1 look for outliers in x and y variables#######
################################################

#subset data to give only data to be used as y or x variables
FD_comp_sub<-FD_comp[,c(1,9:19)]
FD_comp_sub$row_number<-as.numeric(rownames(FD_comp_sub))

FD_comp_sub_melt<-melt(FD_comp_sub,id.vars=c("SiteID","row_number"))

ggplot(FD_comp_sub_melt,aes(x=variable,y=value))+geom_boxplot()
ggplot(FD_comp_sub_melt,aes(x=value,y=row_number))+geom_point()+facet_wrap(~variable,scales = "free")

#loop to compare observed data to hypothetical data with the same mean and sd
Data_test_summary<-NULL
for (i in 9:19){
Mean_val<-mean(FD_comp[,i])
SD_val<-sd(FD_comp[,i])
Data_test<-data.frame(variable=names(FD_comp[i]),
          data=c("observed","simulated"),
           row_number=rep(as.numeric(rownames(FD_comp_sub)),2),
           values=c(FD_comp[,i],
          rnorm(n = 43,mean =Mean_val,sd = SD_val)))
Data_test_summary<-rbind(Data_test_summary,Data_test)
}

ggplot(Data_test_summary,aes(y=row_number,x=values,colour=data))+geom_point(shape=1)+facet_wrap(~variable,scale="free")
ggplot(Data_test_summary,aes(x=values,fill=data))+geom_density(alpha=0.5)+facet_wrap(~variable,scale="free")

#data looks ok - one possible outlier for FDiv

mean(FD_comp$Age)
sd(FD_comp$Age)
Data_test_age<-data.frame(variable="Age",
                      data=c("observed","simulated"),
                      row_number=rep(as.numeric(rownames(FD_comp_sub)),2),
                      values=c(FD_comp$Age,
                               rpois(n = 43,lambda=18.66279)))
ggplot(Data_test_age,aes(y=row_number,x=values,colour=data))+geom_point(shape=1)
#doing the same for age shows that there are two outliers
#this suggests that the x axis should be log transformed to deal with this problem


#######################################
#2 - look at homogeneity of variance###
#for y variables#######################
#######################################

median(FD_comp$Age)#median Age is around 11, so we can chop data into two pieces to check homogeniety of variance
FD_comp$bin<-ifelse(FD_comp$Age>11,">11","<11")
FD_comp_homo<-FD_comp[,c(1,9:19,ncol(FD_comp))]
FD_comp_homo_melt<-melt(FD_comp_homo,id.vars=c("SiteID","bin"))
ggplot(FD_comp_homo_melt,aes(x=value,fill=bin))+geom_histogram()+facet_wrap(~variable,scales = "free")
#broadly speaking variances between the two groups look similar



#########################################
#3 - normal districution of y variables##
#########################################
#########################################

ggplot(FD_comp_homo_melt,aes(x=value))+geom_histogram()+facet_wrap(~variable,scales = "free")
ggplot(FD_comp_homo_melt,aes(x=value))+geom_density()+facet_wrap(~variable,scales = "free")
#some signs of skew in the data, but mostly fairly normally distributed

########################################
#4- zeros in data#######################
#data contains very few zeros###########
########################################

########################################
#5 - check colliniearity of variables###
########################################
FD_comp$logAge<-log(FD_comp$Age)

ggpairs(FD_comp[,c(ncol(FD_comp),9:19)],
        lower=list(continuous="smooth", params=c(colour="black",alpha=0.2)))

#there is a lot of co-linearity between some y variables
#however, this is not a big issue as we are really only concerned
#with colinearity between x variables

#########################################
#6 - are response variables observations#
#independent?############################
#########################################

#for this work we need to check if the data are spatially correlated

FD_comp$Lat<-FD_comp$Lat+(rnorm(length(FD_comp$Lat),0,0.00001)) 
#carry out analyses to see if there is any evidence of spatial autocorrelation for effect sizes
spat_cor_summary<-NULL
for (i in seq(grep("SpR", (colnames(FD_comp))),grep("Rao", colnames(FD_comp)))){
spat_cor<-spline.correlog(FD_comp$Lat, FD_comp$Long, FD_comp[,i], latlon=T,resamp=1000,quiet=T,na.rm=T)
spat_cor2<-data.frame(Dist=spat_cor$boot$boot.summary$predicted$x[1,],
                      Cor=spat_cor$boot$boot.summary$predicted$y[6,],
                      UCI=spat_cor$boot$boot.summary$predicted$y[2,],
                      LCI=spat_cor$boot$boot.summary$predicted$y[10,],
                      variable=names(FD_comp)[i])
spat_cor_summary<-rbind(spat_cor_summary,spat_cor2)
print(i)
}

theme_set(theme_bw(base_size=12))
Moran_plot1<-ggplot(spat_cor_summary,aes(x=Dist,y=Cor,ymax=LCI,ymin=UCI))+geom_ribbon(alpha=0.2)+geom_line(size=0.5,colour="black")+facet_wrap(~variable,scale="free")
Moran_plot2<-Moran_plot1+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))+ theme(legend.position="none")
Moran_plot2+geom_hline(y=0,lty=2)+xlab("Distance between sites (km)")+ylab("Moran's I correlation")+scale_size_continuous(range = c(1,3))
ggsave("Figures/Spatial_autocorrelation.pdf",width = 16,height=8,dpi=400,units="in")

#there appears to be some spatial autocorrelation for evenness and shannon diversity indicies, but little for anything else

