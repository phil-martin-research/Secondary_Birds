#this script calculates species richness and functional diversity metrics needed for paper

#notes - need to calculate changes in evenness and data on individual traits (e.g. Body size)

rm(list=ls())

library(FD)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(GGally)

#load in data
P_ab<-read.csv("Data/Site_PresAb.csv")
Abun<-read.csv("Data/Site_Abun3.csv")
Traits<-read.csv("Data/Bird_traits.csv")

#first calculate presence/absence statistics
P_ab2<-P_ab[-1]#remove column indicate site number
row.names(Traits)<-gsub(" ", ".", Traits$Scientific.Name, fixed = TRUE)#put dot in between species and genus name
Traits2<-Traits[-c(1:4)]#remove columns that are not needed from trait file
Traits3<-data.matrix(Traits2)#convert trait data to a data.matrix
Trait_sp<-data.frame(Species=row.names(Traits3))#create a dataframe with one column containing species names for which we have traits
Trait_sp$Match<-row.names(Traits3) %in% names(P_ab2)#mark species as "TRUE" if we have details of them in sites and "FALSE" if we do not
remove_sp<-subset(Trait_sp,Match=="FALSE")[,1]#produce vector of species to remove from dataset
Traits4<-Traits3[-which(rownames(Traits3) %in% remove_sp), ]#remove species from trait dataset
FD_calc<-dbFD(Traits4, P_ab2, corr="cailliez",w = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,#produce fd metrics, giving all four traits a similar weight
                                                    1/7,1/7,1/7,1/7,1/7,1/7,1/7,1,1))

#now put data into sites
FD_site<-data.frame(Site.ID=P_ab$Site.ID,FD_calc)
Studies<-read.csv("Data/Studies.csv")
FD_site2<-merge(FD_site,Studies,by.x="Site.ID",by.y="SiteID")
SF_FD<-subset(FD_site2,PF_SF=="SF")
PF_FD<-subset(FD_site2,PF_SF=="PF")
FD_comp<-merge(SF_FD,PF_FD,by="StudyID")
colnames(FD_comp)
FD_comp$SpR_comp<-log(FD_comp$nbsp.x)-log(FD_comp$nbsp.y)
FD_comp$FR_comp<-log(FD_comp$FRic.x)-log(FD_comp$FRic.y)
FD_comp$FE_comp<-log(FD_comp$FEve.x)-log(FD_comp$FEve.y)
FD_comp$FDiv_comp<-log(FD_comp$FDiv.x)-log(FD_comp$FDiv.y)
FD_comp$FDis_comp<-log(FD_comp$FDis.x)-log(FD_comp$FDis.y)


FD_comp2<-(FD_comp[-c(2:29,31,37,39:40,42,44,48:77)])
colnames(FD_comp2)<-c("StudyID","Age","Cont_F","Point_obs","Mist_nets","Transect","Vocal","No_Methods","nspb_P","FRic_P","FEve_P","FDiv_P",
                      "FDis_P","SpR_comp","FR_comp","FE_comp","FDiv_comp","FDis_comp")

write.csv(FD_comp2,"Data/FD_comp.csv",row.names = F)


#now calculate abundance statistics
Abun2<-subset(Abun,SiteID!=3&SiteID!=4&SiteID!=11&SiteID!=12&SiteID!=24&
                SiteID!=25&SiteID!=26&SiteID!=36&SiteID!=47&SiteID!=48&
                SiteID!=53&SiteID!=54&SiteID!=64&SiteID!=65&SiteID!=37&
                SiteID!=5&SiteID!=6&SiteID!=11&SiteID!=10)#subset to remove sites with incomplete data
Abun3<-Abun2[-2]#remove column with species codes
head(Abun3)
Abun3$Species<-gsub(" ", ".", Abun3$Species, fixed = TRUE)#put dot in between species and genus name
Abun4<-spread(Abun3,Species,Abundance)#spread data so that each species has a column
Abun5<-Abun4[-c(1:2)]#remove site number and #NA column

Trait_sp2<-data.frame(Species=row.names(Traits3))#create a dataframe with one column containing species names for which we have traits
Trait_sp2$Match<-row.names(Traits3) %in% names(Abun5)#mark species as "TRUE" if we have details of them in sites and "FALSE" if we do not
remove_sp2<-subset(Trait_sp2,Match=="FALSE")[,1]#produce vector of species to remove from dataset
Trait_ab2<-Traits3[-which(rownames(Traits3) %in% remove_sp2), ]#remove species from trait dataset
Trait_ab2<-Trait_ab2[order(rownames(Trait_ab2)), ]#order trait dataset so it has the smae order as the species dataset
#produce fd metrics using relative abundances, giving all four traits a similar weight
FD_calc<-dbFD(Trait_ab2, Abun5, corr="cailliez",w = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
                                                    1/7,1/7,1/7,1/7,1/7,1/7,1/7,1,1))

#now put data into sites
FD_site<-data.frame(Site.ID=P_ab$Site.ID,FD_calc)
Studies<-read.csv("Data/Studies.csv")
FD_site2<-merge(FD_site,Studies,by.x="Site.ID",by.y="SiteID")
SF_FD<-subset(FD_site2,PF_SF=="SF")
PF_FD<-subset(FD_site2,PF_SF=="PF")
FD_comp<-merge(SF_FD,PF_FD,by="StudyID")
colnames(FD_comp)
FD_comp$SpR_comp<-log(FD_comp$nbsp.x)-log(FD_comp$nbsp.y)
FD_comp$FR_comp<-log(FD_comp$FRic.x)-log(FD_comp$FRic.y)
FD_comp$FE_comp<-log(FD_comp$FEve.x)-log(FD_comp$FEve.y)
FD_comp$FDiv_comp<-log(FD_comp$FDiv.x)-log(FD_comp$FDiv.y)
FD_comp$FDis_comp<-log(FD_comp$FDis.x)-log(FD_comp$FDis.y)

ggplot(FD_comp,aes(x=nbsp.y,y=SpR_comp))+geom_point()+geom_smooth(method="lm",group=NA)
ggplot(FD_comp,aes(x=Age.x,y=FR_comp))+geom_point()+geom_smooth(method="lm",group=NA)
ggplot(FD_comp,aes(x=Age.x,y=FE_comp))+geom_point()+geom_smooth(method="lm",group=NA)
ggplot(FD_comp,aes(x=Age.x,y=FDiv_comp))+geom_point()+geom_smooth(method="lm",group=NA)
ggplot(FD_comp,aes(x=Age.x,y=FDis_comp))+geom_point()+geom_smooth(method="lm",group=NA)

colnames(FD_comp)

FD_comp2<-(FD_comp[-c(2:29,31,37,39:40,42,44,48:77)])
colnames(FD_comp2)<-c("StudyID","Age","Cont_F","Point_obs","Mist_nets","Transect","Vocal","No_Methods","nspb_P","FRic_P","FEve_P","FDiv_P",
                      "FDis_P","SpR_comp","FR_comp","FE_comp","FDiv_comp","FDis_comp")
