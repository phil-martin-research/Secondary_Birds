#this script calculates species richness and functional diversity metrics needed for paper

#notes - need to calculate changes in evenness and data on individual traits (e.g. Body size)

rm(list=ls())

library(FD)
library(fundiv)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(GGally)
library(reshape)

#load in data
P_ab<-read.csv("Data/Site_PresAb2.csv")
Traits<-read.csv("Data/Bird_traits.csv")

#first calculate presence/absence statistics

#produce a loop to calculate statistics for each study one at a time
#for this I need to calculate the FD statistics as well as Petchy and Gaston's FD

row.names(Traits)<-gsub(" ", ".", Traits$Scientific.Name, fixed = TRUE)#put dot in between species and genus name
Traits2<-Traits[-c(1:4)]#remove columns that are not needed from trait file
Traits3<-data.matrix(Traits2)#convert trait data to a data.matrix

#creat a list of all the studies we are using
Unique_study<-unique(P_ab$Study)

#first calculate presence/absence statistics
i<-1
P_ab2<-P_ab[-c(1,2)]#remove column indicate site number
ncol(P_ab)
head(P_ab)[1:3]
PA2<-subset(P_ab,Study==Unique_study[i])
PA2_2<-PA2[-c(1,2)]#remove column indicate site number
keeps<-colSums(PA2_2)>0#get rid of species which are not present in local species pool
PA2_3<-PA2_2[,keeps]
PA2_3<-PA2_3[ , order(names(PA2_3))]
PA2_4<-PA2_3[,-c(1:20)]

names(PA2_4)

row.names(Traits)<-gsub(" ", ".", Traits$Scientific.Name, fixed = TRUE)#put dot in between species and genus name
Traits2<-Traits[-c(1:4)]#remove columns that are not needed from trait file
Traits3<-data.matrix(Traits2)#convert trait data to a data.matrix
Trait_sp<-data.frame(Species=row.names(Traits3))#create a dataframe with one column containing species names for which we have traits
Trait_sp$Match<-row.names(Traits3) %in% names(PA2_4)#mark species as "TRUE" if we have details of them in sites and "FALSE" if we do not
remove_sp<-subset(Trait_sp,Match=="FALSE")[,1]#produce vector of species to remove from dataset
Traits4<-Traits3[-which(rownames(Traits3) %in% remove_sp), ]#remove species from trait dataset
Traits4<-Traits4[order(rownames(Traits4)), ]
FD_calc<-dbFD(Traits4, PA2_4, corr="cailliez",w = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,#produce fd metrics, giving all four traits a similar weight
                                                    1/7,1/7,1/7,1/7,1/7,1/7,1/7,1,1))
Unique_study<-unique(P_ab$Study)

for (i in 1:length(Unique_site)){
  i<-1
  PA_sub<-subset(P_ab,Study==Unique_study[i])
  Sites<-unique(PA_sub$Site.ID)
  Study<-Unique_study[i]
  PA_sub2<-PA_sub[-c(1:2)]#remove columns that indicate site and study number
  keeps<-colSums(PA_sub2)>0#get rid of species which are not present in local species pool
  PA_sub3<-PA_sub2[,keeps]
  PA_sub3<-PA_sub3[ , order(names(PA_sub3))]#order columns for species alphabetically
  Trait_sp2<-data.frame(Species=row.names(Traits3))#create a dataframe with one column containing species names for which we have traits
  Trait_sp2$Match<-row.names(Traits3) %in% names(PA_sub3)#mark species as "TRUE" if we have details of them in sites and "FALSE" if we do not
  remove_sp2<-subset(Trait_sp2,Match=="FALSE")[,1]#produce vector of species to remove from dataset
  Trait_ab2<-Traits3[-which(rownames(Traits3) %in% remove_sp2), ]#remove species from trait dataset
  Trait_ab2<-Trait_ab2[order(rownames(Trait_ab2)), ]#order trait dataset so it has the same order as the species dataset
  FD_calc<-dbFD(Trait_ab2, PA_sub3, corr="cailliez",w = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
                                                            1/7,1/7,1/7,1/7,1/7,1/7,1/7,1,1))
  FD_dendro_summary<-FD_dendro(S=Trait_ab2, A=PA_sub3,Cluster.method = "average", ord = "podani",Weigthedby = "abundance")
  
  FD_site<-data.frame(Site=Unique_site[i],Study=Study,SpR=FD_calc$nbsp,FRic=FD_calc$FRic,
                      FEve=FD_calc$FEve,FDiv=FD_calc$FDiv,FDis=FD_calc$FDis,
                      FDpg=FD_dendro_summary$FDpg,FDw=FD_dendro_summary$FDw,
                      FDwcomm=FD_dendro_summary$FDwcomm,qual.FD=FD_dendro_summary$qual.FD)
  FD_summary<-rbind(FD_site,FD_summary)
  print(i)
}


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

colnames(FD_comp2)<-c("StudyID","Age","Cont_F","Point_obs","Mist_nets","Transect","Vocal","No_Methods","nspb_P","FRic_P","FEve_P","FDiv_P",
                      "FDis_P","SpR_comp","FR_comp","FE_comp","FDiv_comp","FDis_comp")

write.csv(FD_comp2,"Data/FD_comp.csv",row.names = F)