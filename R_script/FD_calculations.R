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
library(vegan)

#load in data
P_ab<-read.csv("Data/Site_PresAb2.csv")
Traits<-read.csv("Data/Bird_traits.csv")
Abun<-read.csv("Data/Site_abun5.csv")

#first calculate presence/absence statistics

#produce a loop to calculate statistics for each study one at a time
#for this I need to calculate the FD statistics as well as Petchy and Gaston's FD
row.names(Traits)<-gsub(" ", ".", Traits$Scientific.Name, fixed = TRUE)#put dot in between species and genus name
Traits2<-Traits[-c(1:4)]#remove columns that are not needed from trait file
Traits3<-data.matrix(Traits2)#convert trait data to a data.matrix

#create a list of all the studies we are using
Unique_study<-unique(P_ab$Study)

#first calculate presence/absence statistics
#for the moment this is only using Petchy and GAston's method
#it might eb useful to come up with another method to calculate functional diversity
FD_summary<-NULL
Unique_study<-unique(P_ab$Study)
for (i in 1:length(Unique_study)){
  PA_sub<-subset(P_ab,Study==Unique_study[i])#subset data so that it only represents one study
  Sites<-unique(PA_sub$Site.ID)#get site ID numbers
  PF_SF<-PA_sub$PF_SF#Get info on whether sites are secondary or primary forest
  Age<-PA_sub$Age#get info on site ages
  Methods<-PA_sub[,5:8]#get info on the methods used in studies
  PA_sub2<-PA_sub[-c(1:8)]#remove columns that indicate site, study number, whether they are primary or secondary, and methods used in study
  keeps<-colSums(PA_sub2)>0#identify species not present in local species pool
  PA_sub3<-PA_sub2[,keeps]#get rid of species which are not present in local species pool
  PA_sub3<-PA_sub3[ , order(names(PA_sub3))]#order columns for species alphabetically
  Trait_sp2<-data.frame(Species=row.names(Traits3))#create a dataframe with one column containing species names for which we have traits
  Trait_sp2$Match<-row.names(Traits3) %in% names(PA_sub3)#mark species as "TRUE" if we have details of them in sites and "FALSE" if we do not
  remove_sp2<-subset(Trait_sp2,Match=="FALSE")[,1]#produce vector of species to remove from dataset
  Trait_ab2<-Traits3[-which(rownames(Traits3) %in% remove_sp2), ]#remove species from trait dataset
  Trait_ab2<-Trait_ab2[order(rownames(Trait_ab2)), ]#order trait dataset so it has the same order as the species dataset
  FD_dendro_summary<-FD_dendro(S=Trait_ab2, A=PA_sub3,Cluster.method = "average", ord = "podani",Weigthedby = "abundance")
  FD_summary_study<-dbFD(Trait_ab2, PA_sub3, corr="sqrt",w = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,#produce fd metrics, giving all four traits a similar weight
                  1/7,1/7,1/7,1/7,1/7,1/7,1/7,1,1),w.abun = FALSE,calc.FRic=T,m=10)
  FD_site<-data.frame(Site=Sites,Study=Unique_study[i],PF_SF=PF_SF,Age=Age,Methods,SpR=FD_dendro_summary$n_sp,FDpg=FD_dendro_summary$FDpg,
                      FRic=FD_summary_study$FRic,qual_FRic=FD_summary_study$qual.FRic,FEve=FD_summary_study$FEve,
                      FDiv=FD_summary_study$FDiv,FDis=FD_summary_study$FDis,RaoQ=FD_summary_study$RaoQ)
  FD_site_PF<-subset(FD_site,PF_SF=="PF")#subset to give only primary forest sites
  FD_site_SF<-subset(FD_site,PF_SF=="SF")#subset to give only secondary forest sites
  #include diversity metrics for primary forest reference sites
  FD_site_SF$PF_SpR<-FD_site_PF$SpR
  FD_site_SF$PF_FDpg<-FD_site_PF$FDpg
  FD_site_SF$PF_FRic<-FD_site_PF$FRic
  FD_site_SF$PF_FEve<-FD_site_PF$FEve
  FD_site_SF$PF_FDiv<-FD_site_PF$FDiv
  FD_site_SF$PF_FDis<-FD_site_PF$FDis
  FD_site_SF$PF_RaoQ<-FD_site_PF$RaoQ
  #now calculate the log response ratio effect size as a measure of difference between secondary and primary sites
  FD_site_SF$SpR_comp<-log(FD_site_SF$SpR)-log(FD_site_PF$SpR)
  FD_site_SF$FDpg_comp<-log(FD_site_SF$FDpg)-log(FD_site_PF$FDpg)
  FD_site_SF$FR_comp<-log(FD_site_SF$FRic)-log(FD_site_PF$FRic)
  FD_site_SF$FE_comp<-log(FD_site_SF$FEve)-log(FD_site_PF$FEve)
  FD_site_SF$FDiv_comp<-log(FD_site_SF$FDiv)-log(FD_site_PF$FDiv)
  FD_site_SF$FDis_comp<-log(FD_site_SF$FDis)-log(FD_site_PF$FDis)
  FD_site_SF$Rao_comp<-log(FD_site_SF$RaoQ)-log(FD_site_PF$RaoQ)
  FD_summary<-rbind(FD_site_SF,FD_summary)
  print(i)
}

write.csv(FD_summary,"Data/FD_summary_comp.csv",row.names=F)


#now calculate statistics for data including species abundance
#for this I need to produce FD statistics for one site at a time
#probably best to use a loop to achieve this

Abun3<-Abun[-9]#remove column with species codes
Abun3$Species<-gsub(" ", ".", Abun3$Species, fixed = TRUE)#put dot in between species and genus name

#now sort out trait data
row.names(Traits)<-gsub(" ", ".", Traits$Scientific.Name, fixed = TRUE)#put dot in between species and genus name
Traits2<-Traits[-c(1:4)]#remove columns that are not needed from trait file
Traits3<-data.matrix(Traits2)#convert trait data to a data.matrix


#for each unique site run this loop once
#remove sites NA values for species abundance
Abun3<-subset(Abun3,Species!="#N/A")
Abun3<-Abun3[complete.cases(Abun3),]
FD_site_SF_Summary<-NULL
FD_site_PF_Summary<-NULL
Unique_study<-unique(Abun3$Study)
for (i in 1:length(Unique_study)){
  Abun_sub<-subset(Abun3,Study==Unique_study[i])#subset data so that it is only from one study
  Study_info<-unique(Abun_sub[,1:8])#Store info on study ID
  Abun_sub2<-Abun_sub[-c(3:8)]#remove columns that indicate site, study number, whether they are primary or secondary, and methods used in study
  Abun_sub2$Abundance<-Abun_sub2$Abundance*100
  Abun_sub2<-spread(Abun_sub2,Species,Abundance)#spread data so that each species has a column
  keeps<-colSums(Abun_sub2,na.rm = T)>0#identify species not present in local species pool
  Abun_sub3<-Abun_sub2[,keeps]#get rid of species which are not present in local species pool
  Abun_sub4<-Abun_sub3[,-c(1:2)]#remove site and study IDs
  Abun_sub4[is.na(Abun_sub4)] <- 0#replace any NA values for abundance with zero
  Trait_sp2<-data.frame(Species=row.names(Traits3))#create a dataframe with one column containing species names for which we have traits
  Trait_sp2$Match<-row.names(Traits3) %in% names(Abun_sub3)#mark species as "TRUE" if we have details of them in sites and "FALSE" if we do not
  remove_sp2<-subset(Trait_sp2,Match=="FALSE")[,1]#produce vector of species to remove from dataset
  Trait_ab2<-Traits3[-which(rownames(Traits3) %in% remove_sp2), ]#remove species from trait dataset
  Trait_ab2<-Trait_ab2[order(rownames(Trait_ab2)), ]#order trait dataset so it has the same order as the species dataset
  
  FD_dendro_summary<-FD_dendro(S=Trait_ab2, A=Abun_sub4,Cluster.method = "average", ord = "podani",Weigthedby = "abundance")
  FD_summary_study<-dbFD(Trait_ab2, Abun_sub4, corr="sqrt",w = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,#produce fd metrics, giving all four traits a similar weight
                                                               1/7,1/7,1/7,1/7,1/7,1/7,1/7,1,1),w.abun = T,calc.FRic=T,m=10)
  #Shannon diversity
  Shan_div<-diversity(Abun_sub4,index="shannon")
  #Pielou's evenness
  Even<-Shan_div/log(specnumber(Abun_sub4))
  FD_site<-data.frame(Study_info,SpR=FD_dendro_summary$n_sp,Shan_div,Even,FDpg=FD_dendro_summary$FDpg,FDw=FD_dendro_summary$FDw,
                      FRic=FD_summary_study$FRic,qual_FRic=FD_summary_study$qual.FRic,FEve=FD_summary_study$FEve,
                      FDiv=FD_summary_study$FDiv,FDis=FD_summary_study$FDis,RaoQ=FD_summary_study$RaoQ,FD_summary_study$CWM)
  #convert all values to numeric in new dataframe
  for (y in 9:ncol(FD_site)){
    FD_site[,y]<-as.numeric(as.character(FD_site[,y]))
  }
  
  FD_site_PF<-subset(FD_site,PF_SF=="PF")#subset to give only primary forest sites
  FD_site_SF<-subset(FD_site,PF_SF=="SF")#subset to give only secondary forest sites
  
  FD_site_PF_Summary<-rbind(FD_site_PF,FD_site_PF_Summary)
  FD_site_SF_Summary<-rbind(FD_site_SF,FD_site_SF_Summary)
  print(i)
}


#then calculate proportional change in variables between seconday and primary forests
SF_prop_summary<-NULL
Unique_study<-unique(data.frame(FD_site_SF_Summary$SiteID,FD_site_SF_Summary$Study))
for (i in 1:nrow(Unique_study)){
  SF_sub<-subset(FD_site_SF_Summary,Study==Unique_study[i,2]&SiteID==Unique_study[i,1])
  PF_sub<-subset(FD_site_PF_Summary,Study==Unique_study[i,2])
  SF_prop_sub<-data.frame(SF_sub[,c(1:8)],(SF_sub[,c(9:ncol(SF_sub))]-PF_sub[,c(9:ncol(PF_sub))])/PF_sub[,c(9:ncol(PF_sub))])
  SF_prop_summary<-rbind(SF_prop_sub,SF_prop_summary)
}


SF_prop_check<-SF_prop_summary[,c(1:20)]
SF_prop_check$row_number<-as.numeric(rownames(SF_prop_check))
SF_prop_check2<-subset(SF_prop_check,FRic<400)
SF_prop_check_melt<-melt(SF_prop_check2,id.vars=c("SiteID","Study","Age","PF_SF","Point_obs","Mist_nets","Transect","Vocal","row_number"))
SF_prop_check_melt$value2<-plogis(((SF_prop_check_melt$value+1)/7))

plogis(6.5/7)
ggplot(SF_prop_check_melt,aes(y=row_number,x=value2))+geom_point()+facet_wrap(~variable,scales = "free")





write.csv(FD_summary_abun,"Data/FD_abun_summary_comp.csv",row.names=F)
