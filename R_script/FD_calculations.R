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

#create a list of all the studies we are using
Unique_study<-unique(P_ab$Study)

#first calculate presence/absence statistics
#for the moment this is only using Petchy and GAston's method
#it might eb useful to come up with another method to calculate functional diversity
FD_summary<-NULL
Unique_study<-unique(P_ab$Study)
for (i in 1:length(Unique_study)){
  PA_sub<-subset(P_ab,Study==Unique_study[i])
  Sites<-unique(PA_sub$Site.ID)
  PF_SF<-PA_sub$PF_SF
  Age<-PA_sub$Age
  PA_sub2<-PA_sub[-c(1:4)]#remove columns that indicate site and study number
  keeps<-colSums(PA_sub2)>0#get rid of species which are not present in local species pool
  PA_sub3<-PA_sub2[,keeps]
  PA_sub3<-PA_sub3[ , order(names(PA_sub3))]#order columns for species alphabetically
  ncol(PA_sub3)
  names(PA_sub3)
  Trait_sp2<-data.frame(Species=row.names(Traits3))#create a dataframe with one column containing species names for which we have traits
  Trait_sp2$Match<-row.names(Traits3) %in% names(PA_sub3)#mark species as "TRUE" if we have details of them in sites and "FALSE" if we do not
  remove_sp2<-subset(Trait_sp2,Match=="FALSE")[,1]#produce vector of species to remove from dataset
  Trait_ab2<-Traits3[-which(rownames(Traits3) %in% remove_sp2), ]#remove species from trait dataset
  Trait_ab2<-Trait_ab2[order(rownames(Trait_ab2)), ]#order trait dataset so it has the same order as the species dataset
  FD_dendro_summary<-FD_dendro(S=Trait_ab2, A=PA_sub3,Cluster.method = "average", ord = "podani",Weigthedby = "abundance")
  FD_summary_study<-dbFD(Trait_ab2, PA_sub3, corr="sqrt",w = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,#produce fd metrics, giving all four traits a similar weight
                  1/7,1/7,1/7,1/7,1/7,1/7,1/7,1,1),w.abun = FALSE,calc.FRic=T,m=10)
  FD_site<-data.frame(Site=Sites,Study=Unique_study[i],PF_SF=PF_SF,Age=Age,SpR=FD_dendro_summary$n_sp,FDpg=FD_dendro_summary$FDpg,
                      FRic=FD_summary_study$FRic,qual_FRic=FD_summary_study$qual.FRic,FEve=FD_summary_study$FEve,
                      FDiv=FD_summary_study$FDiv,FDis=FD_summary_study$FDis,RaoQ=FD_summary_study$RaoQ)
  FD_summary<-rbind(FD_site,FD_summary)
  print(i)
}

FD_summary_melt<-melt(FD_summary,id.vars=c("Site","Study","SpR","PF_SF"))
write.csv(FD_summary,"Data/FD_summary.csv",row.names=F)


ggplot(FD_summary_melt,aes(x=PF_SF,y=value,colour=PF_SF))+geom_point()+facet_wrap(~variable,scales = "free")+geom_boxplot()


#now put data into sites
SF_FD<-subset(FD_summary,PF_SF=="SF")
PF_FD<-subset(FD_summary,PF_SF=="PF")
FD_comp<-merge(SF_FD,PF_FD,by="Study")
colnames(FD_comp)
FD_comp$SpR_comp<-log(FD_comp$SpR.x)-log(FD_comp$SpR.y)
FD_comp$FR_comp<-log(FD_comp$FRic.x)-log(FD_comp$FRic.y)
FD_comp$FE_comp<-log(FD_comp$FEve.x)-log(FD_comp$FEve.y)
FD_comp$FDiv_comp<-log(FD_comp$FDiv.x)-log(FD_comp$FDiv.y)
FD_comp$FDis_comp<-log(FD_comp$FDis.x)-log(FD_comp$FDis.y)
FD_comp$Rao_comp<-log(FD_comp$RaoQ.x)-log(FD_comp$RaoQ.y)


colnames(FD_comp2)<-c("StudyID","Age","Cont_F","Point_obs","Mist_nets","Transect","Vocal","No_Methods","nspb_P","FRic_P","FEve_P","FDiv_P",
                      "FDis_P","SpR_comp","FR_comp","FE_comp","FDiv_comp","FDis_comp")

write.csv(FD_comp2,"Data/FD_comp.csv",row.names = F)


#maybe I can try to calculate the abundance based metrics?

#now calculate abundance statistics
#for this I need to produce FD statistics for one site at a time
#probably best to use a loop to achieve this

Abun2<-subset(Abun,SiteID!=3&SiteID!=4&SiteID!=11&SiteID!=12&SiteID!=24&
                SiteID!=25&SiteID!=26&SiteID!=36&SiteID!=47&SiteID!=48&
                SiteID!=53&SiteID!=54&SiteID!=64&SiteID!=65&SiteID!=37&
                SiteID!=5&SiteID!=6&SiteID!=11&SiteID!=10)#subset to remove sites with incomplete data
Abun3<-Abun2[-3]#remove column with species codes
Abun3$Species<-gsub(" ", ".", Abun3$Species, fixed = TRUE)#put dot in between species and genus name
head(Abun3)
#for each unique site run this loop once
Abun3<-subset(Abun3,Species!="#N/A")
FD_summary_abun<-NULL
Unique_study<-unique(Abun3$Study)
for (i in 1:length(Unique_study)){
  i<-1
  Abun_sub<-subset(Abun3,Study==Unique_site[i])
  Study<-unique(Abun_sub$Study)
  Abun_sub$Abundance<-Abun_sub$Abundance*100
  Abun_sub2<-spread(Abun_sub,Species,Abundance)#spread data so that each species has a column
  Abun_sub3<-Abun_sub2[,-c(1,2)]#remove site and study number columns
  Trait_sp2<-data.frame(Species=row.names(Traits3))#create a dataframe with one column containing species names for which we have traits
  Trait_sp2$Match<-row.names(Traits3) %in% names(Abun_sub3)#mark species as "TRUE" if we have details of them in sites and "FALSE" if we do not
  remove_sp2<-subset(Trait_sp2,Match=="FALSE")[,1]#produce vector of species to remove from dataset
  Trait_ab2<-Traits3[-which(rownames(Traits3) %in% remove_sp2), ]#remove species from trait dataset
  Trait_ab2<-Trait_ab2[order(rownames(Trait_ab2)), ]#order trait dataset so it has the same order as the species dataset
  
  FD_dendro_summary<-FD_dendro(S=Trait_ab2, A=PA_sub3,Cluster.method = "average", ord = "podani",Weigthedby = "abundance")
  FD_summary_study<-dbFD(Trait_ab2, PA_sub3, corr="sqrt",w = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,#produce fd metrics, giving all four traits a similar weight
                                                               1/7,1/7,1/7,1/7,1/7,1/7,1/7,1,1),w.abun = T,calc.FRic=T,m=10)
  FD_site<-data.frame(Site=Sites,Study=Unique_study[i],PF_SF=PF_SF,Age=Age,SpR=FD_dendro_summary$n_sp,FDpg=FD_dendro_summary$FDpg,
                      FRic=FD_summary_study$FRic,qual_FRic=FD_summary_study$qual.FRic,FEve=FD_summary_study$FEve,
                      1/7,1/7,1/7,1/7,1/7,1/7,1/7,1,1))
  
  FD_dendro(S=Trait_ab2, A=Abun_sub3,Cluster.method = "average", ord = "podani",Weigthedby = "abundance")
  
  FD_site<-data.frame(Site=Unique_site[i],Study=Study,SpR=FD_calc$nbsp,FRic=FD_calc$FRic,FEve=FD_calc$FEve,FDiv=FD_calc$FDiv,FDis=FD_calc$FDis)
  FD_summary<-rbind(FD_site,FD_summary)
  print(i)
}