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
library(ggplot2)
library(lme4)
library(MuMIn)
library(picante)

#load in data
Traits<-read.csv("Data/Bird_traits_Oct_2016.csv")
#remove traits we don't want to use for this analysis
Traits<-Traits[,c(1:14,17:21,23,24,26)]
Abun<-read.csv("Data/Site_abun5.csv")
#subset to remove the studies of de Lima et al and Slik and Van Balen
Abun<-subset(Abun,Study!=28&Study!=15)

#produce a loop to calculate statistics for each study one at a time
#for this I need to calculate the FD statistics as well as Petchy and Gaston's FD
row.names(Traits)<-gsub(" ", ".", Traits$Scientific.Name, fixed = TRUE)#put dot in between species and genus name
Traits2<-Traits[-c(1:4)]#remove columns that are not needed from trait file
Traits3<-data.matrix(Traits2)#convert trait data to a data.matrix

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
  Abun_sub2<-Abun_sub[-c(3:8,11)]#remove columns that indicate site, study number, whether they are primary or secondary, and methods used in study
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
  #calculate standardised effect size FD
  #calculate distance based on trait differences of species
  D<-vegdist(Trait_ab2,"gower")
  #calculate tree based on distances
  tree<-hclust(D,"average")
  #covert this to a phylogeny
  tree.p<-as.phylo(tree)
  SES_FD<-ses.pd(Abun_sub4,tree.p,null.model = "richness",runs=999,iterations = 1000)
  FD_dendro_summary<-FD_dendro(S=Trait_ab2, A=Abun_sub4,Cluster.method = "average", ord = "podani",Weigthedby = "abundance")
  FD_summary_study<-dbFD(Trait_ab2, Abun_sub4, corr="sqrt",w = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,#produce fd metrics, giving all four traits a similar weight
                                                                 1/5,1/5,1/5,1/5,1/5,1,1,1),w.abun = T,calc.FRic=T,m=10)
  
  FD_site<-data.frame(Study_info,SpR=FD_dendro_summary$n_sp,FDpg=FD_dendro_summary$FDpg,FDw=FD_dendro_summary$FDw,
                      FRic=FD_summary_study$FRic,qual_FRic=FD_summary_study$qual.FRic,FEve=FD_summary_study$FEve,
                      FDiv=FD_summary_study$FDiv,FDis=FD_summary_study$FDis,RaoQ=FD_summary_study$RaoQ,FD_summary_study$CWM,FD_SES=SES_FD$pd.obs.z)

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

str(FD_site_PF_Summary)


SES_FD_summary<-merge(x=FD_site_SF_Summary[,c(1:8,ncol(FD_site_SF_Summary))],y=FD_site_PF_Summary[,c(1:8,ncol(FD_site_PF_Summary))],by=c("Study","Point_obs","Mist_nets","Transect","Vocal"))
SES_FD_summary2<-SES_FD_summary[,-c(6,8,10:12)]
names(SES_FD_summary2)<-c("study","Point_obs","Mist_nets","Transect","Vocal","Age","SES_Sec","SES_P")
SES_FD_summary2$SES_Diff<-SES_FD_summary2$SES_Sec-SES_FD_summary2$SES_P

write.csv(SES_FD_summary2,"Data/SES_summary.csv",row.names=F)

SES_FD_summary2<-read.csv("Data/SES_summary.csv")
ggplot(SES_FD_summary2,aes(x=Age,y=SES_Diff))+geom_point()+scale_x_log10()+geom_smooth(method="lm")
ggplot(SES_FD_summary2,aes(x=Age,y=SES_Sec))+geom_point()+scale_x_log10()+geom_smooth(method="lm")+geom_hline(yintercept=mean(SES_FD_summary2$SES_P))


M1<-lmer(SES_Sec~log(Age)+(1|study),data=SES_FD_summary2)
M2<-lmer(SES_Diff~log(Age)+(1|study),data=SES_FD_summary2)

summary(M1)
r.squaredGLMM(M1)
summary(M2)
r.squaredGLMM(M2)


#then calculate proportional change in variables between secondary and primary forests
SF_prop_summary<-NULL
SF_prop_summary2<-NULL
Unique_study<-unique(data.frame(FD_site_SF_Summary$SiteID,FD_site_SF_Summary$Study))
for (i in 1:nrow(Unique_study)){
  SF_sub<-subset(FD_site_SF_Summary,Study==Unique_study[i,2]&SiteID==Unique_study[i,1])
  SF_sub<-(SF_sub[,c(1:16)])
  PF_sub<-subset(FD_site_PF_Summary,Study==Unique_study[i,2])
  PF_sub<-(PF_sub[,c(1:16)])
  SF_prop_sub<-data.frame(SF_sub[,c(1:8)],log(SF_sub[,c(9:ncol(SF_sub))]/PF_sub[,c(9:ncol(PF_sub))]))
  SF_prop_summary<-rbind(SF_prop_sub,SF_prop_summary)
}

write.csv(SF_prop_summary,"Data/FD_abun_summary_comp.csv",row.names=F)
