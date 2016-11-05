#this script is to calculate the standard effect size of FD as per the reviewer's comments

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
library(picante)

#load in data
Traits<-read.csv("Data/Bird_traits.csv")
Abun<-read.csv("Data/Site_abun5.csv")

#Since FD is known to be correlated with species richness, we also calculated sesFD, 
#which adjusts FD for species richness by comparing the observed FD at each point to 999 
#scenarios in which the number of species at each point is held constant but species identity 
#is randomly drawn from the community. Communities that are more functionally diverse than would 
#be expected by chance (given their species richness) have positive values of sesFD, and those which are 
#less functionally diverse than expected by chance have nega- tive values

#tidy data
Abun3<-Abun[-9]#remove column with species codes
Abun3$Species<-gsub(" ", ".", Abun3$Species, fixed = TRUE)#put dot in between species and genus name
#now sort out trait data
row.names(Traits)<-gsub(" ", ".", Traits$Scientific.Name, fixed = TRUE)#put dot in between species and genus name
Traits2<-Traits[-c(1:4)]#remove columns that are not needed from trait file
Traits3<-data.matrix(Traits2)#convert trait data to a data.matrix
Abun3<-subset(Abun3,Species!="#N/A")
Abun3<-Abun3[complete.cases(Abun3),]
FD_site_SF_Summary<-NULL
FD_site_PF_Summary<-NULL

#create a loop in which the number of species is constant for each loop
#but the species identity varies at random

Unique_study<-unique(Abun3$Study)

i<-1

Abun_sub<-subset(Abun3,Study==Unique_study[i])#subset data so that it is only from one study
Study_info<-unique(Abun_sub[,1:8])#Store info on study ID
Abun_sub2<-Abun_sub[-c(3:8,11)]#remove columns that indicate site, study number, whether they are primary or secondary, and methods used in study
sum_abundance<-sum(Abun_sub2$Abundance)

for (i in 1:1000){
  Abun_sub2$Abundance2<-0.01
  for 
  sum(Abun_sub2$Abundance2)
}


Abun_sub2<-spread(Abun_sub2,Species,Abundance)#spread data so that each species has a column
keeps<-colSums(Abun_sub2,na.rm = T)>0#identify species not present in local species pool
Abun_sub3<-Abun_sub2[,keeps]#get rid of species which are not present in local species pool
Abun_sub4<-Abun_sub3[,-c(1:2)]#remove site and study IDs
ses.pd(Abun_sub4)
Abun_sub4[is.na(Abun_sub4)] <- 0#replace any NA values for abundance with zero
Trait_sp2<-data.frame(Species=row.names(Traits3))#create a dataframe with one column containing species names for which we have traits
Trait_sp2$Match<-row.names(Traits3) %in% names(Abun_sub3)
