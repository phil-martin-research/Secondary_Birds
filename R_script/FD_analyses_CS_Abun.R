##################################################################################
#script to analyse changes in functional diversity in secondary tropical forests##
##################################################################################

rm(list=ls())
library(FD)
library(RCurl)
library(plyr)
library(reshape)

#traits are defined as each food type as separate binary trait

#sites - abundance matrix

# species - contains duplicate rows of species and traits in binary form
species<-read.csv("C:/Users/Catherine/Dropbox/Phil & Catherine - Secondary forests/Analysis - Dec 2013/Feeding guilds 11.12.13 v3.csv",header=T)
head(species)

#create species by site matrix, species as columns, sites as row names
Sp_melt<-melt(species,id.vars =c("Scientific.name","SiteID"))
head(Sp_melt)
cast(Sp_melt,value~Scientific.name,)

# traits2 - contains only unique rows (1785 species) and traits in binary form
traits2<-subset(species, !duplicated(Taxon.ID))
names(traits2)
rownames(traits2)<-gsub(" ", ".", traits2$Scientific.name)#set TaxonID as row names

#select just columns that have trait data in
traits3<-traits2[,c(4:10)] # EDIT - 4:7 in original - data in 4:10
traits3<-traits3[order(row.names(traits3)),]

# EDIT - remove Sc and P diet groups
traits4<-traits2[,c(5:8,10)]
traits4<-traits4[order(row.names(traits4)),]

#create abundance matrix with site as a row and species as a column

#sites1 - abundance data
sites1<-read.csv("C:/Users/Catherine/Dropbox/Phil & Catherine - Secondary forests/Analysis - Dec 2013/Site_Abun2.csv",header=T)
head(sites1)
rownames(sites1)<-sites1[,1] #set SiteID as row names
sites1[,1]<-NULL #remove SiteID column
#reorder these so that the columns are in ascending alphabetical order 
sites1<-sites1[,order(names(sites1))]
#and do the same for traits
# traits3<-traits3[,order(names(traits3))] # EDIT - removed as done above?

#convert sites1 from list to numeric data matrix
sites1<-data.matrix(sites1)
mode(sites1)

test<-dbFD(traits4,sites1) # without Sc and P dietary groups
# test1<-dbFD(traits3,sites2) - with all dietary groups - error

setwd("C:/Users/Catherine/Dropbox/Phil & Catherine - Secondary forests/Analysis - Dec 2013/R Outputs")
write.csv(test$FRic,"Abun_FRic.csv")
write.csv(test$FEve,"Abun_FEve.csv")
write.csv(test$FDiv,"Abun_FDiv.csv")
write.csv(test$nbsp,"Abun_nbsp.csv")
write.csv(test$sing.sp,"Abun_singsp.csv")
write.csv(test$qual.FRic,"Abun_qual.csv")
