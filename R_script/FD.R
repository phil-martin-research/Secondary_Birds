rm(list=ls())
library(FD)
dummy$trait #dataframe of 8 functional traits for 8 species
dummy$abun # matrix of abundance of 8 species in 10 communities

#traits
#1. diets in form of factors
#OR
#2. each food type as separate trait - binary

#sites
#A. abundance 
#OR
#B. presence/absence matrix

setwd("C:/Users/Catherine/Dropbox/Phil & Catherine - Secondary forests/Analysis - Dec 2013")
# species - contains duplicate rows of species and traits in binary form
species<-read.csv("Feeding guilds 25.02.14.csv",header=T,check.names=F)
names(species)
# traits2 - contains only unique rows (1785 species) and traits in binary form
traits2<-subset(species, !duplicated(TaxonID))
names(traits2)
rownames(traits2)<-traits2[,1] #set TaxonID as row names
traits2[,1] <- NULL #remove Taxon ID column

#traits2=x with 7 traits (=food types) with binary code (form 2 - see above)

#sites2 - presence/absence data
sites2<-read.csv("Site_PresAb.csv",header=T,check.names=F)
rownames(sites2)<-sites2[,1] #set SiteID as row names
sites2[,1]<-NULL #remove SiteID column

#sites2=a with presence abundance data (form B - see above)
# for pres/ab, w.abun=FALSE
# asym.bin?
#stand.FRic - standardise by global richness - constrain between 0 and 1

#FRic = functional richness
#FEve = functional eveness
#FDiv = functional diversity - set calc.FDiv=TRUE.

#2B.1
dbFD(x=traits2,a=sites2,w.abun=FALSE,calc.FRic=TRUE,calc.FDiv=TRUE)

# gives error - add asym.bin (unsure how to apply)
#2B.2
dbFD(x=traits2,a=sites2,w.abun=FALSE,asym.bin=c(1,2,3,4,5,6,7),calc.FRic=TRUE,calc.FDiv=TRUE)
# OR are binary traits asymmetric? If not, add asym.bin as null
dbFD(x=traits2,a=sites2,w.abun=FALSE,asym.bin=NULL,calc.FRic=TRUE,calc.FDiv=TRUE)

#?cannot calculate FDiv or FEve as do not have abundance?
dbFD(x=traits2,a=sites2,w.abun=FALSE,calc.FRic=TRUE,calc.FDiv=FALSE)
