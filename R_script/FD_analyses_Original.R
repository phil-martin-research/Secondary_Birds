##################################################################################
#script to analyse changes in functional diversity in secondary tropical forests##
##################################################################################


#I have modified the script a bit while playing around with it. I have tracked the changes and they are here
#https://github.com/PhilAMartin/Secondary_Birds

rm(list=ls())
library(FD)
library(RCurl)
library(plyr)
library(reshape)


dummy$trait #dataframe of 8 functional traits for 8 species
dummy$abun # matrix of abundance of 8 species in 10 communities

ex<-dbFD(dummy$trait,dummy$abun)
?dbFD


#traits are defined as each food type as separate buinary trait

#sites are either (1) abundance or (2) presence/absence matrix

# species - contains duplicate rows of species and traits in binary form
x<-getURL("https://raw.githubusercontent.com/PhilAMartin/Secondary_Birds/master/Data/Feeding%20guilds%2011.12.13%20v1.csv",ssl.verifypeer = FALSE)
species<-read.csv(textConnection(x),header=T)
head(species)

#create species by site matrix, species as columns, sites as row names
Sp_melt<-melt(species,id.vars =c("Scientific.name","SiteID"))
head(Sp_melt)
?cast
cast(Sp_melt,value~Scientific.name,)
?cast

# traits2 - contains only unique rows (1785 species) and traits in binary form
traits2<-subset(species, !duplicated(Taxon.ID))
names(traits2)
rownames(traits2)<-gsub(" ", ".", traits2$Scientific.name)#set TaxonID as row names

#select just columns that have trait data in
traits3<-traits2[,c(4:7)]
traits3<-traits3[order(row.names(traits3)),]


#create presence absence matrix with site as a row and species as a column

#sites2 - presence/absence data
x2<-getURL("https://raw.githubusercontent.com/PhilAMartin/Secondary_Birds/master/Data/Site_PresAb.csv",ssl.verifypeer = FALSE)
sites2<-read.csv(textConnection(x2),header=T)
head(sites2)
rownames(sites2)<-sites2[,1] #set SiteID as row names
sites2[,1]<-NULL #remove SiteID column
#reorder these so that the columns are in ascending alphabetical order 
sites2<-sites2[,order(names(sites2))]
#and do the same for traits
sites2<-sites2[,order(names(sites2))]


#sites2=a with presence abundance data (form B - see above)
# for pres/ab, w.abun=FALSE
# asym.bin?
#stand.FRic - standardise by global richness - constrain between 0 and 1

#FRic = functional richness
#FEve = functional eveness
#FDiv = functional diversity - set calc.FDiv=TRUE.

dbFD(x=traits3,a=sites2,w.abun=FALSE,calc.FRic=TRUE,calc.FDiv=F,asym.bin=NULL)

# gives error - add asym.bin (unsure how to apply)
#2B.2
dbFD(x=traits2,a=sites2)

,w.abun=FALSE,asym.bin=c(1,2,3,4,5,6,7),calc.FRic=TRUE,calc.FDiv=TRUE)
# OR are binary traits asymmetric? If not, add asym.bin as null
dbFD(x=traits2,a=sites2,w.abun=FALSE,asym.bin=NULL,calc.FRic=TRUE,calc.FDiv=TRUE)

#?cannot calculate FDiv or FEve as do not have abundance?
dbFD(x=traits2,a=sites2,w.abun=FALSE,calc.FRic=TRUE,calc.FDiv=FALSE)
