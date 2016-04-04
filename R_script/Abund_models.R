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

#load in data
rm(list = ls())
FD_comp<-read.csv("Data/FD_abun_summary_comp.csv")

head(FD_comp)

ggplot(data=FD_comp,aes(x=Age,y=SpR_comp))+geom_point()
ggplot(data=FD_comp,aes(x=Age,y=FDpg_comp))+geom_point()
ggplot(data=FD_comp,aes(x=Age,y=FE_comp))+geom_point()
ggplot(data=FD_comp,aes(x=Age,y=FDiv_comp))+geom_point()
ggplot(data=FD_comp,aes(x=Age,y=FDis_comp))+geom_point()
ggplot(data=FD_comp,aes(x=Age,y=FR_comp))+geom_point()

FD_comp2<-FD_comp[,c(1,26:ncol(FD_comp))]
FD_comp_melt<-melt(FD_comp2,id.vars="SiteID")

ggplot(data=FD_comp_melt,aes(x=variable,y=value))+geom_boxplot()+geom_point()


#Species richness
M0<-lmer(SpR_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M1<-lmer(SpR_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(SpR_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)

#Functional diversity - Petchy & Gaston method
M0<-lmer(FDpg_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M1<-lmer(FDpg_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(FDpg_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)

#Functional diversity - Petchy & Gaston method weighted by abundance
M0<-lmer(FDw_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Study),data=FD_comp)
M1<-lmer(FDw_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Study),data=FD_comp)
M2<-lmer(FDw_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)

#Functional richness
M0<-lmer(FR_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M1<-lmer(FR_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(FR_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)
#functional richness is considerably lower in secondary forests

#Functional richness
M0<-lmer(FR_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M1<-lmer(FR_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(FR_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)
#functional richness is considerably lower in secondary forests

#Functional diversity - Laliberte method
M0<-lmer(FDiv_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M1<-lmer(FDiv_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(FDiv_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)


#Functional dispersion
M0<-lmer(FDis_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M1<-lmer(FDis_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(FDis_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)

#Functional evenesss
M0<-lmer(FE_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M1<-lmer(FE_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(FE_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)

#Rao's Q
M0<-lmer(Rao_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M1<-lmer(Rao_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(Rao_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)
