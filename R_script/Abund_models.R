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

#species level metrics

#Species richness
M0<-lmer(SpR_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M1<-lmer(SpR_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(SpR_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)

head(FD_comp)

#Shannon diversity
M0<-lmer(Shan_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M1<-lmer(Shan_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(Shan_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)
summary(M1)

#Community evenness
M0<-lmer(Even_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M1<-lmer(Even_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(Even_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
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
#functional richness is considerably lower in secondary forests - but the errors are large

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


#now look at community weighted means
head(FD_comp)
FD_comp_CWD<-FD_comp[,c(1:8,59:ncol(FD_comp))]
#remove infinite values
is.na(FD_comp_CWD) <- sapply(FD_comp_CWD, is.infinite)
FD_comp_CWD[is.na(FD_comp_CWD)] <- 0

#reorganise the data
FD_comp_CWD_melt<-melt(FD_comp_CWD,id.vars=c("SiteID","Study","Age","PF_SF","Point_obs","Mist_nets","Transect","Vocal"))
ggplot(FD_comp_CWD_melt, aes(x=variable, y=value)) + stat_summary(fun.y="mean", geom="point")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#this looks interesting, now lets test this statistically
head(FD_comp_CWD)

#Invertebrate feeding species
M0<-lmer(Diet.Inv_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
M1<-lmer(Diet.Inv_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
M2<-lmer(Diet.Inv_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
AICc(M0,M1,M2)
summary(M0)

#Vegetative (?) feeding species
M0<-lmer(Diet.Vend_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
M1<-lmer(Diet.Vend_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
M2<-lmer(Diet.Vend_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
AICc(M0,M1,M2)
summary(M0)

#Body size
M0<-lmer(BodyMass.Value_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
M1<-lmer(BodyMass.Value_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
M2<-lmer(BodyMass.Value_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
AICc(M0,M1,M2)
summary(M0)

#aerial foraging
M0<-lmer(ForStrat.aerial_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
M1<-lmer(ForStrat.aerial_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
M2<-lmer(ForStrat.aerial_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
AICc(M0,M1,M2)
summary(M0)


head(FD_comp_CWD)
