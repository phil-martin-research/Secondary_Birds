#script to calculate species richness of forest dependant speceis library(FD)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape)
library(lme4)
library(MuMIn)

#load in data
Forest_dep<-read.csv("Data/Site_Abun5.csv")
head(Forest_dep)
Forest_dep2<-subset(Forest_dep,For_dep=="High"|For_dep=="Medium")
Forest_dep2<-subset(Forest_dep,For_dep=="High")

For_rich_summary<-ddply(Forest_dep2,.(SiteID,Study,Age,PF_SF,Point_obs,Mist_nets,Transect,Vocal),summarise,For_rich=length(Vocal))

#create loop to compare richness in primary and secondary forest
SF_summary<-NULL
Un_study<-unique(Forest_dep2$Study)
for (i in 1:length(Un_study)){
  SF_sub<-subset(For_rich_summary,Study==Un_study[i]&PF_SF=="SF")
  PF_sub<-subset(For_rich_summary,Study==Un_study[i]&PF_SF=="PF")
  SF_sub$Prop_rich<-log(SF_sub$For_rich)-log(PF_sub$For_rich)
  SF_summary<-rbind(SF_sub,SF_summary)
}

ggplot(SF_summary,aes(x=Age,y=Prop_rich))+geom_point()+scale_x_log10()+geom_smooth(method="lm")

#produce models of this relationship
M0<-lmer(Prop_rich~1+(1|Study),data=SF_summary)
M1<-lmer(Prop_rich~log(Age)+(1|Study),data=SF_summary)

AICc(M0,M1)
summary(M1)
