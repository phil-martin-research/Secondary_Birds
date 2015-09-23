rm(list=ls())

# load data
setwd("C:/Users/Catherine/Dropbox/Phil & Catherine - Secondary forests/Analysis - Dec 2013/Data")
data<-read.csv("FD Data To Use.csv",header=T)
names(data)

# plot data
plot(data$Age,data$Prop.SR) # SF age vs. proportion primary forest species richness
plot(data$Age,data$Prop.SSR) # SF age vs. proportion primary forest singular species richness
plot(data$Age,data$Prop.Fric) # SF age vs. proportion primary forest functional richness
plot(data$Age,data$Prop.Feve) # SF age vs. proportion primary forest functional eveness
plot(data$Age,data$Prop.Fdiv) # SF age vs. proportion primary forest functional diversity

library(nlme)

# remove sites without continuity with PF data
data2<-subset(data,data$PFCont!="")

#FEve model

Feve1<-lme(Prop.Feve~Age+Dist2+PFCont,random=(~1|StudyID),method="ML",data=data2)
plot(Feve1) # good
summary(Feve1)

# FEve - remove Dist2

Feve2<-lme(Prop.Feve~Age+PFCont,random=(~1|StudyID),method="ML",data=data2)
plot(Feve2) # good
summary(Feve2)

#FEve - remove PFCont

Feve3<-lme(Prop.Feve~Age,random=(~1|StudyID),method="ML",data=data2)
plot(Feve3) # good
summary(Feve3)

# FEve - null model

Feve4<-lme(Prop.Feve~1,random=(~1|StudyID),method="ML",data=data2)
plot(Feve4) # ok
summary(Feve4)

#Feve - try log(Age)

Feve5<-lme(Prop.Feve~log(Age)+Dist2+PFCont,random=(~1|StudyID),method="ML",data=data2)
plot(Feve5) # good
summary(Feve5)

# not significant

# null model (Feve4) is best - Intercept= 1.176092, SE= 0.05328591

# no difference between FEve in PF and SF

# null model with all data (no exclusion due to PFCont)

Feve6<-lme(Prop.Feve~1,random=(~1|StudyID),method="ML",data=data)
plot(Feve6) # ok
summary(Feve6)

# null model (Feve6) with all data - Intercept = 1.146455, SE= 0.04862022