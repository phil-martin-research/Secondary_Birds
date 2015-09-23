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

#FDiv model

Fdiv1<-lme(Prop.Fdiv~Age+Dist2+PFCont,random=(~1|StudyID),method="ML",data=data2)
plot(Fdiv1) # good
summary(Fdiv1) # -ve AIC, +Ve logLik

#FDiv - remove Dist2

Fdiv2<-lme(Prop.Fdiv~Age+PFCont,random=(~1|StudyID),method="ML",data=data2)
plot(Fdiv2) # ok
summary(Fdiv2) # -ve AIC, +Ve logLik

#FDiv - remove Age

Fdiv3<-lme(Prop.Fdiv~PFCont,random=(~1|StudyID),method="ML",data=data2)
plot(Fdiv3) # ok
summary(Fdiv3) # -ve AIC, +Ve logLik

# FDiv - null model (data2)

Fdiv4<-lme(Prop.Fdiv~1,random=(~1|StudyID),method="ML",data=data2)
plot(Fdiv4) # ok
summary(Fdiv4) # -ve AIC, +Ve logLik

#FDiv - try log(Age)

Fdiv5<-lme(Prop.Fdiv~log(Age)+Dist2+PFCont,random=(~1|StudyID),method="ML",data=data2)
plot(Fdiv5) # good
summary(Fdiv5) # -ve AIC, +Ve logLik

# not significant

# null model (Fdiv4) is best, Intercept=0.9977873, SE=0.008632664

# no difference between FDiv in PF and SF

# null model with all data (no exclusion due to PFCont)

Fdiv6<-lme(Prop.Fdiv~1,random=(~1|StudyID),method="ML",data=data)
plot(Fdiv6) # ok
summary(Fdiv6)

# null model (Fdiv6) with all data - Intercept=0.9974589, SE=0.007544835