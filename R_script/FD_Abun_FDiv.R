rm(list=ls())

# load data
setwd("C:/Users/Catherine/Dropbox/Phil & Catherine - Secondary forests/Analysis - Dec 2013/Data")
data<-read.csv("FD_Abun Data To Use.csv",header=T)
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

AbFdiv1<-lme(Prop.Fdiv~SF.Age+Dist2+PFCont,random=(~1|StudyID),method="ML",data=data2)
plot(AbFdiv1) # ok
summary(AbFdiv1) # -ve AIC, +Ve logLik

#FDiv - remove SF.Age

AbFdiv2<-lme(Prop.Fdiv~Dist2+PFCont,random=(~1|StudyID),method="ML",data=data2)
plot(AbFdiv2) # ok
summary(AbFdiv2) # -ve AIC, +Ve logLik

#FDiv - remove Dist2

AbFdiv3<-lme(Prop.Fdiv~PFCont,random=(~1|StudyID),method="ML",data=data2)
plot(AbFdiv3) # ok
summary(AbFdiv3) # -ve AIC, +Ve logLik

#FDiv - null model

AbFdiv4<-lme(Prop.Fdiv~1,random=(~1|StudyID),method="ML",data=data2)
plot(AbFdiv4) # ok
summary(AbFdiv4) # -ve AIC, +Ve logLik

# try full model with all data

AbFdiv5<-lme(Prop.Fdiv~SF.Age+Dist2,random=(~1|StudyID),method="ML",data=data)
plot(AbFdiv5) # ok
summary(AbFdiv5) # -ve AIC, +Ve logLik

#FDiv - remove Age (all data)

AbFdiv6<-lme(Prop.Fdiv~Dist2,random=(~1|StudyID),method="ML",data=data)
plot(AbFdiv6) # ok
summary(AbFdiv6) # -ve AIC, +Ve logLik

#FDiv - null model - all data

AbFdiv7<-lme(Prop.Fdiv~1,random=(~1|StudyID),method="ML",data=data)
plot(AbFdiv7) # ok
summary(AbFdiv7) # -ve AIC, +Ve logLik

# FDiv - try log(Age)

AbFdiv8<-lme(Prop.Fdiv~log(SF.Age)+Dist2,random=(~1|StudyID),method="ML",data=data)
plot(AbFdiv8) # ok
summary(AbFdiv8) # -ve AIC, +Ve logLik

# FDiv - remove Dist2

AbFdiv9<-lme(Prop.Fdiv~log(SF.Age),random=(~1|StudyID),method="ML",data=data)
plot(AbFdiv9) # ok
summary(AbFdiv9) # -ve AIC, +Ve logLik

---
  # BEST MODEL
  
  #FDiv - null model - all data
  
AbFdiv7<-lme(Prop.Fdiv~1,random=(~1|StudyID),method="ML",data=dataAb)
plot(AbFdiv7) # ok
summary(AbFdiv7) # -ve AIC, +Ve logLik

# Intercept=0.9972861, SE=0.01028214