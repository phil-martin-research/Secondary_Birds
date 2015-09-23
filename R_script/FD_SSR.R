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

# Singular species richness model

SSR1<-lme(Prop.SSR~Age+log(Age)+Dist2+PFCont,random=(~1|StudyID),method="ML",data=data2)
plot(SSR1) # good
summary(SSR1) # -ve AIC, +Ve logLik

# SSR - remove Dist2

SSR2<-lme(Prop.SSR~Age+log(Age)+PFCont,random=(~1|StudyID),method="ML",data=data2)
plot(SSR2) # good
summary(SSR2) # -ve AIC, +Ve logLik

# SSR - remove PFCont

SSR3<-lme(Prop.SSR~Age+log(Age),random=(~1|StudyID),method="ML",data=data2)
plot(SSR3) # good
summary(SSR3) # -ve AIC, +Ve logLik

# SSR - remove Age

SSR4<-lme(Prop.SSR~log(Age),random=(~1|StudyID),method="ML",data=data2)
plot(SSR4) # ok
summary(SSR4) # -ve AIC, +Ve logLik

# SSR - null model

SSR5<-lme(Prop.SSR~1,random=(~1|StudyID),method="ML",data=data2)
plot(SSR5) # ok
summary(SSR5) # -ve AIC, +Ve logLik

# Intercept=0.971765, SE=0.04029374, null model is best

plot(data2$Age,data2$Prop.SSR, xlab="Secondary Forest Age", ylab="Singular species richness relative to primary forest")
abline(h=1,lty=2)

# mean

abline(h=0.971765,lty=2, col="red")

# confidence intervals

abline(h=0.971765+(1.96*0.04029374),lty=2, col="green")
abline(h=0.971765-(1.96*0.04029374),lty=2, col="green")

# no difference between SSR in PF and SR

# SSR - null model (all data)

SSR6<-lme(Prop.SSR~1,random=(~1|StudyID),method="ML",data=data)
plot(SSR6) # ok
summary(SSR6) # -ve AIC, +Ve logLik

# Intercept=0.984084, SE=0.03589116, null model is best