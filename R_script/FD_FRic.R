rm(list=ls())

# load data
setwd("C:/Users/Catherine/Dropbox/Phil & Catherine - Secondary forests/Analysis - Dec 2013/Data")
setwd("C:/Users/Phil/Dropbox/Work/Active projects/Catherine Sayer/Phil & Catherine - Secondary forests (1)/Analysis - Dec 2013/Data")

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
data2$PFCont<-factor(data2$PFCont)


# FRic model

FRic1<-lme(Prop.Fric~Age+log(Age)+Dist2+PFCont,random=(~1|StudyID),method="ML",data=data2)
plot(FRic1) # ok
summary(FRic1)

# FRic - remove Age

FRic2<-lme(Prop.Fric~log(Age)+Dist2+PFCont,random=(~1|StudyID),method="ML",data=data2)
plot(FRic2) # ok
summary(FRic2)

# FRic - remove log(Age)

FRic3<-lme(Prop.Fric~Dist2+PFCont,random=(~1|StudyID),method="ML",data=data2)
plot(FRic3) # ok
summary(FRic3)

# FRic - remove Dist2

FRic4<-lme(Prop.Fric~PFCont,random=(~1|StudyID),method="ML",data=data2)
plot(FRic4) # ok
summary(FRic4)

# FRic -null model

FRic5<-lme(Prop.Fric~1,random=(~1|StudyID),method="ML",data=data2)
plot(FRic5) # ok
summary(FRic5)

# intercept=1.085704, SE=0.1583511, null model is best

plot(data2$Age,data2$Prop.Fric, xlab="Secondary Forest Age", ylab="Functional richness relative to primary forest")
abline(h=1,lty=2)

# mean

abline(h=1.085704,lty=2, col="red")

# confidence intervals

abline(h=1.085704+(1.96*0.1583511),lty=2, col="green")
abline(h=1.085704-(1.96*0.1583511),lty=2, col="green")

# no difference between FRic in PF and SR

# FRic -null model with all data

FRic6<-lme(Prop.Fric~1,random=(~1|StudyID),method="ML",data=data)
plot(FRic6) # ok
summary(FRic6)

# intercept=1.094123, SE=0.1362707, null model is best
