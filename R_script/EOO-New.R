rm(list=ls())
setwd("C:/Users/Catherine/Dropbox/Phil & Catherine - Secondary forests/Analysis - Dec 2013/EOO Analysis")
library(MASS)
library(nlme)

#EOO

# original analysis - date files (ALL)

dataALL<-read.csv("SF Data SR (ALL).csv",header=T)
plot(dataALL$Age,dataALL$PEOO)
plot(dataALL$Age,log(dataALL$PEOO))
plot(dataALL$Age,1/dataALL$PEOO)

dataALLC<-subset(dataALL,dataALL$PFCont!="")
names(dataALLC)

# best choice from original analysis

mod24ALL<-lme((1/PEOO)~1,random=(~1|StudyID),data=dataALL,method="ML")
plot(mod24ALL)
summary(mod24ALL)

# full model - reciprocal transformation

mod1ALL<-lme((1/PEOO)~Age+log(Age)+PFCont+Dist2,random=(~1|StudyID),data=dataALLC,method="ML")
plot(mod1ALL)
summary(mod1ALL)

# data file - excluding all estimate data (NONE)

dataNONE<-read.csv("SF Data SR (NONE).csv",header=T)
plot(dataNONE$Age,dataNONE$PEOO)
plot(dataNONE$Age,log(dataNONE$PEOO))
plot(dataNONE$Age,1/dataNONE$PEOO)

dataNONEC<-subset(dataNONE,dataNONE$PFCont!="")
names(dataNONEC)

#with dataNONEC dataset - exclude datasets without PFCont Info

mod1NONE<-lme((1/PEOO)~Age+log(Age)+PFCont+Dist2,random=(~1|StudyID),data=dataNONEC,method="ML")
plot(mod1NONE)
summary(mod1NONE)

# remove Dist2

mod2NONE<-lme((1/PEOO)~Age+log(Age)+PFCont,random=(~1|StudyID),data=dataNONEC,method="ML")
plot(mod2NONE)
summary(mod2NONE)

# remove log(Age)

mod3NONE<-lme((1/PEOO)~Age+PFCont,random=(~1|StudyID),data=dataNONEC,method="ML")
plot(mod3NONE)
summary(mod3NONE)

# remove Age

mod4NONE<-lme((1/PEOO)~PFCont,random=(~1|StudyID),data=dataNONEC,method="ML")
plot(mod4NONE)
summary(mod4NONE)

# null model (StudyID as RE)

mod5NONE<-lme((1/PEOO)~1,random=(~1|StudyID),data=dataNONEC,method="ML")
plot(mod5NONE)
summary(mod5NONE)

# with dataNONE dataset

mod6NONE<-lme((1/PEOO)~Age+log(Age)+Dist2,random=(~1|StudyID),data=dataNONE,method="ML")
plot(mod6NONE)
summary(mod6NONE)

# remove Dist2

mod7NONE<-lme((1/PEOO)~Age+log(Age),random=(~1|StudyID),data=dataNONE,method="ML")
plot(mod7NONE)
summary(mod7NONE)

# remove log(Age)

mod8NONE<-lme((1/PEOO)~Age,random=(~1|StudyID),data=dataNONE,method="ML")
plot(mod8NONE)
summary(mod8NONE)

# null model

mod9NONE<-lme((1/PEOO)~1,random=(~1|StudyID),data=dataNONE,method="ML")
plot(mod9NONE)
summary(mod9NONE)

# CWEOO

# best choice from original analysis
mod31ALL<-lme((1/PCWEOO)~1,random=(~1|StudyID),data=dataALL,na.action=na.omit,method="ML")
plot(mod31ALL)
summary(mod31ALL)

# full model - reciprocal transformation

mod1ALLCW<-lme((1/PEOO)~Age+log(Age)+PFCont+Dist2,random=(~1|StudyID),data=dataALLC,method="ML")
plot(mod1ALLCW)
summary(mod1ALLCW)

# data file - excluding all estimate data (NONE)

# remove sites with no CWEOO values
dataNONECCW<-subset(dataNONEC,dataNONEC$PCWEOO!="")

plot(dataNONE$Age,dataNONE$PCWEOO)
plot(dataNONE$Age,log(dataNONE$PCWEOO))
plot(dataNONE$Age,1/dataNONE$PCWEOO)

mod1NONECW<-lme((1/PCWEOO)~Age+log(Age)+PFCont+Dist2,random=(~1|StudyID),data=dataNONECCW,method="ML")
plot(mod1NONECW)
summary(mod1NONECW)

# remove Dist2

mod2NONECW<-lme((1/PCWEOO)~Age+log(Age)+PFCont,random=(~1|StudyID),data=dataNONECCW,method="ML")
plot(mod2NONECW)
summary(mod2NONECW)

# remove PFCont

mod3NONECW<-lme((1/PCWEOO)~Age+log(Age),random=(~1|StudyID),data=dataNONECCW,method="ML")
plot(mod3NONECW)
summary(mod3NONECW)

# remove log(Age)

mod4NONECW<-lme((1/PCWEOO)~Age,random=(~1|StudyID),data=dataNONECCW,method="ML")
plot(mod4NONECW)
summary(mod4NONECW)

# null model

mod5NONECW<-lme((1/PCWEOO)~1,random=(~1|StudyID),data=dataNONECCW,method="ML")
plot(mod5NONECW)
summary(mod5NONECW)

# with dataNONECW dataset
dataNONECW<-subset(dataNONE,dataNONE$PCWEOO!="")

mod6NONECW<-lme((1/PCWEOO)~Age+log(Age)+Dist2,random=(~1|StudyID),data=dataNONECW,method="ML")
plot(mod6NONECW)
summary(mod6NONECW)

# remove Dist2

mod7NONECW<-lme((1/PCWEOO)~Age+log(Age),random=(~1|StudyID),data=dataNONECW,method="ML")
plot(mod7NONECW)
summary(mod7NONECW)

# remove log(Age)

mod8NONECW<-lme((1/PCWEOO)~Age,random=(~1|StudyID),data=dataNONECW,method="ML")
plot(mod8NONECW)
summary(mod8NONECW)

# null model

mod9NONECW<-lme((1/PCWEOO)~1,random=(~1|StudyID),data=dataNONECW,method="ML")
plot(mod9NONECW)
summary(mod9NONECW)

# diagrams

# CWEOO - mod7NONECW

mod7NONECW<-lme((1/PCWEOO)~Age+log(Age),random=(~1|StudyID),data=dataNONECW,method="ML")
plot(mod7NONECW)
summary(mod7NONECW)

par(mfrow=c(1,1))
par(mar=c(5,6,3,2))
plot(dataNONECW$Age,dataNONECW$PCWEOO,xlab="Age of secondary forest", ylab="Proportion of primary forest \ncommunity-weighted mean global range size",ylim=c(0,7))
range(dataNONECW$Age)
Ages<-seq(1,100,0.1)
summary(mod7NONECW)
preds<-0.5873429+(Ages*-0.0115546)+((log(Ages))*0.2193131)
lines(Ages,preds)
abline(h=1,lty=2)

# CWEOO - mod3NONECW

mod3NONECW<-lme((1/PCWEOO)~Age+log(Age),random=(~1|StudyID),data=dataNONECCW,method="ML")
plot(mod3NONECW)
summary(mod3NONECW)

plot(dataNONECCW$Age,dataNONECCW$PCWEOO,xlab="Age of secondary forest", ylab="Proportion of primary forest \ncommunity-weighted mean global range size",ylim=c(0,7))
range(dataNONECCW$Age)
Ages<-seq(1,100,0.1)
summary(mod3NONECW)
preds<-0.6077158+(Ages*-0.0110606)+((log(Ages))*0.2032250)
lines(Ages,preds)
abline(h=1,lty=2)