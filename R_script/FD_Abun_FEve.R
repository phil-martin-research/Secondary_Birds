rm(list=ls())

# load data
setwd("C:/Users/Catherine/Dropbox/Phil & Catherine - Secondary forests/Analysis - Dec 2013/Data")
data<-read.csv("FD_Abun Data To Use.csv",header=T)
names(data)

library(nlme)

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

AbFeve1<-lme(Prop.Feve~SF.Age+Dist2+PFCont,random=(~1|StudyID),method="ML",data=data2)
plot(AbFeve1) # good
summary(AbFeve1)

# FEve model - remove Dist2

AbFeve2<-lme(Prop.Feve~SF.Age+PFCont,random=(~1|StudyID),method="ML",data=data2)
plot(AbFeve2) # good
summary(AbFeve2)

#FEve model - remove PFCont

AbFeve3<-lme(Prop.Feve~SF.Age,random=(~1|StudyID),method="ML",data=data2)
plot(AbFeve3) # good
summary(AbFeve3)

# Feve - null model

AbFeve4<-lme(Prop.Feve~1,random=(~1|StudyID),method="ML",data=data2)
plot(AbFeve4) # ok
summary(AbFeve4)

# FEve - try with log(Age)

AbFeve5<-lme(Prop.Feve~log(SF.Age),random=(~1|StudyID),method="ML",data=data2)
plot(AbFeve5) # good
summary(AbFeve5)

plot(log(data$SF.Age),data$Prop.Feve) # log(SF age) vs. proportion primary forest functional eveness

# lower AIC than null model but log(Age) is not signigicant (0.0579) - relationship seems to be based on single point

# try with all data (data)

#FEve model (data)

AbFeve6<-lme(Prop.Feve~SF.Age+Dist2,random=(~1|StudyID),method="ML",data=data)
plot(AbFeve6) # good
summary(AbFeve6)

#FEve model - remove Dist2

AbFeve7<-lme(Prop.Feve~SF.Age,random=(~1|StudyID),method="ML",data=data)
plot(AbFeve7) # good
summary(AbFeve7)

# Feve model - mull model - all data

AbFeve8<-lme(Prop.Feve~1,random=(~1|StudyID),method="ML",data=data)
plot(AbFeve8) # ok
summary(AbFeve8)

# FEve - try with log(Age) - all data

AbFeve9<-lme(Prop.Feve~log(SF.Age),random=(~1|StudyID),method="ML",data=data)
plot(AbFeve9) # good
summary(AbFeve9)

# log(Age) significant and low AIC

# try full model with log(Age)

AbFeve10<-lme(Prop.Feve~log(SF.Age)+Dist2,random=(~1|StudyID),method="ML",data=data)
plot(AbFeve10) # good
summary(AbFeve10)

# not significant

---
# BEST MODEL
  
AbFeve9<-lme(Prop.Feve~log(SF.Age),random=(~1|StudyID),method="ML",data=data)
plot(AbFeve9) # good
summary(AbFeve9)

# R2 MuMIn

library(MuMIn)

r.squaredGLMM(AbFeve9)

# R2 = 0.1804323

# deviance OLD METHOD
null.model<-lme(Prop.Feve~1,random=(~1|StudyID),method="ML",data=data)
null.dev<--2*logLik(null.model)[1]
best.dev<--2*logLik(AbFeve9)[1]
null.dev
best.dev
1-(best.dev/null.dev)

plot(log(data$SF.Age),data$Prop.Feve) # log(SF age) vs. proportion primary forest functional eveness
plot(data$SF.Age,data$Prop.Feve)

#                 Value  Std.Error DF   t-value p-value
#(Intercept)  1.3539224 0.10975260 19 12.336131  0.0000
#log(SF.Age) -0.1061108 0.03916426 11 -2.709378  0.0203

plot(data$SF.Age,data$Prop.Feve,xlab="Age of Secondary Forest", ylab="Functional Evenness relative to primary forest",ylim=c(0,2))
range(data$SF.Age)
Ages<-seq(1,100,0.1)
preds<-1.3539224+((log(Ages)*-0.1061108))
lines(Ages,preds)
abline(h=1,lty=2)
predsUP<-(1.3539224+(0.10975260*1.96)+((log(Ages)*(-0.1061108+(0.03916426*1.96)))))
lines(Ages,predsUP,lty=2,col="green")
predsLOW<-(1.3539224-(0.10975260*1.96)+((log(Ages)*(-0.1061108-(0.03916426*1.96)))))
lines(Ages,predsLOW,lty=2,col="green")
