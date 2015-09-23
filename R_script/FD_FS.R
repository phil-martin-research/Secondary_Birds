rm(list=ls())
data<-read.csv("SF Data Forest 2 No Species Groups.csv",header=T)
names(data)

#forest specialists
data$Sloss<-(data$SFS-data$PFS)/data$PFS
data$Sloss2<-qlogis((data$Sloss+1)/2)
data2<-subset(data,data$SFS!=0)

names(data2)
plot(data2$Age,(data2$SFS/data2$PFS))

library(nlme)

m5FS<-lme(Sloss2~log(Age)+Dist2,random=(~1|StudyID), method="ML",data=data2)
plot(m5FS)
summary(m5FS)

m6FS<-lme(Sloss2~log(Age),random=(~1|StudyID), method="ML",data=data2)
plot(m6FS)
summary(m6FS)
# lowest AIC and Log(Age) significant

m7FS<-lme(Sloss2~1,random=(~1|StudyID), method="ML",data=data2)
plot(m7FS)
summary(m7FS)

# plot m6FS - NEED TO TRANSFORM BACK
plot(data$Age,(data$SFS/data$PFS),xlab="Secondary forest age (years)",ylab="Proportion of forest specialists relative to primary forest",ylim=c(0,1.2))
Ages<-seq(1,100,0.1)
preds<--1.1424143+((log(Ages))*0.2322352)
values<-plogis(preds)*2
lines(Ages,values)
abline(h=1,lty=2)

# new R2

library(MuMIn)

r.squaredGLMM(m6FS)

# r2 = 0.1824868

# old deviance

null.model<-lme(Sloss2~1,random=(~1|StudyID),method="ML",data=data2)
null.dev<--2*logLik(null.model)[1]
best.dev<--2*logLik(m6FS)[1]
null.dev
best.dev
1-(best.dev/null.dev)

# 0.1373726 of deviance not explained by null model