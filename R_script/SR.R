DataSR<-read.csv("SF Data SR No Species Groups.csv",header=T)
names(DataSR) # different than DataSR in original analyses - groups removed

library(nlme)

plot(DataSR$Age,DataSR$Prop.SR)
mSR1<-lme(Prop.SR~Age+log(Age)+Dist2,random=(~1|StudyID),method="ML",data=DataSR)
plot(mSR1) # plot ok
summary(mSR1)

mSR2<-lme(Prop.SR~Age+log(Age),random=(~1|StudyID),method="ML",data=DataSR)
plot(mSR2) # plot ok
summary(mSR2)

mSR3<-lme(Prop.SR~Age,random=(~1|StudyID),method="ML",data=DataSR)
plot(mSR3) # plot ok
summary(mSR3) # lowest AIC but Age non-significant

mSR4<-lme(Prop.SR~1,random=(~1|StudyID),method="ML",data=DataSR)
plot(mSR4)
summary(mSR4)

# plot mSR3
max(DataSR$Prop.SR)
par(mar=c(5,6,2,2))
plot(DataSR$Age,DataSR$Prop.SR,ylim=c(0,1.4),xlab="Secondary forest age (years)",ylab="Species richness relative to primary forest")
Ages<-seq(1,100,0.1)
preds<-0.8673732+(Ages*0.0027354)
lines(Ages,preds)
abline(h=1,lty=2)

# R2 MuMIn

library(MuMIn)

r.squaredGLMM(mSR3)

# R2= 0.0453442

# OLD deviance
null.model<-lme(Prop.SR~1,random=(~1|StudyID),method="ML",data=DataSR)
null.dev<--2*logLik(null.model)[1]
best.dev<--2*logLik(mSR3)[1]
null.dev
best.dev
1-(best.dev/null.dev) # deviance=0.5199802
