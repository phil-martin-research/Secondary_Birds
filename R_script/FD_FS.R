rm(list=ls())
sf_data<-read.csv("Data/SF Data Forest 2 No Species Groups.csv",header=T)
names(sf_data)

#forest specialists
sf_data$Sloss<-(sf_data$SFS-sf_data$PFS)/sf_data$PFS
sf_data$Sloss2<-qlogis((sf_data$Sloss+1)/2)
sf_data2<-subset(sf_data,sf_data$SFS!=0)

names(sf_data2)
plot(sf_data2$Age,(sf_data2$SFS/sf_data2$PFS))

library(nlme)

m5FS<-lme(Sloss2~log(Age)+Dist2,random=(~1|StudyID), method="ML",data=sf_data2)
plot(m5FS)
summary(m5FS)

m6FS<-lme(Sloss2~log(Age),random=(~1|StudyID), method="ML",data=sf_data2)
plot(m6FS)
summary(m6FS)
# lowest AIC and Log(Age) significant

m7FS<-lme(Sloss2~1,random=(~1|StudyID), method="ML",data=sf_data2)
plot(m7FS)
summary(m7FS)

# plot m6FS - NEED TO TRANSFORM BACK
plot(sf_data$Age,(sf_data$SFS/sf_data$PFS),xlab="Secondary forest age (years)",ylab="Proportion of forest specialists relative to primary forest",ylim=c(0,1.2))
df<-data.frame(Age=seq(1,100,0.1))
df$Preds<-plogis(predict(m6FS,newdata=df,level=0))*2
preds<--1.1424143+((log(Ages))*0.2322352)
values<-plogis(preds)*2
lines(Ages,values)
abline(h=1,lty=2)

P1<-ggplot(sf_data,aes(Age,SFS/PFS))+geom_point(size=3,shape=1)
P2<-P1+geom_line(data=df,aes(y=Preds),size=1.5)+geom_hline(lty=2,y=1)+ylab("Proportion of forest specialists \nrelative to primary forest")+xlab("Secondary forest age (years)")
P2+theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_rect(size=1.5,colour="black",fill=NA))

ggsave(filename = "Figures/Forest_specialists.pdf",width = 6,height=4,units='in',dpi=400)
ggsave(filename = "Figures/Forest_specialists.png",width = 6,height=4,units='in',dpi=400)
# new R2

library(MuMIn)

r.squaredGLMM(m6FS)

# r2 = 0.1824868

# old deviance

null.model<-lme(Sloss2~1,random=(~1|StudyID),method="ML",data=sf_data2)
null.dev<--2*logLik(null.model)[1]
best.dev<--2*logLik(m6FS)[1]
null.dev
best.dev
1-(best.dev/null.dev)

# 0.1373726 of deviance not explained by null model