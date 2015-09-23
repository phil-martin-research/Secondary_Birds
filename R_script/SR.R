DataSR<-read.csv("Data/SF Data SR No Species Groups.csv",header=T)
names(DataSR) # different than DataSR in original analyses - groups removed

library(nlme)
library(ggplot2)

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
df<-data.frame(Age=seq(1,100,0.1))
df$Prop.SR<-predict(mSR3,newdata=df,level=0)
P1<-ggplot(DataSR,aes(x=Age,y=Prop.SR))+geom_point(size=3,shape=1)+xlab("Secondary forest age (years)")+ylab("Species richness relative to primary forest")
P2<-P1+geom_line(data=df,size=1.5)+geom_hline(y=1,lty=2)
P2+theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_rect(size=1.5,colour="black",fill=NA))
ggsave(filename = "Figures/Species_richness.pdf",width = 6,height=4,units='in',dpi=400)
ggsave(filename = "Figures/Species_richness.png",width = 6,height=4,units='in',dpi=400)



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
