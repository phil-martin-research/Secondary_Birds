#this script is for analysis of differences in species richness and functional diversity metrics
#between primary and secondary forest

library(FD)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(GGally)
library(reshape)
library(lme4)
library(MuMIn)

#load in data
rm(list = ls())
FD_comp<-read.csv("Data/FD_comp.csv")

head(FD_comp)

ggplot(data=FD_comp,aes(x=Age,y=SpR_comp))+geom_point()
ggplot(data=FD_comp,aes(x=Age,y=FE_comp))+geom_point()
ggplot(data=FD_comp,aes(x=Age,y=FDiv_comp))+geom_point()
ggplot(data=FD_comp,aes(x=Age,y=FDis_comp))+geom_point()
ggplot(data=FD_comp,aes(x=Age,y=FR_comp))+geom_point()

head(FD_comp)

#Species richness
M0<-lmer(SpR_comp~1+(1|No_Methods)+(1|StudyID),data=FD_comp)
M1<-lmer(SpR_comp~Age+(1|No_Methods)+(1|StudyID),data=FD_comp)
M2<-lmer(SpR_comp~log(Age)+(1|No_Methods)+(1|StudyID),data=FD_comp)
M3<-lmer(SpR_comp~nspb_P+(1|No_Methods)+(1|StudyID),data=FD_comp)
AICc(M0,M1,M2,M3)
summary(M0)


#Functional evenness
M0<-lmer(FE_comp~1+(1|No_Methods)+(1|StudyID),data=FD_comp)
M1<-lmer(FE_comp~Age+(1|No_Methods)+(1|StudyID),data=FD_comp)
M2<-lmer(FE_comp~nspb_P+(1|No_Methods)+(1|StudyID),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)

#Functional diversity
M0<-lmer(FDiv_comp~1+(1|No_Methods)+(1|StudyID),data=FD_comp)
M1<-lmer(FDiv_comp~Age+(1|No_Methods)+(1|StudyID),data=FD_comp)
M2<-lmer(FDiv_comp~nspb_P+(1|No_Methods)+(1|StudyID),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)

#Functional dispersion
M0<-lmer(FDis_comp~1+(1|No_Methods)+(1|StudyID),data=FD_comp)
M1<-lmer(FDis_comp~Age+(1|No_Methods)+(1|StudyID),data=FD_comp)
M2<-lmer(FDis_comp~nspb_P+(1|No_Methods)+(1|StudyID),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)

#Functional richness
M0<-lmer(FR_comp~1+(1|No_Methods)+(1|StudyID),data=FD_comp)
M1<-lmer(FR_comp~Age+(1|No_Methods)+(1|StudyID),data=FD_comp)
M2<-lmer(FR_comp~nspb_P+(1|No_Methods)+(1|StudyID),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)


#summarise results
Results<-data.frame(Variable=c("Species richness","Functional evenness","Functional diversity","Functional dispersion"),
           Mean=c(-0.18026,-0.0004565,-0.003702,-0.01389),
           SE=c(0.05371,0.0007043,0.006916,0.01876))

Results$UCI<-Results$Mean+(Results$SE*1.96)
Results$LCI<-Results$Mean-(Results$SE*1.96)

#plot results
theme_set(theme_bw(base_size=12))
P1<-ggplot(Results,aes(x=Variable,y=Mean,ymax=UCI,ymin=LCI))+geom_pointrange(shape=1)+ylab("log Response Ratio")+geom_hline(yintercept=0,lty=2)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
P1+theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_rect(size=1.5,colour="black",fill=NA))
ggsave("Figures/Pa_diversity.pdf",width = 6,height = 4,units = "in",dpi = 400)
