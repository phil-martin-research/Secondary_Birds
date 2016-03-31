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
FD_comp<-read.csv("Data/FD_summary_comp.csv")

head(FD_comp)

ggplot(data=FD_comp,aes(x=Age,y=SpR_comp))+geom_point()
ggplot(data=FD_comp,aes(x=Age,y=FDpg_comp))+geom_point()
ggplot(data=FD_comp,aes(x=Age,y=FE_comp))+geom_point()
ggplot(data=FD_comp,aes(x=Age,y=FDiv_comp))+geom_point()
ggplot(data=FD_comp,aes(x=Age,y=FDis_comp))+geom_point()
ggplot(data=FD_comp,aes(x=Age,y=FR_comp))+geom_point()

head(FD_comp)

#Species richness
M0<-lmer(SpR_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M1<-lmer(SpR_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(SpR_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)

#Functional diversity - Petchy & Gaston method
M0<-lmer(FDpg_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M1<-lmer(FDpg_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(FDpg_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)

#Functional diversity - Laliberte method
M0<-lmer(FDiv_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M1<-lmer(FDiv_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(FDiv_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)


#Functional dispersion
M0<-lmer(FDis_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M1<-lmer(FDis_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(FDis_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)

#Functional evenesss
M0<-lmer(FE_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M1<-lmer(FE_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(FE_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)

#Rao's Q
M0<-lmer(Rao_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M1<-lmer(Rao_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(Rao_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)

head(FD_comp)
ggplot(FD_comp,aes(x=SpR,y=RaoQ))+geom_point()

#summarise results
Results<-data.frame(Variable=c("Species richness","Functional diversity - P&G","Functional diversity - Laliberte","Functional dispersion","Functional eveness","Rao's Q"),
           Mean=c(-0.14367,-0.07152,-0.005805,-0.03854,-0.02691,0.13977),
           SE=c(0.06171,0.05778,0.004896,0.03497,0.01779,0.09889))

Results$UCI<-Results$Mean+(Results$SE*1.96)
Results$LCI<-Results$Mean-(Results$SE*1.96)

#plot results
theme_set(theme_bw(base_size=12))
P1<-ggplot(Results,aes(x=Variable,y=Mean,ymax=UCI,ymin=LCI))+geom_point(shape=1,size=3)+ylab("Difference between secondary \nand primary forest sites (log response ratio)")+geom_hline(yintercept=0,lty=2)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
P2<-P1+theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_rect(size=1.5,colour="black",fill=NA))
P2+scale_x_discrete(limits=c("Species richness","Functional diversity - P&G","Functional diversity - Laliberte","Functional dispersion","Functional eveness","Rao's Q"))+geom_errorbar(width=0.3)
ggsave("Figures/Pa_diversity.pdf",width = 6,height =6,units = "in",dpi = 400)
