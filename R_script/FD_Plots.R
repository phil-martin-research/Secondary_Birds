rm(list=ls())

# load data
setwd("C:/Users/Catherine/Dropbox/Phil & Catherine - Secondary forests/Analysis - Dec 2013/Data")
sf_data<-read.csv("Data/FD Data To Use.csv",header=T)
dataAb<-read.csv("Data/FD_Abun Data To Use.csv",header=T)

library(nlme)

-----
# PLOTS with data2 - some sites removed
  
# remove sites without continuity with PF data
data2<-subset(sf_data,PFCont!="")

# SSR - null model

SSR5<-lme(Prop.SSR~1,random=(~1|StudyID),method="ML",data=data2)
plot(SSR5) # ok
summary(SSR5) # -ve AIC, +Ve logLik

# FRic -null model

FRic5<-lme(Prop.Fric~1,random=(~1|StudyID),method="ML",data=data2)
plot(FRic5) # ok
summary(FRic5)

# plots

library(ggplot2)
#create dataframe
FD_df<-data.frame(
  metric=c("SSR","FRic"),
  mean=c(0.971765,1.085704),
  se=c(0.04029374,0.1583511))
#define top and bottom of error bars
limits <- aes(ymax = mean+(1.96*se), ymin=mean-(1.96*se))
p <- ggplot(FD_df, aes(y=mean,x=metric))
p + geom_pointrange(limits) + labs(x = "Metric", y="Mean value relative to primary forest") + 
  ylim(0.6,1.6) + geom_hline(linetype="dashed",aes(yintercept=1)) +theme_bw()

----
# PLOTS with all data
  
  
# SSR - null model (all data)

SSR6<-lme(Prop.SSR~1,random=(~1|StudyID),method="ML",data=data)
plot(SSR6)
summary(SSR6)

# FRic -null model with all data

FRic6<-lme(Prop.Fric~1,random=(~1|StudyID),method="ML",data=data)
plot(FRic6)
summary(FRic6)

# FDiv - null model with all ABUNDANCE data

AbFdiv7<-lme(Prop.Fdiv~1,random=(~1|StudyID),method="ML",data=dataAb)
plot(AbFdiv7) # ok
summary(AbFdiv7) # -ve AIC, +Ve logLik

# plots

library(ggplot2)
#create dataframe
FD_df<-data.frame(
  metric=c("Functional diversity","Functional richness","Single species richness"),
  mean=c(0.9972861,1.094123,0.984084),
  se=c(0.01028214,0.1362707,0.03589116))
#define top and bottom of error bars
limits <- aes(ymax = mean+(1.96*se), ymin=mean-(1.96*se))

theme_set(theme_bw(base_size=12))
p <- ggplot(FD_df, aes(y=mean,x=metric))
p2<-p + geom_pointrange(limits,shape=1) + labs(x = "Metric", y="Value in secondary forest relative\n to primary forest") + 
  ylim(0.6,1.6) + geom_hline(linetype="dashed",aes(yintercept=1))
p2+theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size=1.5,colour="black",fill=NA))
ggsave(filename = "Figures/Functional_div.pdf",width = 6,height=4,units='in',dpi=400)
ggsave(filename = "Figures/Functional_div.png",width = 6,height=4,units='in',dpi=400)
