#script to look at models for standardised FD effect size

#load R packages
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(lme4)
library(MuMIn)
library(afex)
library(cowplot)
library(influence.ME)

#load data on standardised FD and site data
SES_FD<-read.csv("Data/SES_summary.csv")
Sites<-read.csv("Data/Site_data.csv")
#merge data
SES_FD_sites<-merge(SES_FD,Sites,by="SiteID")
FD_comp_sites<-subset(SES_FD_sites,StudyID!=15&StudyID!=28&Disturbance_type!="Plantation"&Disturbance_type!="Fire")

#test different models against each other
M0<-lmer(SES_Diff~1+(1|study),data=FD_comp_sites)
M1<-lmer(SES_Diff~log(Age)+(1|study),data=FD_comp_sites)
M2<-lmer(SES_Diff~Disturbance_type+(1|study),data=FD_comp_sites)
AICc(M0,M1,M2)

M1_FDS<-lmer(SES_Diff~log(Age)+(1|study),data=FD_comp_sites)
FDS_coefs <- data.frame(coef(summary(M1_FDS)))
mixed(SES_Diff~log(Age)+(1|study),data=FD_comp_sites,test.intercept =T)
FDS_coefs$p.z<-c(0.0001,0.0001)
FDS_coefs$variable<-"FD_standardised"
write.csv(SpR_coefs,"Tables/Standardised_FD_Coefs_summary.csv")

#now plot results
new.data_FDS<-data.frame(Age=seq(min(FD_comp_sites$Age),max(FD_comp_sites$Age)))
FDS_preds<-predict(M1_FDS,new.data_FDS,re.form=NA,se.fit=F)
new.data_FDS$FDS<-FDS_preds
new.data_FDS$SE<-predict(M1_FDS,new.data_FDS,re.form=NA,se.fit=T)$se.fit

#plot results
theme_set(theme_cowplot())
FDS_P1<-ggplot(FD_comp_sites,aes(x=Age,y=SES_Diff))+geom_point(shape=1)+scale_x_log10()
FDS_P2<-FDS_P1+geom_line(data=new.data_FDS,aes(x=Age,y=FDS),size=1)+xlab("Time since last disturbance (Years)")
FDS_P3<-FDS_P2+geom_hline(yintercept=0,lty=2)+ylab("Relative secondary forest \nstandardised functional diversity")
FDS_P4<-FDS_P3+geom_ribbon(data=new.data_FDS,aes(x=Age,y=FDS,ymax=FDS+(2*SE),ymin=FDS-(2*SE)),alpha=0.2)
FDS_P4
ggsave("Figures/SES_FD_age.pdf",width=6,height=6,dpi=400,units="in")
ggsave("Figures/SES_FD_age.png",width=6,height=6,dpi=400,units="in")
