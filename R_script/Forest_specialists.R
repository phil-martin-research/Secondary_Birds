#script to calculate species richness of forest dependant speceis library(FD)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape)
library(lme4)
library(MuMIn)
library(cowplot)
library(afex)

#load in data
Forest_dep<-read.csv("Data/Site_Abun5.csv")
head(Forest_dep)
Forest_dep2<-subset(Forest_dep,For_dep=="High"|For_dep=="Medium")
Forest_dep2<-subset(Forest_dep,For_dep=="High")

For_rich_summary<-ddply(Forest_dep2,.(SiteID,Study,Age,PF_SF,Point_obs,Mist_nets,Transect,Vocal),summarise,For_rich=length(Vocal))

#create loop to compare richness in primary and secondary forest
SF_summary<-NULL
Un_study<-unique(Forest_dep2$Study)
for (i in 1:length(Un_study)){
  SF_sub<-subset(For_rich_summary,Study==Un_study[i]&PF_SF=="SF")
  PF_sub<-subset(For_rich_summary,Study==Un_study[i]&PF_SF=="PF")
  SF_sub$Prop_rich<-log(SF_sub$For_rich)-log(PF_sub$For_rich)
  SF_summary<-rbind(SF_sub,SF_summary)
}

#look at the different random effects structures

#create a variable to quantify the number of methods used
SF_summary_new<-NULL
for (i in 1:nrow(FD_comp)){
  FD_comp_sub<-SF_summary[i,]
  Point<-ifelse(FD_comp_sub$Point_obs=="Yes",1,0)
  Mist<-ifelse(FD_comp_sub$Mist_nets=="Yes",1,0)
  Tran<-ifelse(FD_comp_sub$Transect=="Yes",1,0)
  Vocal<-ifelse(FD_comp_sub$Vocal=="Yes",1,0)
  FD_comp_sub$Methods<-sum(Point,Mist,Tran,Vocal)
  SF_summary_new<-rbind(SF_summary_new,FD_comp_sub)
}


M0_1<-lmer(Prop_rich~log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=SF_summary_new,REML=T)
M0_2<-lmer(Prop_rich~log(Age)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=SF_summary_new,REML=T)
M0_3<-lmer(Prop_rich~log(Age)+(1|Mist_nets)+(1|Transect)+(1|Study),data=SF_summary_new,REML=T)
M0_4<-lmer(Prop_rich~log(Age)+(1|Mist_nets)+(1|Study),data=SF_summary_new,REML=T)
M0_5<-lmer(Prop_rich~log(Age)+(1|Study)+(1|Methods),data=SF_summary_new,REML=T)
M0_6<-lmer(Prop_rich~log(Age)+(1|Study),data=SF_summary_new,REML=T)
AICc(M0_1,M0_2,M0_3,M0_4,M0_5,M0_6)



#produce models of this relationship
M0<-lmer(Prop_rich~1+(1|Study),data=SF_summary)
M1<-lmer(Prop_rich~log(Age)+(1|Study),data=SF_summary)

AICc(M0,M1)
summary(M1)

mixed(Prop_rich~log(Age)+(1|Study),data=SF_summary,test.intercept =T)

r.squaredGLMM(M1)

#plot results
new.data<-data.frame(Age=seq(min(SF_summary$Age),max(SF_summary$Age),length.out = 100))
new.data$Prop_rich<-predict(M1,new.data,re.form=NA,se.fit=T)$fit
new.data$SE<-predict(M1,new.data,re.form=NA,se.fit=T)$se.fit


head(new_data_preds)

theme_set(theme_bw(base_size=12))
P1<-ggplot(SF_summary,aes(x=Age,y=Prop_rich))+geom_point(shape=1)+scale_x_log10()
P2<-P1+geom_line(data=new.data,aes(x=Age,y=Prop_rich),size=1)+xlab("Age (years - log scale)")
P3<-P2+geom_hline(yintercept=0,lty=2)+ylab("Difference between secondary \nand primary forest sites (log response ratio)")
<-P3+geom_ribbon(data=new.data,aes(x=Age,y=Prop_rich,ymax=Prop_rich+(2*SE),ymin=Prop_rich-(2*SE)),alpha=0.5)
ggsave("Figures/Age_models_specialist.pdf",width=8,height=4,dpi=400,units="in")
ggsave("Figures/Age_models_specialist.png",width=8,height=4,dpi=400,units="in")
