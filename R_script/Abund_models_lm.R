#this script is for analysis of differences in species richness and functional diversity metrics
#between primary and secondary forest using metrics accounting for species abundance

library(FD)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(GGally)
library(reshape)
library(lme4)
library(MuMIn)
library(ncf)
library(nlme)
library(afex)
library(cowplot)
library(influence.ME)
library(ggcorrplot)


#todo - tidy up this script

#load data
rm(list = ls())
FD_comp<-read.csv("Data/FD_abun_summary_comp.csv")
str(FD_comp)
FD_comp<-FD_comp[,-c(11,13)]
Sites<-read.csv("Data/Site_data.csv")
FD_comp_sites<-merge(FD_comp,Sites,by="SiteID")
head(FD_comp_sites)
FD_comp_sites<-FD_comp_sites[,-c(16,18)]
FD_comp_sites<-subset(FD_comp_sites,StudyID!=15&StudyID!=28&Disturbance_type!="Plantation"&Disturbance_type!="Fire")

#before any analyses assess whether there is any correlation between the different measures of
#diversity we are using here
head(FD_comp_sites)
FD_comp_sites_sub<- FD_comp_sites[,c(9:14)]
head(FD_comp_sites_sub)
cormat <- round(cor(FD_comp_sites_sub),2)
ggcorr(FD_comp_sites_sub,label = T)
#there is strong correlation between Spr, FDpg and Fric.
#THis suggests that we should only use one of these metrics
ggsave("Figures/Correlation_matrix.png",width=6,height=6,dpi=100,units="in")


#create a variable to quantify the number of methods used
FD_comp_sites_new<-NULL
for (i in 1:nrow(FD_comp_sites)){
  FD_comp_sites_sub<-FD_comp_sites[i,]
  Point<-ifelse(FD_comp_sites_sub$Point_obs=="Yes",1,0)
  Mist<-ifelse(FD_comp_sites_sub$Mist_nets=="Yes",1,0)
  Tran<-ifelse(FD_comp_sites_sub$Transect=="Yes",1,0)
  Vocal<-ifelse(FD_comp_sites_sub$Vocal=="Yes",1,0)
  FD_comp_sites_sub$Methods<-sum(Point,Mist,Tran,Vocal)
  FD_comp_sites_new<-rbind(FD_comp_sites_new,FD_comp_sites_sub)
}
FD_comp_sites<-FD_comp_sites_new

#now test which arrangement of random effects is best
AICc_summary<-NULL
for (i in seq(grep("SpR", (colnames(FD_comp_sites))),grep("FDis", colnames(FD_comp_sites)))){
  #run null models to check which random effects structure has the best fit
  M0_1<-lmer(FD_comp_sites[[i]]~log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_sites,REML=T)
  M0_2<-lmer(FD_comp_sites[[i]]~log(Age)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_sites,REML=T)
  M0_3<-lmer(FD_comp_sites[[i]]~log(Age)+(1|Mist_nets)+(1|Transect)+(1|Study),data=FD_comp_sites,REML=T)
  M0_4<-lmer(FD_comp_sites[[i]]~log(Age)+(1|Mist_nets)+(1|Study),data=FD_comp_sites,REML=T)
  M0_5<-lmer(FD_comp_sites[[i]]~log(Age)+(1|Study)+(1|Methods),data=FD_comp_sites,REML=T)
  M0_6<-lmer(FD_comp_sites[[i]]~log(Age)+(1|Study),data=FD_comp_sites,REML=T)
  AICc_sum<-data.frame(Variable=names(FD_comp_sites[i]),
                       Random_effects=c("Point obs+Mist nets+Transect+Vocal+Study",
                                        "Mist nets+Transect+Vocal+Study",
                                        "Mist nets+ Transect+Vocal+Study",
                                        "Mist_nets+Transect+Study",
                                        "Mist_nets+Study",
                                        "Study"),
                       No_Ran_effects=c(5,4,3,2,2,1),
                       AICc=AICc(M0_1,M0_2,M0_3,M0_4,M0_5,M0_6)$AICc)
  AICc_summary<-rbind(AICc_summary,AICc_sum) 
}


#without exception the model with the lowest number of random effects comes out on top
AICc_summary<-AICc_summary[with(AICc_summary, order(Variable, AICc)), ]
AICc_summary$Rank<-c(1,2,3,4,5,6)
write.csv(AICc_summary,"Tables/Random_model_AICc.csv")

#########################################################################
#following reviewers comments look at the influence of data points 
#when age is logged or unlogged
#########################################################################


FD_comp_sites<-FD_comp_sites[complete.cases(FD_comp_sites),]
M1<-lmer(SpR~1+log(Age)+(1|Study),data=FD_comp_sites,REML=F)
M2<-lmer(SpR~1+Age+(1|Study),data=FD_comp_sites,REML=F)
infl1 <- influence(M1, obs = TRUE)
infl2 <- influence(M2, obs = TRUE)
Log_plot<-qplot(FD_comp_sites$Age,cooks.distance(infl1))+geom_smooth(se=F,method="lm")+xlab("Time since last disturbance (Years)")+ylab("Cook's distance")
Unlogged_plot<-qplot(FD_comp_sites$Age,cooks.distance(infl2))+geom_smooth(se=F,method="lm")+xlab("Time since last disturbance (Years)")+ylab("Cook's distance")
both_plots<-plot_grid(Unlogged_plot, Log_plot, labels = c("A", "B"))
save_plot("figures/unlogged_log.png", both_plots,
         ncol = 2, # we're saving a grid plot of 2 columns
         nrow = 1, # and 2 rows
         # each individual subplot should have an aspect ratio of 1.3
         base_aspect_ratio = 1.3
)


##########################################################################
#comparison of alternative models for species richness and functional ####
#diversity measures#######################################################
##########################################################################

#This loop runs all alternative models and then produces 
#a model selection tables and parameter estimates as an output

FD_comp_sites_subset<-subset(FD_comp_sites,Disturbance_type!="")
model_sel_summary<-NULL
for (i in seq(grep("SpR", (colnames(FD_comp_sites_subset))),grep("FDis", colnames(FD_comp_sites_subset)))){
  #run null models to check which random effects structure has the best fit
  M0<-lmer(FD_comp_sites_subset[[i]]~1+(1|Study),data=FD_comp_sites_subset,REML=F)
  M1<-lmer(FD_comp_sites_subset[[i]]~1+log(Age)+(1|Study),data=FD_comp_sites_subset,REML=F)
  M2<-lmer(FD_comp_sites_subset[[i]]~1+Disturbance_type+(1|Study),data=FD_comp_sites_subset,REML=F)
  model_sel<-data.frame(Variable=names(FD_comp_sites_subset[i]),
                        x_var=c("Null model","log(Age)","Disturbance type"),
                        AICc=AICc(M0,M1,M2)$AICc)
  model_sel$delta<-(model_sel$AICc)-min(model_sel$AICc)
  model_sel$R2<-round(c(r.squaredGLMM(M0)[1],r.squaredGLMM(M1)[1],r.squaredGLMM(M2)[1]),2)
  model_sel_summary<-rbind(model_sel_summary,model_sel)
}

write.csv(model_sel_summary,"Tables/model_sel_summary.csv")

#model comparison shows that models inclduing disturbance never perform better than
#those models containing age, or null models
#Species richness and Functional dispersal respond to age, but everything else fits a null model

################################################################
#produce coefficients for the most parsimoniuous models#########
################################################################

#Species richness
M1_SPR<-lmer(SpR~log(Age)+(1|Study),data=FD_comp_sites_subset)
SpR_coefs <- data.frame(coef(summary(M1_SPR)))
mixed(SpR~log(Age)+(1|Study),data=FD_comp_sites_subset,test.intercept =T)
SpR_coefs$p.z<-c(0.03,0.09)
SpR_coefs$variable<-"SpR"

#FDis
M1_FDis<-lmer(FDis~log(Age)+(1|Study),data=FD_comp_sites_subset)
FDiv_coefs <- data.frame(coef(summary(M1_FDis)))
mixed(FDis~log(Age)+(1|Study),data=FD_comp_sites_subset,test.intercept =T)
FDiv_coefs$p.z<-c(0.37,0.05)
FDiv_coefs$variable<-"FDiv"

#plot results
new.data_FDis<-data.frame(Age=seq(min(FD_comp_sites_subset$Age),max(FD_comp_sites_subset$Age)))
FDis_preds<-predict(M1_FDis,new.data_FDis,re.form=NA,se.fit=F)
new.data_FDis$FDis<-FDis_preds
new.data_FDis$SE<-predict(M1_FDis,new.data_FDis,re.form=NA,se.fit=T)$se.fit


new.data_SpR<-data.frame(Age=seq(min(FD_comp_sites_subset$Age),max(FD_comp_sites_subset$Age)))
new.data_SpR$SpR<-predict(M1_SPR,new.data_SpR,re.form=NA)
new.data_SpR$SE<-predict(M1_SPR,new.data_SpR,re.form=NA,se.fit=T)$se.fit

#plot for species richness
theme_set(theme_cowplot())
SPR_P1<-ggplot(FD_comp_sites_subset,aes(x=Age,y=SpR))+geom_point(shape=1)+scale_x_log10()
SPR_P2<-SPR_P1+geom_line(data=new.data_SpR,aes(x=Age,y=SpR),size=1)+xlab("Time since last disturbance (Years)")
SPR_P3<-SPR_P2+geom_hline(yintercept=0,lty=2)+ylab("Relative secondary forest \nspecies richness (response ratio)")
SPR_P4<-SPR_P3+geom_ribbon(data=new.data_SpR,aes(x=Age,y=SpR,ymax=SpR+(2*SE),ymin=SpR-(2*SE)),alpha=0.5)

#plot for FDis
FDis_P1<-ggplot(FD_comp_sites_subset,aes(x=Age,y=FDis))+geom_point(shape=1)+scale_x_log10()
FDis_P2<-FDis_P1+geom_line(data=new.data_FDis,aes(x=Age,y=FDis),size=1)+xlab("Time since last disturbance (Years)")
FDis_P3<-FDis_P2+geom_hline(yintercept=0,lty=2)+ylab("Relative secondary forest \nfunctional dispersion (response ratio)")
FDis_P4<-FDis_P3+geom_ribbon(data=new.data_FDis,aes(x=Age,y=FDis,ymax=FDis+(2*SE),ymin=FDis-(2*SE)),alpha=0.5)

#################################################
#now look at variables that don't respond to age#
#################################################

#variables that don't respond to age - FDpg, FE, FDiv, Fric
head(FD_comp_sites_subset)
coefs_summary<-NULL
var_list<-c("FDpg","FRic","FEve", "FDis")
for (i in 1:length(var_list)){
  j<-grep(var_list[i], colnames(FD_comp_sites_subset))
  M1<-lmer(FD_comp_sites_subset[[j]]~1+(1|Study),data=FD_comp_sites_subset)
  coefs <- data.frame(coef(summary(M1)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs$variable<-var_list[i]
  coefs_summary<-rbind(coefs,coefs_summary)
}

Coefs_summary<-rbind(coefs_summary,FDiv_coefs,SpR_coefs)

coefs_summary$var2<-c("Functional \ndispersal \n(FDiv)","Functional \nEvenness \n(FEve)","Functional \nrichness \n(FRic)","Functional \nDiversity \n(FD)")

write.csv(Coefs_summary,"Tables/Coefs_summary.csv")

Plot1<-ggplot(coefs_summary,aes(x=var2,y=Estimate,ymax=Estimate+(1.96*Std..Error),ymin=Estimate-(1.96*Std..Error)))+geom_point(shape=1,size=3)+geom_errorbar(width=0.5)
Plot2<-Plot1+geom_hline(yintercept=0,lty=2)+ylab("Difference between secondary \nand primary forest sites (log response ratio)")
Plot2+xlab("Diversity variable")
ggsave("Figures/Null_models_abun.pdf",width=6,height=6,dpi=400,units="in")
ggsave("Figures/Null_models_abun.png",width=6,height=6,dpi=400,units="in")


#############################################################
#this section is for modelling of specialist species richness#
##############################################################

#load in data
Forest_dep<-read.csv("Data/Site_Abun5.csv")
head(Forest_dep)
Forest_dep2<-subset(Forest_dep,For_dep=="High"|For_dep=="Medium")
Forest_dep2<-subset(Forest_dep,For_dep=="High")
Forest_dep2<-subset(Forest_dep2,Study!=15&Study!=28&Study!=23)

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


SF_summary2<-merge(SF_summary,Sites,by = "SiteID")
SF_summary3<-subset(SF_summary2,StudyID!=15&StudyID!=28&Disturbance_type!="Plantation"&Disturbance_type!="Fire")


#produce models of this relationship
M0<-lmer(Prop_rich~1+(1|Study),data=SF_summary3)
M1<-lmer(Prop_rich~log(Age)+(1|Study),data=SF_summary3)
M2<-lmer(Prop_rich~Disturbance_type+(1|Study),data=SF_summary3)

AICc(M0,M1,M2)
summary(M1)

r.squaredGLMM(M1)

#plot results
new.data<-data.frame(Age=seq(min(SF_summary3$Age),max(SF_summary3$Age),length.out = 1000))
new.data$Prop_rich<-predict(M1,new.data,re.form=NA,se.fit=T)$fit
new.data$SE<-predict(M1,new.data,re.form=NA,se.fit=T)$se.fit


ForSpec_P1<-ggplot(SF_summary3,aes(x=Age,y=Prop_rich))+geom_point(shape=1)+scale_x_log10()
ForSpec_P2<-ForSpec_P1+geom_line(data=new.data,aes(x=Age,y=Prop_rich),size=1)+xlab("Time since last disturbance (Years)")
ForSpec_P3<-ForSpec_P2+geom_hline(yintercept=0,lty=2)+ylab("Relative secondary forest \nspecialist richness (response ratio)")
ForSpec_P4<-ForSpec_P3+geom_ribbon(data=new.data,aes(x=Age,y=Prop_rich,ymax=Prop_rich+(2*SE),ymin=Prop_rich-(2*SE)),alpha=0.5)

time_plots<-plot_grid(SPR_P4, ForSpec_P4,FDis_P4, labels = c("(a)", "(b)","(c)"), align = "h",ncol = 1)
save_plot("Figures/time_plots.pdf", time_plots,
          nrow = 3, # we're saving a grid plot of 2 columns
          # each individual subplot should have an aspect ratio of 1.3
          base_height = 4,
          base_width=5
)
save_plot("Figures/time_plots.png", time_plots,
          nrow = 3, # we're saving a grid plot of 2 columns
          # each individual subplot should have an aspect ratio of 1.3
          base_height =4,
          base_width=5
)
