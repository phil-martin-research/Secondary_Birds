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

#load in data
rm(list = ls())
FD_comp<-read.csv("Data/FD_abun_summary_comp.csv")
Location<-read.csv("Data/Study_location.csv")
FD_comp<-merge(FD_comp,Location,by="SiteID")
FD_comp$Lat<-FD_comp$Lat+(rnorm(length(FD_comp$Lat),0,0.00001)) 

#create a variable to quantify the number of methods used
FD_comp_new<-NULL
for (i in 1:nrow(FD_comp)){
  FD_comp_sub<-FD_comp[i,]
  Point<-ifelse(FD_comp_sub$Point_obs=="Yes",1,0)
  Mist<-ifelse(FD_comp_sub$Mist_nets=="Yes",1,0)
  Tran<-ifelse(FD_comp_sub$Transect=="Yes",1,0)
  Vocal<-ifelse(FD_comp_sub$Vocal=="Yes",1,0)
  FD_comp_sub$Methods<-sum(Point,Mist,Tran,Vocal)
  FD_comp_new<-rbind(FD_comp_new,FD_comp_sub)
}

FD_comp<-FD_comp_new

#species level metrics

#test which arrangement of random effects is best
AICc_summary<-NULL
for (i in seq(grep("SpR", (colnames(FD_comp))),grep("Rao", colnames(FD_comp)))){
  #run null models to check which random effects structure has the best fit
  M0_1<-lmer(FD_comp[[i]]~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp,REML=F)
  M0_2<-lmer(FD_comp[[i]]~1+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp,REML=F)
  M0_3<-lmer(FD_comp[[i]]~1+(1|Mist_nets)+(1|Transect)+(1|Study),data=FD_comp,REML=F)
  M0_4<-lmer(FD_comp[[i]]~1+(1|Mist_nets)+(1|Study),data=FD_comp,REML=F)
  M0_5<-lmer(FD_comp[[i]]~1+(1|Study)+(1|Methods),data=FD_comp,REML=F)
  M0_6<-lmer(FD_comp[[i]]~1+(1|Study),data=FD_comp,REML=F)
  AICc_sum<-data.frame(Variable=names(FD_comp[i]),Ran_effects=c(5,4,3,2,2,1),AICc=AICc(M0_1,M0_2,M0_3,M0_4,M0_5,M0_6)$AICc)
  AICc_summary<-rbind(AICc_summary,AICc_sum) 
}
#without exception the model with the lowest number of random effects comes out on top

#produce a loop that runs all models and then gives model selection tables and parameter estimates as an output
#and then puts the residuals into dataframe for which a spatial correlogram is run
pdf("Figures/Abundance_residuals.pdf")
model_sel_summary<-NULL
for (i in seq(grep("SpR", (colnames(FD_comp))),grep("Rao", colnames(FD_comp)))){
  #run null models to check which random effects structure has the best fit
  M0<-lmer(FD_comp[[i]]~1+(1|Study),data=FD_comp,REML=F)
  M1<-lmer(FD_comp[[i]]~1+Age+(1|Study),data=FD_comp,REML=F)
  model_sel<-data.frame(Variable=names(FD_comp[i]),x_var=c("Null model","log(Age)"),AICc=AICc(M0,M1)$AICc)
  mod_resid<-data.frame(resids=c(resid(M0),resid(M1)),fitted=c(fitted(M0),fitted(M1)),model=rep(x = c("M0","M1"),each=43),title=names(FD_comp[i]))
  Title<-names(FD_comp[i])
  print(ggplot(mod_resid,aes(x=resids))+geom_histogram()+facet_wrap(~model)+ ggtitle(Title))
  print(Resid_plots<-ggplot(mod_resid,aes(x=fitted,y=resids))+geom_point()+facet_wrap(~model)+ ggtitle(Title))
  model_sel_summary<-rbind(model_sel_summary,model_sel)
}
dev.off()

#variables that respond to age - SpR, FDiv
#variables that don't respond to age - Even_comp, FDpg_comp, FDw_comp, FE_comp, FDis_comp, Rao_comp


ggplot(FD_comp,aes(x=log(Age),y=SpR))+geom_point()+geom_smooth(se=F,method="lm")
ggplot(FD_comp,aes(x=log(Age),y=Shan_div))+geom_point()+geom_smooth(se=F,method="lm")
ggplot(FD_comp,aes(x=log(Age),y=FRic))+geom_point()+geom_smooth(se=F,method="lm")
ggplot(FD_comp,aes(x=log(Age),y=FDiv))+geom_point()+geom_smooth(se=F,method="lm")


#first look at variables that respond to age
#Species richness
M1_SPR<-lmer(SpR~log(Age)+(1|Study),data=FD_comp)
summary(M1_SPR)
r.squaredGLMM(M1_SPR)

coefs <- data.frame(coef(summary(M1_SPR)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#FDiv
M1_FDiv<-lmer(FDiv~log(Age)+(1|Study),data=FD_comp)
summary(M1_FDiv)
r.squaredGLMM(M1_FDiv)
coefs <- data.frame(coef(summary(M1_FDiv)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#plot results
new.data<-data.frame(Age=seq(min(FD_comp$Age),max(FD_comp$Age)))
new.data$FDiv<-predict(M1_FDiv,new.data,re.form=NA)
new.data$SpR<-predict(M1_SPR,new.data,re.form=NA)
new_data_preds<-melt(new.data,id.vars="Age")
head()

head(FD_comp)
FD_comp_plots<-FD_comp[,c(3,9,16)]
FD_comp_plots_melt<-melt(FD_comp_plots,id.vars="Age")
P1<-ggplot(FD_comp_plots_melt,aes(x=Age,y=value))+geom_point(shape=1)+facet_wrap(~variable,scales = "free",ncol=1)+scale_x_log10()
P2<-P1+geom_line(data=new_data_preds,aes(x=Age,y=value),size=1)
P2+geom_hline(yintercept=0,lty=2)+ylab("Difference between secondary \nand primary forset sites (log response ratio)")
ggsave("Figures/Age_models_abun.pdf",width=4,height=4,dpi=400,units="in")

#################################################
#now look at variables that don't respond to age#
#################################################

#variables that don't respond to age - FDpg, FDw, FE, FDis, Rao
head(FD_comp)
coefs_summary<-NULL
var_list<-c("FDpg","FDw","FEve","FDis")
for (i in 1:length(var_list)){
  j<-grep(var_list[i], colnames(FD_comp))
  M1<-lmer(FD_comp[[j]]~1+(1|Study),data=FD_comp)
  coefs <- data.frame(coef(summary(M1)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs$variable<-var_list[i]
  coefs_summary<-rbind(coefs,coefs_summary)
}


Plot1<-ggplot(coefs_summary,aes(x=variable,y=Estimate,ymax=Estimate+(1.96*Std..Error),ymin=Estimate-(1.96*Std..Error)))+geom_point(shape=1,size=3)+geom_errorbar(width=0.5)
Plot2<-Plot1+geom_hline(yintercept=0,lty=2)+ylab("Difference between secondary \nand primary forset sites (log response ratio)")
Plot2
ggsave("Figures/Null_models_abun.pdf",width=6,height=4,dpi=400,units="in")

