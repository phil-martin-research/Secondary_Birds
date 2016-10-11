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


#load in data
rm(list = ls())
FD_comp<-read.csv("Data/FD_abun_summary_comp.csv")
Location<-read.csv("Data/Study_location.csv")
FD_comp<-merge(FD_comp,Location,by="SiteID")
FD_comp$Lat<-FD_comp$Lat+(rnorm(length(FD_comp$Lat),0,0.00001)) 


#before any analyses assess whether there is any correlation between the different measures of
#diversity we are using here

head(FD_comp)
FD_comp_sub<- FD_comp[, c(9,12:17)]
head(FD_comp_sub)
cormat <- round(cor(FD_comp_sub),2)
ggcorr(FD_comp_sub,label = T)
ggsave("Figures/Correlation_matrix.png",width=6,height=6,dpi=100,units="in")
ggpairs(FD_comp_sub)

#bootstrap to see which model is best for species richness
Mod_sel<-NULL
Mod_sel_summary<-NULL
for (i in 1:1000){
  Ran_samp<-FD_comp[sample(nrow(FD_comp),size=nrow(FD_comp),replace=TRUE),]
  for (j in seq(grep("SpR_comp", (colnames(Ran_samp))),grep("Rao_comp", colnames(Ran_samp)))){
  M0<-lm(Ran_samp[[j]]~1,data=Ran_samp)
  M1<-lm(Ran_samp[[j]]~log(Age),data=Ran_samp)
  Mod_AICc<-AICc(M0,M1)[,2]
  Mod_sel_sub<-data.frame(AICc=Mod_AICc,R2=c(0,summary(M1)$r.squared),model=c("M0","M1"))
  Mod_sel_sub$Rank1<-ifelse(Mod_sel_sub$AICc==min(Mod_sel_sub$AICc),1,0)
  Mod_sel_sub$variable<-names(Ran_samp[j])
  Mod_sel<-rbind(Mod_sel_sub,Mod_sel)
  }
  Mod_sel_summary<-rbind(Mod_sel_summary,Mod_sel)
  print(i)
}

ddply(Mod_sel_summary,.(variable,model),summarise,M_AICc=mean(AICc),rank=sum(Rank1)/500000)


lm(SpR_comp~log(Age),data=FD_comp)

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
  M0_1<-lmer(FD_comp[[i]]~log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp,REML=T)
  M0_2<-lmer(FD_comp[[i]]~log(Age)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp,REML=T)
  M0_3<-lmer(FD_comp[[i]]~log(Age)+(1|Mist_nets)+(1|Transect)+(1|Study),data=FD_comp,REML=T)
  M0_4<-lmer(FD_comp[[i]]~log(Age)+(1|Mist_nets)+(1|Study),data=FD_comp,REML=T)
  M0_5<-lmer(FD_comp[[i]]~log(Age)+(1|Study)+(1|Methods),data=FD_comp,REML=T)
  M0_6<-lmer(FD_comp[[i]]~log(Age)+(1|Study),data=FD_comp,REML=T)
  AICc_sum<-data.frame(Variable=names(FD_comp[i]),
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

#now look at the influence of data points when age is logged or unlogged
#this is for the reviewers comments
M1<-lmer(SpR~1+log(Age)+(1|Study),data=FD_comp,REML=F)
M2<-lmer(SpR~1+Age+(1|Study),data=FD_comp,REML=F)

infl1 <- influence(M1, obs = TRUE)
infl2 <- influence(M2, obs = TRUE)

plot(infl1, which = "cook")
plot(infl2, which = "cook")

Log_plot<-qplot(FD_comp$Age,cooks.distance(infl1))+geom_smooth(se=F,method="lm")+xlab("Time since last disturbance (Years)")+ylab("Cook's distance")+coord_cartesian(xlim=c(0,105),ylim=c(-0.01,0.25),expand = F)
Unlogged_plot<-qplot(FD_comp$Age,cooks.distance(infl2))+geom_smooth(se=F,method="lm")+xlab("Time since last disturbance (Years)")+ylab("Cook's distance")+coord_cartesian(xlim=c(0,105),ylim=c(-0.01,0.25),expand = F)
both_plots<-plot_grid(Unlogged_plot, Log_plot, labels = c("A", "B"))
save_plot("figures/unlogged_log.png", both_plots,
         ncol = 2, # we're saving a grid plot of 2 columns
         nrow = 1, # and 2 rows
         # each individual subplot should have an aspect ratio of 1.3
         base_aspect_ratio = 1.3
)

#produce a loop that runs all models and then gives model selection tables and parameter estimates as an output
#and then puts the residuals into dataframe for which a spatial correlogram is run
pdf("Figures/Abundance_residuals.pdf")
model_sel_summary<-NULL
for (i in seq(grep("SpR", (colnames(FD_comp))),grep("Rao", colnames(FD_comp)))){
  #run null models to check which random effects structure has the best fit
  M0<-lmer(FD_comp[[i]]~1+(1|Study),data=FD_comp,REML=F)
  M1<-lmer(FD_comp[[i]]~1+log(Age)+(1|Study),data=FD_comp,REML=F)
  model_sel<-data.frame(Variable=names(FD_comp[i]),
                        x_var=c("Null model","log(Age)"),
                        AICc=AICc(M0,M1)$AICc)
  model_sel$delta<-max(model_sel$AICc)-min(model_sel$AICc)
  model_sel$R2<-c(r.squaredGLMM(M0)[1],r.squaredGLMM(M1)[1])
  mod_resid<-data.frame(resids=c(resid(M0),resid(M1)),fitted=c(fitted(M0),fitted(M1)),model=rep(x = c("M0","M1"),each=43),title=names(FD_comp[i]))
  Title<-names(FD_comp[i])
  print(ggplot(mod_resid,aes(x=resids))+geom_histogram()+facet_wrap(~model)+ ggtitle(Title))
  print(Resid_plots<-ggplot(mod_resid,aes(x=fitted,y=resids))+geom_point()+facet_wrap(~model)+ ggtitle(Title))
  model_sel_summary<-rbind(model_sel_summary,model_sel)
}
dev.off()

write.csv(model_sel_summary,"Tables/model_sel_summary.csv")

#variables that respond to age - SpR, FDiv
#variables that don't respond to age - Even_comp, FDpg_comp, FDw_comp, FE_comp, FDis_comp, Rao_comp



ggplot(FD_comp,aes(x=log(Age),y=FRic))+geom_point()+geom_smooth(se=F,method="lm")
ggplot(FD_comp,aes(x=log(Age),y=FDiv))+geom_point()+geom_smooth(se=F,method="lm")


#first look at variables that respond to age
#Species richness
M1_SPR<-lmer(SpR~log(Age)+(1|Study),data=FD_comp)
summary(M1_SPR)
r.squaredGLMM(M1_SPR)

SpR_coefs <- data.frame(coef(summary(M1_SPR)))
mixed(SpR~log(Age)+(1|Study),data=FD_comp,test.intercept =T)
summary(SPR_P)
mixed(FDiv~log(Age)+(1|Study),data=FD_comp,test.intercept =T)


# use normal distribution to approximate p-value
SpR_coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
SpR_coefs$variable<-"SpR"

#FDiv
M1_FDiv<-lmer(FDiv~log(Age)+(1|Study),data=FD_comp)
summary(M1_FDiv)
r.squaredGLMM(M1_FDiv)
FDiv_coefs <- data.frame(coef(summary(M1_FDiv)))
# use normal distribution to approximate p-value
FDiv_coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
FDiv_coefs$variable<-"FDiv"

#plot results
new.data_FDiv<-data.frame(Age=seq(min(FD_comp$Age),max(FD_comp$Age)))
new.data_FDiv$FDiv<-predict(M1_FDiv,new.data_FDiv,re.form=NA)
new.data_FDiv$SE<-predict(M1_FDiv,new.data_FDiv,re.form=NA,se.fit=T)$se.fit


new.data_SpR<-data.frame(Age=seq(min(FD_comp$Age),max(FD_comp$Age)))
new.data_SpR$SpR<-predict(M1_SPR,new.data_SpR,re.form=NA)
new.data_SpR$SE<-predict(M1_SPR,new.data_SpR,re.form=NA,se.fit=T)$se.fit

#plot for species richness
theme_set(theme_cowplot())
SPR_P1<-ggplot(FD_comp,aes(x=Age,y=SpR))+geom_point(shape=1)+scale_x_log10()
SPR_P2<-SPR_P1+geom_line(data=new.data_SpR,aes(x=Age,y=SpR),size=1)+xlab("Time since last disturbance (Years)")
SPR_P3<-SPR_P2+geom_hline(yintercept=0,lty=2)+ylab("Relative secondary forest \nspecies richness (response ratio)")
SPR_P4<-SPR_P3+geom_ribbon(data=new.data_SpR,aes(x=Age,y=SpR,ymax=SpR+(2*SE),ymin=SpR-(2*SE)),alpha=0.5)

#plot for FDiv
FDiv_P1<-ggplot(FD_comp,aes(x=Age,y=FDiv))+geom_point(shape=1)+scale_x_log10()
FDiv_P2<-FDiv_P1+geom_line(data=new.data_FDiv,aes(x=Age,y=FDiv),size=1)+xlab("Time since last disturbance (Years)")
FDiv_P3<-FDiv_P2+geom_hline(yintercept=0,lty=2)+ylab("Relative secondary forest \nfunctional divergence (response ratio)")
FDiv_P4<-FDiv_P3+geom_ribbon(data=new.data_FDiv,aes(x=Age,y=FDiv,ymax=FDiv+(2*SE),ymin=FDiv-(2*SE)),alpha=0.5)

#################################################
#now look at variables that don't respond to age#
#################################################

#variables that don't respond to age - FDpg, FDw, FE, FDis, Rao
head(FD_comp)
coefs_summary<-NULL
var_list<-c("FDpg","FRic","FEve","FDis")
for (i in 1:length(var_list)){
  #i<-2
  j<-grep(var_list[i], colnames(FD_comp))
  M1<-lmer(FD_comp[[j]]~1+(1|Study),data=FD_comp)
  print(var_list[i])
  FD_comp$FDpg
  
  coefs <- data.frame(coef(summary(M1)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs$variable<-var_list[i]
  coefs_summary<-rbind(coefs,coefs_summary)
}

Coefs_summary<-rbind(coefs_summary,FDiv_coefs,SpR_coefs)

coefs_summary$var2<-c("Functional \ndispersion \n(FDis)","Functional \nEvenness \n(FEve)","Functional \nrichness \n(FRic)","Functional \nDiversity \n(FD)")

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

#produce models of this relationship
M0<-lmer(Prop_rich~1+(1|Study),data=SF_summary)
M1<-lmer(Prop_rich~log(Age)+(1|Study),data=SF_summary)

AICc(M0,M1)
summary(M1)

r.squaredGLMM(M1)

#plot results
new.data<-data.frame(Age=seq(min(SF_summary$Age),max(SF_summary$Age),length.out = 1000))
new.data$Prop_rich<-predict(M1,new.data,re.form=NA,se.fit=T)$fit
new.data$SE<-predict(M1,new.data,re.form=NA,se.fit=T)$se.fit


ForSpec_P1<-ggplot(SF_summary,aes(x=Age,y=Prop_rich))+geom_point(shape=1)+scale_x_log10()
ForSpec_P2<-ForSpec_P1+geom_line(data=new.data,aes(x=Age,y=Prop_rich),size=1)+xlab("Time since last disturbance (Years)")
ForSpec_P3<-ForSpec_P2+geom_hline(yintercept=0,lty=2)+ylab("Relative secondary forest \nspecialist richness (response ratio)")
ForSpec_P4<-ForSpec_P3+geom_ribbon(data=new.data,aes(x=Age,y=Prop_rich,ymax=Prop_rich+(2*SE),ymin=Prop_rich-(2*SE)),alpha=0.5)


time_plots<-plot_grid(SPR_P4, ForSpec_P4,FDiv_P4, labels = c("(a)", "(b)","(c)"), align = "h",ncol = 1)
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
