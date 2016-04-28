#this script is for analysis of differences in species richness and functional diversity metrics
#between primary and secondary forest using data for which only species presence/absence is available

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


#load in data
rm(list = ls())
FD_comp<-read.csv("Data/FD_summary_comp.csv")
Location<-read.csv("Data/Study_location.csv")
FD_comp<-merge(FD_comp,Location,by.x="Site",by.y="SiteID")

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
for (i in seq(grep("SpR_comp", (colnames(FD_comp))),grep("FDis_comp", colnames(FD_comp)))){
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
write.csv(AICc_summary,"Tables/PA_Random_model_AICc.csv")

#produce a loop that runs all models and then gives model selection tables and parameter estimates as an output
#and then puts the residuals into dataframe for which a spatial correlogram is run
pdf("Figures/PA_residuals.pdf")
model_sel_summary<-NULL
for (i in seq(grep("SpR_comp", (colnames(FD_comp))),grep("FDis_comp", colnames(FD_comp)))){
  #run null models to check which random effects structure has the best fit
  M0<-lmer(FD_comp[[i]]~1+(1|Study),data=FD_comp,REML=F)
  M1<-lmer(FD_comp[[i]]~1+log(Age)+(1|Study),data=FD_comp,REML=F)
  model_sel<-data.frame(Variable=names(FD_comp[i]),
                        x_var=c("Null model","log(Age)"),
                        AICc=AICc(M0,M1)$AICc)
  model_sel$delta<-max(model_sel$AICc)-min(model_sel$AICc)
  model_sel$R2<-c(r.squaredGLMM(M0)[1],r.squaredGLMM(M1)[1])
  mod_resid<-data.frame(resids=c(resid(M0),resid(M1)),fitted=c(fitted(M0),fitted(M1)),model=rep(x = c("M0","M1"),each=44),title=names(FD_comp[i]))
  Title<-names(FD_comp[i])
  print(ggplot(mod_resid,aes(x=resids))+geom_histogram()+facet_wrap(~model)+ ggtitle(Title))
  print(Resid_plots<-ggplot(mod_resid,aes(x=fitted,y=resids))+geom_point()+facet_wrap(~model)+ ggtitle(Title))
  model_sel_summary<-rbind(model_sel_summary,model_sel)
}
dev.off()

write.csv(model_sel_summary,"Tables/PA_model_sel_summary.csv")

#variables that respond to age - FEve, FDiv
#variables that don't respond to age - Even_comp, FDpg_comp, FDw_comp, FE_comp, FDis_comp, Rao_comp


#first look at variables that respond to age
#FDiv
M1_FE<-lmer(FE_comp~log(Age)+(1|Study),data=FD_comp)
summary(M1_FE)
plot(M1_FE)
r.squaredGLMM(M1_FDiv)

FD_comp$Age_trans<-(log(FD_comp$Age)-mean(log(FD_comp$Age)))/sd(FD_comp$Age)

FE_comp_M1<-lme(FE_comp~Age_trans,random=~1|Study,weights=varExp(form=~Study),data=FD_comp)

summary(FE_comp_M1)


r.squaredGLMM(FE_comp_M1)

?lme

par(mfrow=c(1,1))
plot(FD_comp$FE_comp,predict(FE_comp_M1,re.form=NULL))

plot(fitted(FE_comp_M1),resid(FE_comp_M1))

weights=varExp(form=~logRe|Study)


FE_M1<-lm(FDiv_comp~log(Age),data=FD_comp)
summary(FE_M1)
par(mfrow=c(2,2))
plot(FE_M1)

#first look at variables that respond to age
#FDiv
M1_FDiv<-lmer(FDiv_comp~log(Age)+(1|Study),data=FD_comp)
summary(M1_FDiv)
plot(M1_FDiv)
r.squaredGLMM(M1_FDiv)


ggplot(FD_comp,aes(x=log(Age),y=FDiv_comp))+geom_point()+geom_abline(intercept=0.028323,slope=-0.009531)





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


?mcmcsamp

#plot results
new.data<-data.frame(Age=seq(min(FD_comp$Age),max(FD_comp$Age)))
new.data$FDiv<-predict(M1_FDiv,new.data,re.form=NA)
new.data$SpR<-predict(M1_SPR,new.data,re.form=NA)
new_data_preds<-melt(new.data,id.vars=c("Age"))
new_data_preds$SE<-c(predict(M1_FDiv,new.data,re.form=NA,se.fit=T,)$se.fit,
                     predict(M1_SPR,new.data,re.form=NA,se.fit=T)$se.fit)

new_data_preds$var2<-ifelse(new_data_preds$variable=="SpR","Species richness","Functional divergence (FDiv)")


head(FD_comp)
FD_comp_plots<-FD_comp[,c(3,9,16)]
FD_comp_plots_melt<-melt(FD_comp_plots,id.vars="Age")
FD_comp_plots_melt$var2<-ifelse(FD_comp_plots_melt$variable=="SpR","Species richness","Functional divergence (FDiv)")
theme_set(theme_bw(base_size=12))
P1<-ggplot(FD_comp_plots_melt,aes(x=Age,y=value))+geom_point(shape=1)+facet_wrap(~var2,scales = "free_y",)+scale_x_log10()
P2<-P1+geom_line(data=new_data_preds,aes(x=Age,y=value),size=1)+xlab("Age (years - log scale)")
P3<-P2+geom_hline(yintercept=0,lty=2)+ylab("Difference between secondary \nand primary forest sites (log response ratio)")
P3+geom_ribbon(data=new_data_preds,aes(x=Age,y=value,ymax=value+(2*SE),ymin=value-(2*SE)),alpha=0.5)
ggsave("Figures/Age_models_abun.pdf",width=8,height=4,dpi=400,units="in")
ggsave("Figures/Age_models_abun.png",width=8,height=4,dpi=400,units="in")

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
  #M1<-lmer(FD_comp[[j]]~1+(1|Study),data=FD_comp)
  print(var_list[i])
  FD_comp$FDpg
  mixed(FDpg~1+(1|Study),data=FD_comp,test.intercept =T)
  mixed(FRic~1+(1|Study),data=FD_comp,test.intercept =T)
  mixed(FDpg~1+(1|Study),data=FD_comp,test.intercept =T)
  mixed(FDis~1+(1|Study),data=FD_comp,test.intercept =T)
  mixed(FDiv~1+(1|Study),data=FD_comp,test.intercept =T)
  
  #coefs <- data.frame(coef(summary(M1)))
  #coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  #coefs$variable<-var_list[i]
  #coefs_summary<-rbind(coefs,coefs_summary)
}

Coefs_summary<-rbind(coefs_summary,FDiv_coefs,SpR_coefs)

coefs_summary$var2<-c("Functional dispersion (FDis)","Functional Evenness (FEve)","Functional richness (FRic)","Functional Diversity (FD)")

write.csv(Coefs_summary,"Tables/Coefs_summary.csv")

Plot1<-ggplot(coefs_summary,aes(x=var2,y=Estimate,ymax=Estimate+(1.96*Std..Error),ymin=Estimate-(1.96*Std..Error)))+geom_point(shape=1,size=3)+geom_errorbar(width=0.5)
Plot2<-Plot1+geom_hline(yintercept=0,lty=2)+ylab("Difference between secondary \nand primary forset sites (log response ratio)")
Plot2+xlab("Diversity variable")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("Figures/Null_models_abun.pdf",width=6,height=6,dpi=400,units="in")
ggsave("Figures/Null_models_abun.png",width=6,height=6,dpi=400,units="in")

