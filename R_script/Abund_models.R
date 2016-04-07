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
head(FD_comp)
str(FD_comp)


FD_comp$Lat<-FD_comp$Lat+(rnorm(length(FD_comp$Lat),0,0.00001)) 
#carry out analyses to see if there is any evidence of spatial autocorrelation for effect sizes
spat_cor_summary<-NULL
for (i in seq(grep("SpR_comp", (colnames(FD_comp))),grep("Rao_comp", colnames(FD_comp)))){
spat_cor<-spline.correlog(FD_comp$Lat, FD_comp$Long, FD_comp[,i], latlon=T,resamp=1000,quiet=T,na.rm=T)
spat_cor2<-data.frame(Dist=spat_cor$boot$boot.summary$predicted$x[1,],
                      Cor=spat_cor$boot$boot.summary$predicted$y[6,],
                      UCI=spat_cor$boot$boot.summary$predicted$y[2,],
                      LCI=spat_cor$boot$boot.summary$predicted$y[10,],
                      variable=names(FD_comp)[i])
spat_cor_summary<-rbind(spat_cor_summary,spat_cor2)
print(i)
}



theme_set(theme_bw(base_size=12))
Moran_plot1<-ggplot(spat_cor_summary,aes(x=Dist,y=Cor,ymax=LCI,ymin=UCI))+geom_ribbon(alpha=0.2)+geom_line(size=0.5,colour="black")+facet_wrap(~variable,scale="free")
Moran_plot2<-Moran_plot1+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))+ theme(legend.position="none")
Moran_plot2+geom_hline(y=0,lty=2)+xlab("Distance between sites (km)")+ylab("Moran's I correlation")+scale_size_continuous(range = c(1,3))








#species level metrics

#test which arrangement of random effects is best
AICc_summary<-NULL
for (i in seq(grep("SpR_comp", (colnames(FD_comp))),grep("Rao_comp", colnames(FD_comp)))){
  #run null models to check which random effects structure has the best fit
  M0_1<-lmer(FD_comp[[i]]~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp,REML=F)
  M0_2<-lmer(FD_comp[[i]]~1+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp,REML=F)
  M0_3<-lmer(FD_comp[[i]]~1+(1|Mist_nets)+(1|Transect)+(1|Study),data=FD_comp,REML=F)
  M0_4<-lmer(FD_comp[[i]]~1+(1|Mist_nets)+(1|Study),data=FD_comp,REML=F)
  M0_5<-lmer(FD_comp[[i]]~1+(1|Study),data=FD_comp,REML=F)
  AICc_sum<-data.frame(Variable=names(FD_comp[i]),Ran_effects=c(5,4,3,2,1),AICc=AICc(M0_1,M0_2,M0_3,M0_4,M0_5)$AICc)
  AICc_summary<-rbind(AICc_summary,AICc_sum) 
}
#without exception the model with the lowest number of random effects comes out on top

#produce a loop that runs all models and then gives model selection tables and parameter estimates as an output
#and then puts the residuals into dataframe for which a spatial correlogram is run
model_sel_summary<-NULL
for (i in seq(grep("SpR_comp", (colnames(FD_comp))),grep("Rao_comp", colnames(FD_comp)))){
  #run null models to check which random effects structure has the best fit
  M0<-lmer(FD_comp[[i]]~1+(1|Study),data=FD_comp,REML=F)
  M1<-lmer(FD_comp[[i]]~1+Age+(1|Study),data=FD_comp,REML=F)
  M2<-lmer(FD_comp[[i]]~1+log(Age)+(1|Study),data=FD_comp,REML=F)
  model_sel<-data.frame(Variable=names(FD_comp[i]),x_var=c("Null model","Age","log(Age)"),AICc=AICc(M0,M1,M2)$AICc)
  model_sel_summary<-rbind(model_sel_summary,model_sel) 
}
#variables that respond to age - SpR_comp, Shan_comp, FR_comp, FDiv_comp
#variables that don't respond to age - Even_comp, FDpg_comp, FDw_comp, FE_comp, FDis_comp, Rao_comp

M1<-lmer(SpR_comp~1+Age+(1|Study),data=FD_comp)
M2<-lmer(SpR_comp~1+log(Age)+(1|Study),data=FD_comp)
M3<-lmer(SpR_comp~1+PF_SpR+(1|Study),data=FD_comp)
AICc(M1,M2,M3)

summary(M3)
plot(M2)
plot(FD_comp$PF_SpR,FD_comp$SpR_comp)
points(FD_comp$Age,predict(M2,levels=0,re.form=NA),col="red")

plot(FD_comp$Age,FD_comp$Shan_comp)
plot(FD_comp$Age,FD_comp$FDiv_comp)


spat_cor<-spline.correlog(FD_comp$Lat, FD_comp$Long, resid(M2), latlon=T,resamp=1000,quiet=T,na.rm=T)
spat_cor2<-data.frame(Dist=spat_cor$boot$boot.summary$predicted$x[1,],
                      Cor=spat_cor$boot$boot.summary$predicted$y[6,],
                      UCI=spat_cor$boot$boot.summary$predicted$y[2,],
                      LCI=spat_cor$boot$boot.summary$predicted$y[10,],
                      variable=names(FD_comp)[i])
theme_set(theme_bw(base_size=12))
Moran_plot1<-ggplot(spat_cor2,aes(x=Dist,y=Cor,ymax=LCI,ymin=UCI))+geom_ribbon(alpha=0.2)+geom_line(size=0.5,colour="black")
Moran_plot2<-Moran_plot1+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))+ theme(legend.position="none")
Moran_plot2+geom_hline(y=0,lty=2)+xlab("Distance between sites (km)")+ylab("Moran's I correlation")+scale_size_continuous(range = c(1,3))

spat_cor_summary<-rbind(spat_cor_summary,spat_cor2)
print(i)

#Species richness
M0<-lmer(SpR_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M1<-lmer(SpR_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(SpR_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)


#Shannon diversity
M0<-lmer(Shan_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M1<-lmer(Shan_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(Shan_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)
summary(M1)

#Community evenness
M0<-lmer(Even_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M1<-lmer(Even_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(Even_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)


#Functional diversity - Petchy & Gaston method
M0_1<-lmer(FDpg_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp,REML=F)
M0_2<-lmer(FDpg_comp~1+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp,REML=F)
M0_3<-lmer(FDpg_comp~1+(1|Mist_nets)+(1|Transect)+(1|Study),data=FD_comp,REML=F)
M0_4<-lmer(FDpg_comp~1+(1|Mist_nets)+(1|Study),data=FD_comp,REML=F)
M0_5<-lmer(FDpg_comp~1+(1|Study),data=FD_comp,REML=F)

AICc(M0_1,M0_2,M0_3,M0_4,M0_5)


M1<-lmer(FDpg_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(FDpg_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)

#Functional diversity - Petchy & Gaston method weighted by abundance
M0<-lmer(FDw_comp~1+(1|Study),data=FD_comp)
M1<-lmer(FDw_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Study),data=FD_comp)
M2<-lmer(FDw_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)

#Functional richness
M0<-lmer(FR_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),REML=F,data=FD_comp)
M1<-lmer(FR_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
M2<-lmer(FR_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp)
AICc(M0,M1,M2)
summary(M0)
#functional richness is considerably lower in secondary forests - but the errors are large

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


#now look at community weighted means
head(FD_comp)
FD_comp_CWD<-FD_comp[,c(1:8,59:ncol(FD_comp))]
#remove infinite values
is.na(FD_comp_CWD) <- sapply(FD_comp_CWD, is.infinite)
FD_comp_CWD[is.na(FD_comp_CWD)] <- 0

#reorganise the data
FD_comp_CWD_melt<-melt(FD_comp_CWD,id.vars=c("SiteID","Study","Age","PF_SF","Point_obs","Mist_nets","Transect","Vocal"))
ggplot(FD_comp_CWD_melt, aes(x=variable, y=value)) + stat_summary(fun.y="mean", geom="point")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#this looks interesting, now lets test this statistically
head(FD_comp_CWD)

#Invertebrate feeding species
M0<-lmer(Diet.Inv_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
M1<-lmer(Diet.Inv_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
M2<-lmer(Diet.Inv_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
AICc(M0,M1,M2)
summary(M0)

#Vegetative (?) feeding species
M0<-lmer(Diet.Vend_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
M1<-lmer(Diet.Vend_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
M2<-lmer(Diet.Vend_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
AICc(M0,M1,M2)
summary(M0)

#Body size
M0<-lmer(BodyMass.Value_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
M1<-lmer(BodyMass.Value_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
M2<-lmer(BodyMass.Value_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
AICc(M0,M1,M2)
summary(M0)

#aerial foraging
M0<-lmer(ForStrat.aerial_comp~1+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
M1<-lmer(ForStrat.aerial_comp~1+Age+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
M2<-lmer(ForStrat.aerial_comp~1+log(Age)+(1|Point_obs)+(1|Mist_nets)+(1|Transect)+(1|Vocal)+(1|Study),data=FD_comp_CWD)
AICc(M0,M1,M2)
summary(M0)


head(FD_comp_CWD)
