rm(list=ls())

#load packages needed for this work
library(ggplot2)
library(nlme)
library(MuMIn)
library(gridExtra)

# load data
setwd("C:/Users/Phil/Dropbox/Phil & Catherine - Secondary forests/Analysis - Dec 2013/Data")
setwd("C:/Users/Phil/Dropbox/Work/Active projects/Catherine Sayer/Phil & Catherine - Secondary forests (1)/Analysis - Dec 2013/Data")

data<-read.csv("FD Data To Use.csv",header=T)
names(data)

# plot data
plot(data$Age,data$Prop.SR) # SF age vs. proportion primary forest species richness
plot(data$Age,data$Prop.SSR) # SF age vs. proportion primary forest singular species richness
plot(data$Age,data$Prop.FRic) # SF age vs. proportion primary forest functional richness
plot(data$Age,data$Prop.Feve) # SF age vs. proportion primary forest functional eveness
plot(data$Age,data$Prop.Fdiv) # SF age vs. proportion primary forest functional diversity



# remove sites without continuity with PF data
data2<-subset(data,data$PFCont!="")
#reset factor levels (r has a habit of leaving in empty factor levels when you subset)
data2$PFCont<-factor(data2$PFCont)

###############################################################################################
#notes for Catherine###########################################################################
#I have made some edits to your code, while doing my work on secondary forests I realised that
#we shouldn't include both Age and log(Age) in the same model as they are just different transformations
#of the same variable. Sorry I didn't tell you this earlier. Also I have changed how I do model simplification
#and have now moved towards a model averaged approach, this type of modelling is much more robust than 
#stepwise simplification since it doesn't assume that there is one correct model but rather weights
#models by how parsimonous they are (a balance between model simplicity and explanatory value)
#I think it would be useful to use this approach here so I have edited your code. I have also 
#annotated it so it should be fairly easy to understand....


# Singular species richness model
#first we should look at different ways in which random effects can be defined


ggplot(data2,aes(x=Age,y=Prop.SSR))+geom_point()+geom_smooth(method="lm")+facet_wrap(~StudyID)
#to me from this plot it looks like the slopes of relationships with age may vary by study, we can put this in our random effects

#start with a null model (logging the response gives better residuals in this case)
SSR0.1<-lme(log(Prop.SSR)~1,random=(~1|StudyID),method="REML",data=data2)
SSR0.2<-lme(log(Prop.SSR)~1,random=(~Age|StudyID),method="REML",data=data2)
AICc(SSR0.1,SSR0.2)

#the model just with random intercepts seems the best - it has the lowest AICc
#now come up with all the models we want to test - in effect each one represents a hypothesis

SSR0<-lme(log(Prop.SSR)~1,random=(~1|StudyID),method="REML",data=data2)
SSR1<-lme(log(Prop.SSR)~Age+Dist2+PFCont,random=(~1|StudyID),method="ML",data=data2)
SSR2<-lme(log(Prop.SSR)~log(Age)+Dist2+PFCont,random=(~1|StudyID),method="ML",data=data2)
SSR3<-lme(log(Prop.SSR)~log(Age)+PFCont,random=(~1|StudyID),method="ML",data=data2)
SSR4<-lme(log(Prop.SSR)~Age+PFCont,random=(~1|StudyID),method="ML",data=data2)
SSR5<-lme(log(Prop.SSR)~log(Age)+Dist2,random=(~1|StudyID),method="ML",data=data2)
SSR6<-lme(log(Prop.SSR)~Age+Dist2,random=(~1|StudyID),method="ML",data=data2)
SSR7<-lme(log(Prop.SSR)~log(Age),random=(~1|StudyID),method="ML",data=data2)
SSR8<-lme(log(Prop.SSR)~Age,random=(~1|StudyID),method="ML",data=data2)

#plot residuals and qqplots of all these models

grid.arrange(plot(SSR0),plot(SSR1),plot(SSR2),plot(SSR3),plot(SSR4),plot(SSR5),plot(SSR6),plot(SSR7),plot(SSR8),
             qqnorm(SSR0),qqnorm(SSR1),qqnorm(SSR2),qqnorm(SSR3),qqnorm(SSR4),qqnorm(SSR5),qqnorm(SSR6),qqnorm(SSR7),qqnorm(SSR8),ncol=9)
#the resdiduals are not fantastic, models look like the tend to over predict for low predictions and under predict for high predications
#however this is probably as good as we'll get it at the moment


#now we put these models in a list
All_mods<-list(SSR0,SSR1,SSR2,SSR3,SSR4,SSR5,SSR6,SSR7,SSR8)
#and produce a model selection table
ms1<-mod.sel(object=All_mods,rank="AICc",fit=T,trace=T)
#and add the rsquared to this table
ms1$r2<-c(r.squaredGLMM(SSR7)[1],r.squaredGLMM(SSR8)[1],r.squaredGLMM(SSR5)[1],r.squaredGLMM(SSR3)[1],r.squaredGLMM(SSR0)[1],r.squaredGLMM(SSR4)[1],r.squaredGLMM(SSR6)[1],r.squaredGLMM(SSR2)[1],r.squaredGLMM(SSR1)[1])
#you can see from the r squared values these models have very little explanatory power

#now get models with a delta AICc<7, this is the standard way to drop model that have little statistical support
delta7<- get.models(ms1, subset = cumsum(delta) <= 7)
#and average the coefficients from these models
avgm <- model.avg(delta7,se.fit=T)

#now we create a dataset with all variables

summary(data2)

new_data<-expand.grid(Age=seq(1,100),Dist2=c("Y","N"),PFCont=c("No","Yes"))

Preds<-predict(avgm,newdata=new_data,level=0,se.fit=T)
#bind data and predictions together for use in figures
Preds2<-cbind(new_data,Preds)
head(Preds2)


#plots of this result
theme_set(theme_bw(base_size=12))
Plot1<-ggplot(data2,aes(x=Age,y=Prop.SSR))+geom_point()+geom_line(data=Preds2,aes(x=Age,y=exp(fit),group=interaction(Dist2,PFCont),colour=interaction(Dist2,PFCont)))
Plot2<-Plot1+geom_line(data=Preds2,aes(x=Age,y=exp(fit+(1.96*se.fit)),group=interaction(Dist2,PFCont),colour=interaction(Dist2,PFCont)),lty=2)+geom_line(data=Preds2,aes(x=Age,y=exp(fit-(1.96*se.fit)),group=interaction(Dist2,PFCont),colour=interaction(Dist2,PFCont)),lty=2)
Plot3<-Plot2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))+geom_hline(y=1,lty=2,size=2,colour="grey")
Plot3+ scale_colour_discrete(name="Contiguous with forests & ",
                            breaks=c("ctrl", "trt1", "trt2"),
                            labels=c("Control", "Treatment 1", "Treatment 2"))
