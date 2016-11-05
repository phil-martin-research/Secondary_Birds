#this script is to produce analyses specifically to address the concerns of the reviewer
#none of this analysis will be included in the manuscript

library(ggplot2)
library(cowplot)

#first look at correlation of traits to address comments 14 and 15
Bird_traits<-read.csv("Data/Bird_traits_Oct_2016.csv")
head(Bird_traits)

#look at correlation between mass for HBW and EltonTraits
HBW_Elton_P1<-ggplot(Bird_traits,aes(x=BodyMass.Value,y=Body_mass_HBW))+geom_abline(size=2,lty=2,alpha=0.5)+geom_point(alpha=0.1)+scale_x_log10()+scale_y_log10()
HBW_Elton_P1+xlab("Bird body mass (g) taken from EltonTraits database")+ylab("Bird body mass (g) taken from HBW")
ggsave(filename = "Figures/Eltontraits_HBW.png",width = 8,height = 5,dpi = 100,units = "in")

#also look at correlation between Body length and mass
Length_weight_P1<-ggplot(Bird_traits,aes(x=Body_lenth_HBW,y=Body_mass_HBW))+geom_point(alpha=0.1)+scale_y_log10()+scale_x_log10()
Length_weight_P1+xlab("Bird body length (cm) taken from HBW")+ylab("Bird body mass (g) taken from HBW")+geom_smooth(se=F,method="lm")
ggsave(filename = "Figures/Length_mass.png",width = 8,height = 5,dpi = 400,units = "in")
