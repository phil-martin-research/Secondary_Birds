# map using ggplot
install.packages("maps")
install.packages("ggplot2")

library(ggplot2)
library(maps)
site_data<-read.csv("Data/Site Locations.csv",header=T)

theme_set(theme_bw(base_size=12))
world<-map_data("world")
p1<- ggplot()
p2<- p1 + geom_polygon( data=world, aes(x=long, y=lat, group = group),fill="grey")#Create a base plot
p3 <- p2 + geom_point( data=site_data, aes(x=long, y=lat, size =Secondary.Sites), color="navyblue",alpha=0.5) + scale_size(name="No. of secondary forest sites",range=c(3,7))
p4<- p3 + coord_equal(xlim = c(-150,160),ylim=c(-40,40))
p5 <- p4 + theme(legend.title = element_text()) + theme(legend.text = element_text())
p6<-p5+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=1.5,colour="black",fill=NA),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
p6+coord_equal()


ggsave(filename = "Figures/Study_locations.pdf",width = 8,height=4,units='in',dpi=400)
ggsave(filename = "Figures/Study_locations.png",width = 8,height=4,units='in',dpi=400)
