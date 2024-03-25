rm(list=ls())
library(tidyverse)
bacmap <- read.csv("bacmap.csv")
# filter删除metadata中plant为Pine的行
bacmap <- bacmap %>% filter(bacmap$Plant != "Pine")
# 获取20个nature杂志风格的颜色
mycol <- c("#E64B35","#4DBBD5","#00A087","#3C5488",
    "#F39B7F","#4f72d2","#15a382","#DC0000","#c06e26",
    "#B09C85","#b5426a","#F2A900","#789048","#DA70BF",
    "#02468f","#b89a2b","#008856","#E68310","#222222","#96014e")
world_map <- map_data("world")
map <- ggplot() +
          scale_color_manual(values = mycol)+
          geom_polygon(data=world_map, aes(x = long, y = lat, group = group),fill="#c4c4c4", colour = "white",linewidth=0.15) +
          geom_point(data=bacmap,aes(x=Longitude,y = Latitude,colour=Plant,shape=Pathogen),size = 3)+
          theme_light()+theme(axis.title = element_blank())+
        scale_shape_manual(values=c(15,16,17,18,8,11))+
        #去除坐标轴
        theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
              axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
              legend.title=element_blank(),
              legend.position = c(.1,.5),
              legend.background = element_blank(),
              legend.key = element_blank())
map
ggsave("bac_map.pdf",width = 9,height = 6.75,units = "in")
                       