rm(list = ls())
library(tidyverse)
library(readxl)
metadata <- read_excel("metadata.xlsx")
my_color <- c("#E64B35","#4DBBD5","#00A087","#3C5488",
    "#F39B7F","#4f72d2","#0096b0","#DC0000","#c06e26",
    "#B09C85","#b5426a","#F2A900","#222222","#DA70BF",
    "#02468f","#b89a2b","#008856","#E68310","#789048","#96014e")
#filter筛选intervene为缺失值的数据
metadata <- metadata %>% filter(is.na(intervene))
#筛选microbiome为B且parts为RS的行
metadata <- metadata %>% filter(microbiome == "F" & parts == "RS")
# filter删除metadata中plant为Pine的行
metadata <- metadata %>% filter(metadata$plant != "Pine")
#group_by分组统计plant列每个元素的个数
plant <- metadata %>% group_by(plant) %>% summarise(n = n())
plant$plant[which(plant$n<=12)]="Others"
plant <- plant %>% group_by(plant) %>% summarise(n = sum(n))
#添加配色列
plant$color <- my_color[1:nrow(plant)]
labs <- paste0(plant$plant," \n(", round(plant$n/sum(plant$n)*100,2), "%)")
#6X6
pie(plant$n,labels=labs, init.angle=90,col = plant$color,
    border="black")
pathogen <- metadata %>% group_by(pathogen) %>% summarise(n = n())
pathogen$pathogen[which(pathogen$n<=30)]="Others"
pathogen <- pathogen %>% group_by(pathogen) %>% summarise(n = sum(n))
pathogen$color <- my_color[1:nrow(pathogen)]
labs <- paste0(pathogen$pathogen," \n(", round(pathogen$n/sum(pathogen$n)*100,2), "%)")
pie(pathogen$n,labels=labs, init.angle=90,col = pathogen$color,
    border="black")
amplicon <- metadata %>% group_by(amplicon) %>% summarise(n = n())
amplicon$amplicon[which(amplicon$n<=15)]="Others"
amplicon <- amplicon %>% group_by(amplicon) %>% summarise(n = sum(n))
amplicon$color <- my_color[1:nrow(amplicon)]
labs <- paste0(amplicon$amplicon," \n(", round(amplicon$n/sum(amplicon$n)*100,2), "%)")
pie(amplicon$n,labels=labs, init.angle=90,col = amplicon$color,
    border="black")
country <- metadata %>% group_by(country) %>% summarise(n = n())
country$country[which(country$n<=16)]="Others"
country <- country %>% group_by(country) %>% summarise(n = sum(n))
country$color <- my_color[1:nrow(country)]
labs <- paste0(country$country," \n(", round(country$n/sum(country$n)*100,2), "%)")
pie(country$n,labels=labs, init.angle=90,col = country$color,
    border="black")
