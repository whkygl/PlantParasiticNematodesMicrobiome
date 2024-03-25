rm(list = ls())
library(tidyverse)
library(ggforce) # 做桑基图转用包
library(readxl) # 读取excel文件
metadata_all <- read_xlsx("metadata_all.xlsx")
#filter筛选intervene为缺失值的数据
metadata <- metadata_all %>% filter(is.na(intervene))
metadata <- metadata %>% filter(metadata$plant != "Pine")
#筛选microbiome为B且parts为RS的行
metadata <- metadata %>% filter(microbiome == "B" & parts == "RS")
#选取paper,pathogen,plant,plant_status列
metadata <- metadata %>% select(paper,pathogen,plant,plant_status)
#将plant_status列的D改为Disease,H改为Health
metadata$plant_status <- ifelse(metadata$plant_status == "D","Disease","Health")
#group_by对这四列分组求个数
metadata <- metadata %>% group_by(paper,pathogen,plant,plant_status) %>% summarise(count = n())
#统计plant_status列中Disease的count总数
disease_count <- metadata %>% filter(plant_status == "Disease")
disease_count <- sum(disease_count$count)
#统计plant_status列中Health的count总数
health_count <- metadata %>% filter(plant_status == "Health")
health_count <- sum(health_count$count)
#在metadata的plant_status列中的Disease后添加(),并将disease_count的值写入括号内
metadata$plant_status <- ifelse(metadata$plant_status == "Disease",
    paste0("Disease(",disease_count,")"),metadata$plant_status)
#在metadata的plant_status列中的Health后添加(),并将health_count的值写入括号内
metadata$plant_status <- ifelse(metadata$plant_status == "Health",
    paste0("Health(",health_count,")"),metadata$plant_status)
my_color <- c("#E64B35","#4DBBD5","#00A087","#3C5488",
    "#F39B7F","#4f72d2","#15a382","#DC0000","#c06e26",
    "#B09C85","#b5426a","#F2A900","#222222","#DA70BF",
    "#02468f","#b89a2b","#008856","#E68310","#789048","#96014e")
bac_sankey <- metadata[, -1]
data <- gather_set_data(bac_sankey, 1:3)
#调整因子x的顺序
data$x <- factor(data$x, levels = c("3", "1", "2"))
ggplot(data, aes(x, id = id, split = y, value = count)) +
    geom_parallel_sets(aes(fill = pathogen), alpha = 0.7, axis.width = 0.05) +
    geom_parallel_sets_axes(axis.width = 0.05) +
    geom_parallel_sets_labels(colour = "black", angle = 0, position = position_nudge(x = .05, y = 0), hjust = 0) +
    theme_void() +
    guides(fill = "none") +
    scale_fill_manual(values = my_color)
ggsave("bac_sankey.pdf",width = 8.15,height = 6.55,units = "in")
