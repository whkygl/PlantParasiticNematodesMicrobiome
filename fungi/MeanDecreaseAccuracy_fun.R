library(tidyverse)
imp <- read.delim("f_genus_imp.txt")
# filter挑选Phylum列为非空的行
imp <- imp %>% filter(Phylum != "")
mycol <- c("#E64B35","#4DBBD5","#00A087","#3C5488",
           "#F39B7F","#4f72d2","#15a382","#810681","#c06e26",
           "#B09C85","#b5426a","#F2A900","#789048","#DA70BF",
           "#02468f","#b89a2b","#008856","#E68310","#222222","#96014e")
p <- ggplot(imp, aes(x = reorder(Genus,-MeanDecreaseAccuracy), y = MeanDecreaseAccuracy, fill = Phylum)) +
  geom_col() + theme_classic() + scale_y_continuous(expand=c(0, 0)) +
  #将调色板改为mycol
  scale_fill_manual(values = mycol)+
  #去掉x轴标签和x轴标题
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.5,0.8)) +
  #去掉图例标题
  guides(fill = guide_legend(title = NULL))+
  #修改y轴标题为“Mean decrease in accuracy”
  ylab("Mean decrease in accuracy")

p