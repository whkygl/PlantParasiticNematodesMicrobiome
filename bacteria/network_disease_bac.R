rm(list = ls())
library(tidyverse)
library(readxl)
library(vegan)
library(randomForest)
library(ROCR)
library(caret)
library(ggpubr)
library(psych)
library(circlize)
library(viridis)
library(igraph)
library(ggraph)
library(colormap)
library(wesanderson)
metadata <- read_excel("metadata_all.xlsx")
bacteria <- read.delim("bacteria.txt")
metadata <- metadata %>% filter(parts=="RS" & is.na(intervene))
metadata <- metadata %>% filter(plant != "Pine")
bac <- bacteria
id <- intersect(colnames(bac),metadata$NO)
bac <- bac[,c("OTUID",id,"taxonomy")]

rownames(bac) <- bac[,1]
bac <- bac[,-1]
bac <- bac[,-1310]
colSums(bac)
# 按2000抽平
set.seed(123)
bac1 <- as.data.frame(t(rrarefy(t(bac),2000)))
# 删除otu_rare中列和小于2000的列
bac1 <- bac1[,colSums(bac1)>=2000]
bac <- bac1
# 为bac1分配分类
taxonomy <- bacteria %>% select(OTUID, taxonomy)
bac1$OTUID <- row.names(bac1)
bac1 <- left_join(bac1, taxonomy, by = "OTUID")
bac1 <- bac1 %>% select(OTUID, everything())
bac <- bac1
id <- intersect(colnames(bac),metadata$NO)
# filter选取metadata中NO为id的行
metadata <- metadata %>% filter(NO %in% id)
bacteria <- bac
bacteria[,799] <- gsub(",s.*","",bacteria[,799])
# 若bac最后一列包含g:则只保留g:后面的内容，否则改为Unassign
bacteria[,799] <- ifelse(grepl("g:",bacteria[,799]),gsub(".*g:","",bacteria[,799]),"Unassign")
rownames(bacteria) <- bacteria[,1]
bacteria <- bacteria[,-1]
bacteria <- bacteria %>% group_by(taxonomy) %>% summarise_all(sum)
bacteria <- as.data.frame(bacteria)




rownames(bacteria) <- bacteria[,1]
bacteria <- bacteria[,-1]
# arrange将metadata中的NO列按照bacteria的列名重排
idc= colnames(bacteria)
idc
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$NO
metadata <- metadata[idc,]
bac <- bacteria
bac_t <- as.data.frame(t(bac))
bac_t$group <- factor(metadata$plant_status,levels = c("D","H"))
#列名中的-改为_
colnames(bac_t) <- gsub("-","_",colnames(bac_t))
imp <- read.delim("b_genus_imp.txt")
vertices <- read.delim("b_genus_imp.txt")
set.seed(123)
# 将行随机排列
vertices <- vertices[sample(1:nrow(vertices)),]
#filter选取Phylum不为""的行
id2 <- c(imp$Genus,"group")
# 选取bac_t中列名为id2的列
bac_t <- bac_t %>% select(id2)

#差异比较
# 每个值除以该行和
# bac_t <- as.data.frame(bac_t)
# bac_t2 <- bac_t[,1:24]
# bac_t2_prob <- as.data.frame(t(apply(bac_t2, 1, function(x) x / sum(x))))
# bac_t2_prob$group <- bac_t$group
bac_d <- bac_t%>% filter(group == "D")
#删除group列
bac_d <- bac_d[,-25]
#将bac_d的列随机排列
set.seed(123)
bac_d <- bac_d[,sample(1:24)]
occor = corr.test(bac_d,use="pairwise",method="spearman",adjust="BH",alpha=.05) 
occor.r = occor$r
occor.p = occor$p
occor.r[occor.p>0.05 | abs(occor.r)<0.4]=0
#将occor.r的行名变为第一列
occor.r <- as.data.frame(occor.r)
occor.r$from <- rownames(occor.r)
occor.r <- occor.r %>% select(from,everything())
dataUU <- occor.r
dataUU[dataUU==0] <- NA
dataUU[,-1][dataUU[,-1]>0] <- 1
dataUU[dataUU<0] <- -1
connect <- dataUU %>% 
  gather(key="to", value="value", -1) %>%
  mutate(to = gsub("\\.", "",to)) %>%
  na.omit()
connect <- connect %>% 
  mutate(group1 = ifelse(value==1, "p", "n"))
mygraph <- graph_from_data_frame(connect, vertices = vertices, directed = FALSE )
com <- walktrap.community(mygraph)
number_of_bar<-nrow(vertices)
vertices$id = seq(1, nrow(vertices))
angle= 360 * (vertices$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
vertices$hjust<-ifelse(angle>180, 1, 0)
vertices$angle<-ifelse(angle>180, 90-angle+180, 90-angle)
# Create a graph object with igraph
mygraph <- graph_from_data_frame( connect, vertices = vertices, directed = FALSE) 
occor.r2 <- occor.r[,-1]
mygraph2 <- graph.adjacency(as.matrix(occor.r2), mode = "undirected", weighted = TRUE, diag = FALSE)
mycolor = c("#3B4992FF","#008B45FF")
#(b)曲线链接
p <- ggraph(mygraph,layout = 'linear', circular = TRUE) +
  geom_edge_arc(aes(edge_colour=as.factor(group1)), edge_alpha=0.3, edge_width=0.7) +
  geom_node_point(aes(fill=as.factor(color), size = MeanDecreaseAccuracy), shape=21,color='black',alpha=1) +
  scale_size_continuous(range=c(1,7)) +
  scale_fill_manual(values=mycolor) +
  geom_node_text(aes(x = x*1.06, y=y*1.06, label=name, angle=angle,hjust=hjust), color = "black",size=3) +
  scale_color_manual(values=mycolor) +
  scale_edge_color_manual(values=mycolor) +
  expand_limits(x = c(-1.6, 1.6), y = c(-1.6, 1.6))+
  coord_fixed()+
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks =element_blank(),
    axis.text =element_blank(),
    axis.title = element_blank(),
    plot.margin=unit(c(0,0,0,0), "null"),
    panel.spacing=unit(c(0,0,0,0), "null")
  ) + theme(legend.position="none")
p
# 计算拓扑属性
library(ggClusterNet)
dat = net_properties.4(mygraph2,n.hub = F)
write.csv(dat,"net_properties.csv")