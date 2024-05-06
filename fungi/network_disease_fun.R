rm(list = ls())
library(tidyverse)
library(readxl)
library(vegan)
library(randomForest)
library(ROCR)
library("caret")
library(ggpubr)
library(psych)
library(circlize)
library(viridis)
library(igraph)
library(ggraph)
library(colormap)
library(wesanderson)
metadata <- read_excel("metadata_all.xlsx")
fungi <- read.delim("fungi.txt")
metadata <- metadata %>% filter(parts=="RS" & is.na(intervene))
metadata <- metadata %>% filter(plant != "Pine")
fun <- fungi
id <- intersect(colnames(fun),metadata$NO)
fun <- fun[,c("OTUID",id,"taxonomy")]

rownames(fun) <- fun[,1]
fun <- fun[,-1]
fun <- fun[,-1310]
colSums(fun)
# 按2000抽平
set.seed(123)
fun1 <- as.data.frame(t(rrarefy(t(fun),2000)))
# 删除otu_rare中列和小于2000的列
fun1 <- fun1[,colSums(fun1)>=2000]
fun <- fun1
# 为fun1分配分类
taxonomy <- fungi %>% select(OTUID, taxonomy)
fun1$OTUID <- row.names(fun1)
fun1 <- left_join(fun1, taxonomy, by = "OTUID")
fun1 <- fun1 %>% select(OTUID, everything())
fun <- fun1
id <- intersect(colnames(fun),metadata$NO)
# filter选取metadata中NO为id的行
metadata <- metadata %>% filter(NO %in% id)
fungi <- fun
fungi[,547] <- gsub(",s.*","",fungi[,547])
# 若fun最后一列包含g:则只保留g:后面的内容，否则改为Unassign
fungi[,547] <- ifelse(grepl("g:",fungi[,547]),gsub(".*g:","",fungi[,547]),"Unassign")
rownames(fungi) <- fungi[,1]
fungi <- fungi[,-1]
fungi <- fungi %>% group_by(taxonomy) %>% summarise_all(sum)
fungi <- as.data.frame(fungi)




rownames(fungi) <- fungi[,1]
fungi <- fungi[,-1]
# arrange将metadata中的NO列按照fungi的列名重排
idc= colnames(fungi)
idc
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$NO
metadata <- metadata[idc,]
fun <- fungi
fun_t <- as.data.frame(t(fun))
fun_t$group <- factor(metadata$plant_status,levels = c("D","H"))
#列名中的-改为_
colnames(fun_t) <- gsub("-","_",colnames(fun_t))
imp <- read.delim("b_genus_imp.txt")
vertices <- read.delim("b_genus_imp.txt")
set.seed(123)
# 将行随机排列
vertices <- vertices[sample(1:nrow(vertices)),]
#filter选取Phylum不为""的行
id2 <- c(imp$Genus,"group")
# 选取fun_t中列名为id2的列
fun_t <- fun_t %>% select(id2)

#差异比较
# 每个值除以该行和
# fun_t <- as.data.frame(fun_t)
# fun_t2 <- fun_t[,1:24]
# fun_t2_prob <- as.data.frame(t(apply(fun_t2, 1, function(x) x / sum(x))))
# fun_t2_prob$group <- fun_t$group
fun_d <- fun_t%>% filter(group == "D")
#删除group列
fun_d <- fun_d[,-19]
#将fun_d的列随机排列
set.seed(123)
fun_d <- fun_d[,sample(1:18)]
occor = corr.test(fun_d,use="pairwise",method="spearman",adjust="BH",alpha=.05) 
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
library(ggClusterNet)
dat = net_properties.4(mygraph2, n.hub = F)
write.csv(dat, "net_properties.csv")
