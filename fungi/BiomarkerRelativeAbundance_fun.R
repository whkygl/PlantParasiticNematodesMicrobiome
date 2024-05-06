rm(list = ls())
library(tidyverse)
library(readxl)
library(vegan)
library(randomForest)
library(ROCR)
library(caret)
library(ggpubr)
metadata <- read_excel("metadata_all.xlsx")
fungi <- read.delim("fungi.txt")
metadata <- metadata %>% filter(parts=="RS" & is.na(intervene))
metadata <- metadata %>% filter(plant != "Pine")
set.seed(123)
metadata <- metadata %>% group_by(paper) %>% sample_n(if(n()>50) 50 else n())
fun <- fungi
id <- intersect(colnames(fun),metadata$NO)
fun <- fun[,c("OTUID",id,"taxonomy")]

rownames(fun) <- fun[,1]
fun <- fun[,-1]
fun <- fun[,-710]
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
# 将大于1的值改为1
fun[fun > 1] <- 1
# 不同样品的发现率要大于1/3
fun2 <- fun %>% filter(rowSums(fun)/ncol(fun) >= 1/3)
fun3 <- t(fun2)
fun3 <- as.data.frame(fun3)
fun3$NO <- rownames(fun3)
meta3 <- metadata[,c(1,12)]
fun3 <- merge(fun3, meta3, by = "NO")
rownames(fun3) <- fun3[,1]
fun3 <- fun3[,-1]
summary <- fun3 %>%
  group_by(paper) %>%
  summarize_all(sum)
summary <- as.data.frame(summary)
rownames(summary) <- summary[,1]
summary <- summary[,-1]
fun4 <- t(summary)
fun4 <- as.data.frame(fun4)
fun4[fun4>1] <- 1
# 不同研究的发现率要大于1/2
fun5 <- fun4 %>% filter(rowSums(fun4)/ncol(fun4) >= 1/2)
otuname <- rownames(fun5)
fungi <- fun1 %>% filter(OTUID %in% otuname)
colSums(fungi[,c(-1,-299)])
id <- intersect(colnames(fungi),metadata$NO)
# filter选取metadata中NO为id的行
metadata <- metadata %>% filter(NO %in% id)
fungi[,299] <- gsub(",s.*","",fungi[,299])
# 若fun最后一列包含g:则只保留g:后面的内容，否则改为Unassign
fungi[,299] <- ifelse(grepl("g:",fungi[,299]),gsub(".*g:","",fungi[,299]),"Unassign")
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
#filter选取Phylum不为""的行
imp <- imp %>% filter(Phylum != "")
id2 <- c(imp$Genus,"group")
# 选取fun_t中列名为id2的列
fun_t <- fun_t %>% select(id2)

#差异比较
# 每个值除以该行和
fun_t <- as.data.frame(fun_t)
fun_t2 <- fun_t[,1:18]
fun_t2_prob <- as.data.frame(t(apply(fun_t2, 1, function(x) x / sum(x))))
fun_t2_prob$group <- fun_t$group
# 将每一列按group组做wilcox.test检验
wilcox <- apply(fun_t2_prob[,1:18], 2, function(x) wilcox.test(x ~ fun_t2_prob$group))
# 提取wilcox中的p.value
pvalue <- sapply(wilcox, function(x) x$p.value)
pvalue <- as.data.frame(pvalue)
fun_t3 <- fun_t2_prob %>% pivot_longer(cols = 1:18, names_to = "Genus", values_to = "RA")
# 删除RA中的无效值
fun_t3 <- fun_t3 %>% filter(!is.na(RA))
# group_by按group和Genus分组，将RA变为均值
fun_t4 <- fun_t3 %>% group_by(group,Genus) %>% summarize(meanRA = mean(RA))
# 更改因子顺序
fun_t4$Genus <- factor(fun_t4$Genus,levels = rownames(pvalue))
fun_t4$meanRA <- fun_t4$meanRA * 100
# pvalue中pvalue的值如果大于0.01小于0.05则改为*，大于0.001小于0.01则改为**，大于0小于0.001则改为***，否则改为""
pvalue <- pvalue %>% mutate(pvalue = ifelse(pvalue > 0.01 & pvalue <= 0.05,"*",
                                            ifelse(pvalue > 0.001 & pvalue <= 0.01,"**",
                                                   ifelse(pvalue > 0 & pvalue <= 0.001,"***",""))))
mycol = c("#3B4992FF","#008B45FF")
fun_t4$meanRA <- -fun_t4$meanRA
p <- ggplot(fun_t4, aes(x = Genus, y = meanRA, fill = group))+
  geom_col() + theme_classic()+scale_x_discrete(position = "top")+
  scale_y_continuous(expand=c(0, 0),limits = c(-45,0))+
  scale_fill_manual(values = mycol)+
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.8,0.2)) +
  #去掉图例标题
  guides(fill = guide_legend(title = NULL))+
  #修改y轴标题为“Mean decrease in accuracy”
  ylab("Relative abundance (%)")

p