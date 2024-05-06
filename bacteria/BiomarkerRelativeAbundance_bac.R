rm(list = ls())
library(tidyverse)
library(readxl)
library(vegan)
library(randomForest)
library(ROCR)
library(caret)
library(ggpubr)
metadata <- read_excel("metadata_all.xlsx")
bacteria <- read.delim("bacteria.txt")
metadata <- metadata %>% filter(parts=="RS" & is.na(intervene))
metadata <- metadata %>% filter(plant != "Pine")
set.seed(123)
metadata <- metadata %>% group_by(paper) %>% sample_n(if(n()>50) 50 else n())
bac <- bacteria
id <- intersect(colnames(bac),metadata$NO)
bac <- bac[,c("OTUID",id,"taxonomy")]

rownames(bac) <- bac[,1]
bac <- bac[,-1]
bac <- bac[,-710]
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
# 将大于1的值改为1
bac[bac > 1] <- 1
# 不同样品的发现率要大于1/3
bac2 <- bac %>% filter(rowSums(bac)/ncol(bac) >= 1/3)
bac3 <- t(bac2)
bac3 <- as.data.frame(bac3)
bac3$NO <- rownames(bac3)
meta3 <- metadata[,c(1,12)]
bac3 <- merge(bac3, meta3, by = "NO")
rownames(bac3) <- bac3[,1]
bac3 <- bac3[,-1]
summary <- bac3 %>%
  group_by(paper) %>%
  summarize_all(sum)
summary <- as.data.frame(summary)
rownames(summary) <- summary[,1]
summary <- summary[,-1]
bac4 <- t(summary)
bac4 <- as.data.frame(bac4)
bac4[bac4>1] <- 1
# 不同研究的发现率要大于50%
bac5 <- bac4 %>% filter(rowSums(bac4)/ncol(bac4) >= 1/2)
otuname <- rownames(bac5)
bacteria <- bac1 %>% filter(OTUID %in% otuname)
colSums(bacteria[,c(-1,-402)])
id <- intersect(colnames(bacteria),metadata$NO)
# filter选取metadata中NO为id的行
metadata <- metadata %>% filter(NO %in% id)
bacteria[,402] <- gsub(",s.*","",bacteria[,402])
# 若bac最后一列包含g:则只保留g:后面的内容，否则改为Unassign
bacteria[,402] <- ifelse(grepl("g:",bacteria[,402]),gsub(".*g:","",bacteria[,402]),"Unassign")
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
#filter选取Phylum不为""的行
imp <- imp %>% filter(Phylum != "")
id2 <- c(imp$Genus,"group")
# 选取bac_t中列名为id2的列
bac_t <- bac_t %>% select(id2)

#差异比较
# 每个值除以该行和
bac_t <- as.data.frame(bac_t)
bac_t2 <- bac_t[,1:24]
bac_t2_prob <- as.data.frame(t(apply(bac_t2, 1, function(x) x / sum(x))))
bac_t2_prob$group <- bac_t$group
# 将每一列按group组做wilcox.test检验
wilcox <- apply(bac_t2_prob[,1:24], 2, function(x) wilcox.test(x ~ bac_t2_prob$group))
class(wilcox$Gemmatimonas$p.value)
# 提取wilcox中的p.value
pvalue <- sapply(wilcox, function(x) x$p.value)
pvalue <- as.data.frame(pvalue)
bac_t3 <- bac_t2_prob %>% pivot_longer(cols = 1:24, names_to = "Genus", values_to = "RA")
# 删除RA中的无效值
bac_t3 <- bac_t3 %>% filter(!is.na(RA))
# group_by按group和Genus分组，将RA变为均值
bac_t4 <- bac_t3 %>% group_by(group,Genus) %>% summarize(meanRA = mean(RA))
# 更改因子顺序
bac_t4$Genus <- factor(bac_t4$Genus,levels = rownames(pvalue))
bac_t4$meanRA <- bac_t4$meanRA * 100
# pvalue中pvalue的值如果大于0.01小于0.05则改为*，大于0.001小于0.01则改为**，大于0小于0.001则改为***，否则改为""
pvalue <- pvalue %>% mutate(pvalue = ifelse(pvalue > 0.01 & pvalue <= 0.05,"*",
                                            ifelse(pvalue > 0.001 & pvalue <= 0.01,"**",
                                                   ifelse(pvalue > 0 & pvalue <= 0.001,"***",""))))
mycol = c("#3B4992FF","#008B45FF")
bac_t4$meanRA <- -bac_t4$meanRA
p <- ggplot(bac_t4, aes(x = Genus, y = meanRA, fill = group))+
  geom_col() + theme_classic()+scale_x_discrete(position = "top")+
  scale_y_continuous(expand=c(0, 0),limits = c(-35,0))+
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

