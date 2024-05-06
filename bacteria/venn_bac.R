rm(list = ls())
library(tidyverse)
library(readxl)
library(ggpubr)
library(ggvenn)
library(vegan)
library(ggalt)
library(amplicon)
metadata <- read_excel("metadata_all.xlsx")
bacteria <- read.delim("bacteria.txt")
metadata <- metadata %>% filter(parts=="RS" & is.na(intervene))
metadata <- metadata %>% filter(metadata$plant != "Pine")
bacteria <- bacteria %>% select(OTUID,metadata$NO,taxonomy)
bac2 <- bacteria
# 按属水平分类汇总
bac2[,1311]
bac2[,1311] <- gsub(",s.*","",bac2[,1311])
# 若bac2最后一列包含g:则只保留g:后面的内容，否则改为Unassign
bac2[,1311] <- ifelse(grepl("g:",bac2[,1311]),gsub(".*g:","",bac2[,1311]),"Unassign")
bac2 <- bac2[,-1]
# 对bac按最后一列用group_by分组求和
bac2 <- bac2 %>% group_by(taxonomy) %>% summarise_all(sum)
bac_fil <- bac2[,-1][, colSums(bac2[, -1]) != 0]
bac_fil$taxonomy <- bac2$taxonomy
bac_fil[,972]
bac_fil[,971]
# 将bac_fil的最后一列提到第一列
bac_fil <- bac_fil[,c(972,1:971)]
colSums(bac_fil[,-1])
bac_fil <- as.data.frame(bac_fil)
rownames(bac_fil) <- bac_fil[,1]
bac_fil <- bac_fil[,-1]
bac_fil <- bac_fil[rowSums(bac_fil)!=0,]
intername <- intersect(colnames(bac_fil),metadata$NO)
bac_fil <- bac_fil[,intername]
rownames(metadata) <- metadata$NO
metadata <- metadata[intername,]
bac_fil$sum <- rowSums(bac_fil)
bac_fil <- bac_fil[order(bac_fil$sum,decreasing = TRUE),]
bac_fil <- bac_fil[,-972]
rownames(bac_fil) <- gsub("-", "_", rownames(bac_fil))
rownames(bac_fil) <- gsub("^([0-9])", "A\\1", rownames(bac_fil))
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata[,1]
output = compare(data = bac_fil, metadata = metadata,
                 group = "plant_status", compare_pair = "D-H",
                 method = "wilcox", RA = 0.01,
                 pvalue = 0.05, fdr = 0.1)
output$label = row.names(output)
#筛选阈值确定：p＜0.05，|log2FC|＞1
pvalue = 0.05
log2FC = 1
#根据阈值添加上下调分组标签：
output$group <- case_when(
  output$log2FC > log2FC & output$PValue < pvalue ~ "Enriched_in_Disease",
  output$log2FC < -log2FC & output$PValue < pvalue ~ "Enriched_in_Health",
  TRUE ~ 'Non_significant'
)
result <- output %>% select(group)
result$Genus <- rownames(result)
result <- result %>% filter(Genus != "uncultured")
resultD <- result %>% filter(group %in% c("Enriched_in_Disease", "Non_significant")) %>% .$Genus
resultH <- result %>% filter(group %in% c("Enriched_in_Health", "Non_significant")) %>% .$Genus
p1 <- list(Disease=resultD,Health=resultH) %>% 
  ggvenn(show_percentage = T,show_elements = F,label_sep = ",",
         digits = 1,stroke_color = "black",
         fill_color = c("#3B4992FF", "#008B45FF"),
         set_name_color = c("#3B4992FF","#008B45FF"))
p1
ggsave("genus_venn.pdf", width = 5.93, height = 7.21, units = "in")