rm(list = ls())
library(tidyverse)
library(readxl)
library(ggpubr)
library(ggvenn)
library(vegan)
library(ggalt)
library(amplicon)
metadata <- read_excel("metadata_all.xlsx")
fungi <- read.delim("fungi.txt")
metadata <- metadata %>% filter(parts=="RS" & is.na(intervene))
metadata <- metadata %>% filter(metadata$plant != "Pine")
fungi <- fungi %>% select(OTUID,metadata$NO,taxonomy)
fun2 <- fungi
# 按属水平分类汇总
fun2[,1311]
fun2[,1311] <- gsub(",s.*","",fun2[,1311])
# 若fun2最后一列包含g:则只保留g:后面的内容，否则改为Unassign
fun2[,1311] <- ifelse(grepl("g:",fun2[,1311]),gsub(".*g:","",fun2[,1311]),"Unassign")
fun2 <- fun2[,-1]
# 对fun按最后一列用group_by分组求和
fun2 <- fun2 %>% group_by(taxonomy) %>% summarise_all(sum)
fun_fil <- fun2[,-1][, colSums(fun2[, -1]) != 0]
fun_fil$taxonomy <- fun2$taxonomy
fun_fil[,950]
fun_fil[,949]
# 将fun_fil的最后一列提到第一列
fun_fil <- fun_fil[,c(950,1:949)]
colSums(fun_fil[,-1])
fun_fil <- as.data.frame(fun_fil)
rownames(fun_fil) <- fun_fil[,1]
fun_fil <- fun_fil[,-1]
fun_fil <- fun_fil[rowSums(fun_fil)!=0,]
intername <- intersect(colnames(fun_fil),metadata$NO)
fun_fil <- fun_fil[,intername]
rownames(metadata) <- metadata$NO
metadata <- metadata[intername,]
fun_fil$sum <- rowSums(fun_fil)
fun_fil <- fun_fil[order(fun_fil$sum,decreasing = TRUE),]
fun_fil <- fun_fil[,-950]
rownames(fun_fil) <- gsub("-", "_", rownames(fun_fil))
rownames(fun_fil) <- gsub("^([0-9])", "A\\1", rownames(fun_fil))
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata[,1]
output = compare(data = fun_fil, metadata = metadata,
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
