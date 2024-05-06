rm(list = ls())
library(tidyverse)
library(readxl)
library(ggpubr)
library(ggvenn)
library(vegan)
library(ggalt)
metadata <- read_excel("metadata_all.xlsx")
bacteria <- read.delim("bacteria.txt")
my_color <- c("#E64B35","#4DBBD5","#00A087","#3C5488",
              "#F39B7F","#4f72d2","#0096b0","#DC0000","#c06e26",
              "#b5426a","#F2A900")
metadata <- metadata %>% filter(parts=="RS" & is.na(intervene))
metadata <- metadata %>% filter(metadata$plant != "Pine")
bacteria <- bacteria %>% select(OTUID,metadata$NO,taxonomy)
bacteria[,1311]
bac <- bacteria
bac[,1311] <- gsub(",c.*","",bac[,1311])
# 若bac最后一列包含p:则只保留p:后面的内容，否则改为Unassign
bac[,1311] <- ifelse(grepl("p:",bac[,1311]),gsub(".*p:","",bac[,1311]),"Unassign")
bac <- bac[,-1]
# 对bac按最后一列用group_by分组求和
bac <- bac %>% group_by(taxonomy) %>% summarise_all(sum)
# 删除bac中除第一列外行和为0的列
bac_fil <- bac[,-1][, colSums(bac[, -1]) != 0]
bac_fil$taxonomy <- bac$taxonomy
bac_fil[,972]
bac_fil[,971]
# 将bac_fil的最后一列提到第一列
bac_fil <- bac_fil[,c(972,1:971)]
bac_fil
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
bac_fil_t <- t(bac_fil)
bac_fil_t <- as.data.frame(bac_fil_t)
bac_fil_t$NO <- rownames(bac_fil_t)
bac_filt <- merge(bac_fil_t,metadata,by="NO")
bac_filt$Others <- rowSums(bac_filt[, 12:69], na.rm = TRUE)
bac_filt <- bac_filt[,c(2:11,82,73)]
result <- bac_filt %>% 
  group_by(plant_status) %>% 
  summarize(across(everything(), mean, na.rm=TRUE))
result[,2:12] = result[,2:12]/rowSums(result[,2:12])*100
# 宽型数据变为长型数据
resultlong <- result %>% gather(key = "Phylum", value = "RelativeAbundance", -plant_status)
p <- ggplot(resultlong, aes(x = plant_status, y = RelativeAbundance, fill = Phylum))+
  geom_col() + theme_bw() + 
  labs(y = "Relative abundance (%)", x="")+
  scale_fill_manual(values = my_color)+
  theme(text = element_text(colour='black',size=15), 
        axis.text=element_text(colour='black',size=13))+
  scale_x_discrete(labels = c("Disease", "Health"))
p
ggsave("tax_phylum.pdf", width = 5.93, height = 7.21, units = "in")
# 按2000抽平
otu_rare <- as.data.frame(t(rrarefy(t(bac_fil),2000)))
# 删除otu_rare中列和小于2000的列
otu_rare <- otu_rare[,colSums(otu_rare)>=2000]
otu_rare <- t(otu_rare)
# 计算加权bray-curtis距离
dune_dist <- vegdist(otu_rare, method="bray", binary=F)
dune_pcoa <- cmdscale(dune_dist, k=3, eig=T)
dune_pcoa_points <- as.data.frame(dune_pcoa$points)
sum_eig <- sum(dune_pcoa$eig)
eig_percent <- round(dune_pcoa$eig/sum_eig*100,1)
colnames(dune_pcoa_points) <- paste0("PCoA", 1:3)
dune_pcoa_points
dune_pcoa_points$NO <- rownames(dune_pcoa_points)
dune_pcoa_result <- merge(dune_pcoa_points,metadata,by="NO")
head(dune_pcoa_result)
ggplot(dune_pcoa_result, aes(x=PCoA1, y=PCoA2, color=plant_status)) +
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep="")) +
  geom_point(size=4
  ) + stat_ellipse(level=0.6) +
  theme_classic()
# 基于bray-curtis距离进行计算
interID=intersect(rownames(otu_rare),metadata$NO)
rownames(metadata) <- metadata$NO
metadata <- metadata[interID,]
set.seed(1)
dune.div <- adonis2(otu_rare ~ plant_status, data = metadata, permutations = 999, method="bray")
dune_adonis <- paste0("adonis: R2","=",round(dune.div$R2,3), "; p=", dune.div$`Pr(>F)`)
p3 <- ggplot(dune_pcoa_result, aes(x=PCoA1, y=PCoA2, color=plant_status, group = plant_status)) +
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep="")) +
  geom_point(size=5, alpha = 0.6) + 
  theme_test() + coord_fixed(1)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  color_palette() +
  scale_color_manual(values =c("#3B4992FF","#008B45FF"),labels=c("Disease", "Health"))+
  scale_fill_manual(values =c("#3B4992FF","#008B45FF"))+
  theme(text = element_text(colour='black',size=15), 
        axis.text=element_text(colour='black',size=13),
        legend.text = element_text(size = 15))+
  theme(legend.title = element_blank()) +
  annotate("text", x = 0.2, y = 0.4, label = dune_adonis[1], size = 5)
  
p3
ggsave("genus_pcoa.pdf", width = 7.21, height = 7.21, units = "in")
