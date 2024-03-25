library(tidyverse)
library(readxl)
library(vegan)
library(ggplot2)
library(ggpubr)
pacman::p_load(tidyverse,ggtext,camcorder,scales,ggsci,ggdist,gghalves)
metadata <- read_excel("metadata_all.xlsx")
metadata <- metadata %>% filter(metadata$plant != "Pine")
bacteria <- read.delim("bacteria.txt")
bacteria <- bacteria %>% select(OTUID,metadata$NO,taxonomy)
bacteria[,2329]
bac <- bacteria
bac[,2329] <- gsub(",s.*","",bac[,2329])
bac[,2329] <- ifelse(grepl("g:",bac[,2329]),gsub(".*g:","",bac[,2329]),"Unassign")
bac <- bac[,-1]
# 对bac按最后一列用group_by分组求和
bac <- bac %>% group_by(taxonomy) %>% summarise_all(sum)
# 删除bac中除第一列外行和为0的列
bac_fil <- bac[,-1][, colSums(bac[, -1]) != 0]
bac_fil$taxonomy <- bac$taxonomy
bac_fil[,1670]
bac_fil[,1669]
# 将bac_fil的最后一列提到第一列
bac_fil <- bac_fil[,c(1670,1:1669)]
bac_fil
colSums(bac_fil[,-1])
bac_fil <- as.data.frame(bac_fil)
rownames(bac_fil) <- bac_fil[,1]
bac_fil <- bac_fil[,-1]
# 按2000抽平
otu_rare <- as.data.frame(t(rrarefy(t(bac_fil),2000)))
# 删除otu_rare中列和小于2000的列
otu_rare <- otu_rare[,colSums(otu_rare)>=2000]
# 求otu_rare的alpha多样性
otu_rare_t <- t(otu_rare)
shannon <- diversity(otu_rare_t,"shannon")
simpson <- diversity(otu_rare_t,"simpson")
invsimpson <- diversity(otu_rare_t,"inv")
richness=specnumber(otu_rare_t)
pielou=shannon/log(richness)
?estimateR
hehe <- estimateR(otu_rare_t)
hehe <- as.data.frame(t(hehe))
hehe$shannon <- shannon
hehe$simpson <- simpson
hehe$invsimpson <- invsimpson
hehe$richness <- richness
hehe$pielou <- pielou
hehe$NO <- rownames(hehe)
alpha <- left_join(hehe,metadata,by="NO")
alpha2 <- alpha %>% filter(parts=="RS" & is.na(intervene))
alpha2$microbiome
###定义绘图函数
"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x,
                     xmax = x + width / 2)
            
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                             plyr::arrange(transform(data, x = xmaxv), -y))
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1,])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),
          
          required_aes = c("x", "y")
  )

my_comparisons <- list(c("D", "H"))
p1 <- ggplot(alpha2, aes(x=plant_status, y=richness)) +
  scale_color_manual(values =c("#3B4992FF","#008B45FF"))+
  scale_fill_manual(values =c("#3B4992FF","#008B45FF"))+
  scale_x_discrete(labels = c("Disease", "Health"))+
  geom_half_violin(aes(fill=plant_status,color=plant_status),position=position_nudge(x=.2),width=0.6) +
  geom_jitter(aes(color=plant_status), width=0.1,alpha=0.5) +
  geom_boxplot(width=.07,position=position_nudge(x=0.2),fill="white",size=0.2) +
  theme_bw()+theme(panel.grid = element_blank())+
  theme( axis.text = element_text(size=12,colour = "black"),
         axis.title =  element_text(size=15,colour = "black"),
         axis.text.x =element_text(size = 13),
         legend.position="none")+
  labs(x =NULL,y = 'Richness')+
  stat_compare_means(comparisons = my_comparisons,
                     #label = "p.signif",
                     label = "p.format",
                     method = "wilcox.test"
  )
p1
p2 <- ggplot(alpha2, aes(x=plant_status, y=shannon)) +
  scale_color_manual(values =c("#3B4992FF","#008B45FF"))+
  scale_fill_manual(values =c("#3B4992FF","#008B45FF"))+
  scale_x_discrete(labels = c("Disease", "Health"))+
  geom_flat_violin(aes(fill=plant_status,color=plant_status),position=position_nudge(x=.2),width=0.6) +
  geom_jitter(aes(color=plant_status), width=0.1,alpha=0.5) +
  geom_boxplot(width=.07,position=position_nudge(x=0.2),fill="white",size=0.2) +
  theme_bw()+theme(panel.grid = element_blank())+
  theme( axis.text = element_text(size=12,colour = "black"),
         axis.title =  element_text(size=15,colour = "black"),
         axis.text.x =element_text(size = 13),
         legend.position="none")+
  labs(x =NULL,y = 'Shannon')+
  stat_compare_means(comparisons = my_comparisons,
                     #label = "p.signif",
                     label = "p.format",
                     method = "wilcox.test"
  )
p2
p3 <- ggplot(alpha2, aes(x=plant_status, y=S.chao1)) +
  scale_color_manual(values =c("#3B4992FF","#008B45FF"))+
  scale_fill_manual(values =c("#3B4992FF","#008B45FF"))+
  scale_x_discrete(labels = c("Disease", "Health"))+
  geom_flat_violin(aes(fill=plant_status,color=plant_status),position=position_nudge(x=.2),width=0.6) +
  geom_jitter(aes(color=plant_status), width=0.1,alpha=0.5) +
  geom_boxplot(width=.07,position=position_nudge(x=0.2),fill="white",size=0.2) +
  theme_bw()+theme(panel.grid = element_blank())+
  theme( axis.text = element_text(size=12,colour = "black"),
         axis.title =  element_text(size=15,colour = "black"),
         axis.text.x =element_text(size = 13),
         legend.position="none")+
  labs(x =NULL,y = 'Chao1')+
  stat_compare_means(comparisons = my_comparisons,
                     #label = "p.signif",
                     label = "p.format",
                     method = "wilcox.test"
  )
p3
p4 <- ggplot(alpha2, aes(x=plant_status, y=S.ACE)) +
  scale_color_manual(values =c("#3B4992FF","#008B45FF"))+
  scale_fill_manual(values =c("#3B4992FF","#008B45FF"))+
  scale_x_discrete(labels = c("Disease", "Health"))+
  geom_flat_violin(aes(fill=plant_status,color=plant_status),position=position_nudge(x=.2),width=0.6) +
  geom_jitter(aes(color=plant_status), width=0.1,alpha=0.5) +
  geom_boxplot(width=.07,position=position_nudge(x=0.2),fill="white",size=0.2) +
  theme_bw()+theme(panel.grid = element_blank())+
  theme( axis.text = element_text(size=12,colour = "black"),
         axis.title =  element_text(size=15,colour = "black"),
         axis.text.x =element_text(size = 13),
         legend.position="none")+
  labs(x =NULL,y = 'ACE')+
  stat_compare_means(comparisons = my_comparisons,
                     #label = "p.signif",
                     label = "p.format",
                     method = "wilcox.test"
  )
p4
p5 <- ggplot(alpha2, aes(x=plant_status, y=pielou)) +
  scale_color_manual(values =c("#3B4992FF","#008B45FF"))+
  scale_fill_manual(values =c("#3B4992FF","#008B45FF"))+
  scale_x_discrete(labels = c("Disease", "Health"))+
  geom_flat_violin(aes(fill=plant_status,color=plant_status),position=position_nudge(x=.2),width=0.6) +
  geom_jitter(aes(color=plant_status), width=0.1,alpha=0.5) +
  geom_boxplot(width=.07,position=position_nudge(x=0.2),fill="white",size=0.2) +
  theme_bw()+theme(panel.grid = element_blank())+
  theme( axis.text = element_text(size=12,colour = "black"),
         axis.title =  element_text(size=15,colour = "black"),
         axis.text.x =element_text(size = 13),
         legend.position="none")+
  labs(x =NULL,y = 'Pielou evenness')+
  stat_compare_means(comparisons = my_comparisons,
                     #label = "p.signif",
                     label = "p.format",
                     method = "wilcox.test"
  )
p5
p6 <- ggplot(alpha2, aes(x=plant_status, y=simpson)) +
  scale_color_manual(values =c("#3B4992FF","#008B45FF"))+
  scale_fill_manual(values =c("#3B4992FF","#008B45FF"))+
  scale_x_discrete(labels = c("Disease", "Health"))+
  geom_flat_violin(aes(fill=plant_status,color=plant_status),position=position_nudge(x=.2),width=0.6) +
  geom_jitter(aes(color=plant_status), width=0.1,alpha=0.5) +
  geom_boxplot(width=.07,position=position_nudge(x=0.2),fill="white",size=0.2) +
  theme_bw()+theme(panel.grid = element_blank())+
  theme( axis.text = element_text(size=12,colour = "black"),
         axis.title =  element_text(size=15,colour = "black"),
         axis.text.x =element_text(size = 13),
         legend.position="none")+
  labs(x =NULL,y = 'Simpson')+
  stat_compare_means(comparisons = my_comparisons,
                     #label = "p.signif",
                     label = "p.format",
                     method = "wilcox.test"
  )
p6
p7 <- ggplot(alpha2, aes(x=plant_status, y=invsimpson)) +
  scale_color_manual(values =c("#3B4992FF","#008B45FF"))+
  scale_fill_manual(values =c("#3B4992FF","#008B45FF"))+
  scale_x_discrete(labels = c("Disease", "Health"))+
  geom_flat_violin(aes(fill=plant_status,color=plant_status),position=position_nudge(x=.2),width=0.6) +
  geom_jitter(aes(color=plant_status), width=0.1,alpha=0.5) +
  geom_boxplot(width=.07,position=position_nudge(x=0.2),fill="white",size=0.2) +
  theme_bw()+theme(panel.grid = element_blank())+
  theme( axis.text = element_text(size=12,colour = "black"),
         axis.title =  element_text(size=15,colour = "black"),
         axis.text.x =element_text(size = 13),
         legend.position="none")+
  labs(x =NULL,y = 'InvSimpson')+
  stat_compare_means(comparisons = my_comparisons,
                     #label = "p.signif",
                     label = "p.format",
                     method = "wilcox.test"
  )
p7
p <- ggarrange(p1,p2,p3,p4,p5,p6,p7,ncol = 4,nrow = 2,labels = c("A","B","C","D","E","F","G"))
p
g <- ggarrange(p1,p4,p6,p7,ncol = 4)
g
alpha3 <- alpha2 %>% select(plant_status,S.chao1,shannon,pielou)
alpha4 <- alpha3 %>% pivot_longer(cols = c(S.chao1,shannon,pielou),
                                  names_to = "group",
                                  values_to = "value")
alpha4$group <- as.factor(alpha4$group)
alpha4 <- alpha4 %>% mutate(group = recode(group,"pielou" = "Pielou_evenness","S.chao1"="Chao1","shannon"="Shannon"))
g2 <- ggplot(alpha4,aes(x=plant_status, y = value))+
  scale_color_manual(values =c("#3B4992FF","#008B45FF"))+
  scale_fill_manual(values =c("#3B4992FF","#008B45FF"))+
  scale_x_discrete(labels = c("Disease", "Health"))+
  geom_half_violin(aes(fill=plant_status,color=plant_status),position=position_nudge(x=.2),width=0.6,side = "r") +
  geom_jitter(aes(color=plant_status), width=0.1,alpha=0.5) +
  geom_boxplot(width=.07,position=position_nudge(x=0.2),fill="white",size=0.2) +
  theme_bw()+theme(panel.grid = element_blank())+
  theme( axis.text = element_text(size=12,colour = "black"),
         axis.title =  element_text(size=15,colour = "black"),
         axis.text.x =element_text(size = 13),
         legend.position="none",
         strip.text = element_text(size = rel(1.05)))+
  facet_wrap(~group,scales = "free")+
  stat_compare_means(comparisons = my_comparisons,
                     #label = "p.signif",
                     label = "p.format",
                     method = "wilcox.test"
  )+xlab(NULL)+ylab(NULL)
g2
