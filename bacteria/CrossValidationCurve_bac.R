rm(list = ls())
library(tidyverse)
library(readxl)
library(vegan)
library(randomForest)
library(ROCR)
library("caret")
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
# 不同研究的发现率要大于1/2
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
set.seed(123)
bac_t.rf <- randomForest(group ~ ., data = bac_t, importance = TRUE, proximity = TRUE)
print(bac_t.rf)
write.table(bac_t.rf$confusion, file = "b_genus_confusion.txt", sep = "\t", quote = F, row.names = T, col.names = T)
imp = as.data.frame(round(importance(bac_t.rf), 2))
imp=imp[order(imp$MeanDecreaseAccuracy,decreasing = F),]
write.table(imp, file = "b_genus_imp.txt", sep = "\t", quote = F, row.names = T, col.names = T)
n = ncol(bac_t)-1
mybac_t = bac_t[1:n]
set.seed(123)
result = rfcv(mybac_t, bac_t$group, cv.fold = 5, scale = "log", step = 0.9)
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))
result1 = result
error.cv = data.frame(num = result$n.var, error.1 =  result$error.cv)
# 用另外4组随机数来看错误率，然后绘图
for (i in 124:127){
  print(i)
  set.seed(i)
  result= rfcv(mybac_t, bac_t$group, cv.fold=5, scale = "log", step = 0.9)
  error.cv = cbind(error.cv, result$error.cv)
}
## 绘制交叉验证曲线
n.var = error.cv$num
error.cv = error.cv[,2:6]
colnames(error.cv) = paste('err',1:5,sep='.')
err.mean = apply(error.cv,1,mean)
allerr = data.frame(num=n.var,err.mean=err.mean,error.cv)
# number of features selected
optimal = 24

write.table(allerr, file = "b_genus_crosstest.txt", sep = "\t", quote = F, row.names = T, col.names = T)
p = ggplot() + 
  geom_line(aes(x = allerr$num, y = allerr$err.1), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.2), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.3), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.4), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.mean), colour = 'black') + 
  geom_vline(xintercept = optimal, colour='black', lwd=0.36, linetype="dashed") + 
  coord_trans(x = "log2") +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 20, 30, 50, 100, 200)) + # , max(allerr$num)
  labs(title=paste('Training set (n = ', dim(bac_t)[1],')', sep = ''), 
       x='Number of Genus', y='Cross-validation error rate') + 
  annotate("text", x = optimal, y = max(allerr$err.mean), label=paste("optimal = ", optimal, sep="")) + theme_test()
p
