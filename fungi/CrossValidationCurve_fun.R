rm(list = ls())
library(tidyverse)
library(readxl)
library(vegan)
library(randomForest)
library(ROCR)
library("caret")
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
# 不同研究的发现率要大于50%
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
set.seed(123)
fun_t.rf <- randomForest(group ~ ., data = fun_t, importance = TRUE, proximity = TRUE)
print(fun_t.rf)
write.table(fun_t.rf$confusion, file = "b_genus_confusion.txt", sep = "\t", quote = F, row.names = T, col.names = T)
imp = as.data.frame(round(importance(fun_t.rf), 2))
imp=imp[order(imp$MeanDecreaseAccuracy,decreasing = F),]
write.table(imp, file = "b_genus_imp.txt", sep = "\t", quote = F, row.names = T, col.names = T)
n = ncol(fun_t)-1
myfun_t = fun_t[1:n]
set.seed(123)
result = rfcv(myfun_t, fun_t$group, cv.fold = 5, scale = "log", step = 0.9)
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))
result1 = result
error.cv = data.frame(num = result$n.var, error.1 =  result$error.cv)
# 用另外4组随机数来看错误率，然后绘图
for (i in 124:127){
  print(i)
  set.seed(i)
  result= rfcv(myfun_t, fun_t$group, cv.fold=5, scale = "log", step = 0.9)
  error.cv = cbind(error.cv, result$error.cv)
}
## 绘制交叉验证曲线
n.var = error.cv$num
error.cv = error.cv[,2:6]
colnames(error.cv) = paste('err',1:5,sep='.')
err.mean = apply(error.cv,1,mean)
allerr = data.frame(num=n.var,err.mean=err.mean,error.cv)
# number of features selected
optimal = 18

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
  labs(title=paste('Training set (n = ', dim(fun_t)[1],')', sep = ''), 
       x='Number of Genus', y='Cross-validation error rate') + 
  annotate("text", x = optimal, y = max(allerr$err.mean), label=paste("optimal = ", optimal, sep="")) + theme_test()
p
