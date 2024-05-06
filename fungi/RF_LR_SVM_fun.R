rm(list = ls())
library(randomForest)
library(e1071)
library(tidyverse)
library(caret)
library(yardstick)
library(vegan)
fungi <- read.delim("fungi.txt")
# 自定义不同算法的五折交叉验证函数
cv <- function(data, class, method, k = 5, ...){
  folds = createFolds(data[[class]], k = k)
  result = data.frame(matrix(ncol = 5, nrow = nrow(data)))
  colnames(result) = c("times", "true", "pred", "method", "prob")
  formula1 <- as.formula(paste(class, "~ ."))
  for(i in 1:k){
    train = data[-folds[[i]],]
    test = data[folds[[i]],]
    model = method(formula = formula1, data = train, ...)
    if (class(model)[1] == "randomForest.formula") {
      pred = predict(model, type = 'prob', newdata = test)
      pred2 = predict(model, type = 'class', newdata = test)
      result[folds[[i]], "times"] = i
      result[folds[[i]], "true"] = as.character(test[[class]])
      result[folds[[i]], "pred"] = as.character(pred2)
      result[folds[[i]], "method"] = class(model)[1]
      result[folds[[i]], "prob"] = pred[,1]
    } else if (class(model)[1] == "svm.formula") {
      pred = predict(model, probability = TRUE, newdata = test)
      pred = attr(pred, "probabilities")[,1]
      pred2 = predict(model, newdata = test)
      result[folds[[i]], "times"] = i
      result[folds[[i]], "true"] = as.character(test[[class]])
      result[folds[[i]], "pred"] = as.character(pred2)
      result[folds[[i]], "method"] = class(model)[1]
      result[folds[[i]], "prob"] = pred
    } else if (class(model)[1] == "glm") {
      pred = predict(model, type = 'response', newdata = test)
      pred = 1-pred
      pred2 = pred>0.5
      result[folds[[i]], "times"] = i
      result[folds[[i]], "true"] = as.character(test[[class]])
      result[folds[[i]], "pred"] = as.character(pred2)
      result[folds[[i]], "method"] = class(model)[1]
      result[folds[[i]], "prob"] = pred
    }
  }
  return(result)
}
RF <- cv(fungi, 'Genus', randomForest, ntree = 1000)
SVM <- cv(fungi, 'Genus', probability = TRUE, svm, kernel = 'radial')
LR <- cv(fungi, 'Genus', glm, family = binomial)
f_all_genus_method <- rbind(RF, SVM, LR)
f_all_genus_method$true <- as.factor(f_all_genus_method$true)
# 绘制ROC曲线
g1 <- f_all_genus_method %>% 
  group_by(method) %>% 
  roc_curve(true, prob) %>%
  autoplot() + theme_test() + 
  guides(color = guide_legend(title = NULL))+
  scale_color_manual(values=c('#5494cc','#0d898a','#e18283')) +
  theme(legend.position=c(.85,.2))
g1
# 计算AUC值
auc1 <- f_all_genus_method %>% 
  group_by(method) %>% 
  roc_auc(true, prob)