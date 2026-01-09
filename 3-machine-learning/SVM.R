###  Support Vector Machine 基于支持向量机的分类预测(SVM)  ###
library(e1071)

input <- read.csv("./svm/175sig_genes_norm_exp.csv",header = T,row.names = 1)

set.seed(2023)
source("./svm/msvmRFE.R")

nfold = 10
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=15, halve.above=200)
top.features = WriteFeatures(results, input, save=F)
# 这里选择前10的特征实行计算，也可以扩大范围，按照各自情况来选择。
featsweep = lapply(1:10, FeatSweep.wrap, results, input)
# 绘制泛化误差与特征数量的图
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
pdf("svm_feature.pdf",width=10, height=8, bg='white')
PlotErrors(errors, no.info=no.info)
# 添加网格线和最佳特征数标记
grid()
opt_features <- which.min(errors)
points(opt_features, errors[opt_features], col = "red", pch = 19, cex = 1.5)
text(opt_features, errors[opt_features], paste("Opt:", opt_features), pos = 4, col = "red")
dev.off()

# 输出最佳特征数和对应的特征
cat("最佳特征数量:", opt_features, "\n")
if(opt_features > 0) {
  cat("最佳特征列表:\n")
  print(head(top.features, opt_features))
}

write.csv(top.features,"top.features.csv",row.names = F)

