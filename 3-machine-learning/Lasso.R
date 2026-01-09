#LASSO回归,输入只能是矩阵,行为样本,列为基因,分组只能是(0,1)
library(glmnet)

input <- read.csv("175sig_genes_norm_exp.csv",header = T,row.names = 1)

x <- input[,-1]
# levels(groups) <- c(0, 1)
y <- as.integer(as.character(input$group))
fit <- glmnet(x, y, family = "binomial", nlambda = 100, alpha = 1)
plot(fit, xvar = "lambda", label = TRUE)
cvfit <- cv.glmnet(data.matrix(x), y, nfolds = 10)
cvfit$lambda.min
cvfit$lambda.1se

pdf("Binomial Deviance lasso.pdf", height = 8, width = 10)
plot(cvfit, ylab = "Binomial Deviance",cex.lab = 1.25)
dev.off()
coef <- coef(cvfit, s = "lambda.min")
factors <- as.matrix(coef)
factors[factors == 0] <- NA
factors <- na.omit(factors)
write.csv(factors, "lasso.csv")
