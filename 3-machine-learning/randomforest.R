# 随机森林，行为样本，列为基因，搞一列为分类 
library(randomForest)
library(ggplot2)
library(RColorBrewer)

input <- read.csv("../lasso/175sig_genes_norm_exp.csv",header = T,row.names = 1)
# 数据预处理
infr <- as.data.frame(input)
colnames(infr) <- make.names(colnames(infr))
infr$group <- as.factor(infr[, 1])  # 统一转换为因子，避免重复操作

# 1. 优化mtry参数筛选
results <- data.frame(m = 1:175, error = numeric(175))
for(i in 1:175) {
  model <- randomForest(group ~ ., data = infr, mtry = i, ntree = 1000)
  results$error[i] <- model$err.rate[nrow(model$err.rate), 1]
}
results1 <- results[1:150, ]  # 截取前150个mtry值

# 绘制mtry-误差图
p1 <- ggplot(na.omit(results1), aes(x = m, y = error)) + 
  geom_line(linewidth = 0.7) + 
  geom_point(data = results1[which.min(results1$error), , drop = FALSE], 
             aes(x = m, y = error), 
             shape = 2, size = 4, color = "red", stroke = 0.5) + 
  labs(x = "Number of predictors", y = "OOB Error") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.9))
print(p1)
ggsave("mtry_error_plot.pdf", p1, width = 8, height = 6, dpi = 300)

# 2. 树数量收敛性分析（使用最优mtry的模型）
best_mtry <- results1$m[which.min(results1$error)]
model_converge <- randomForest(group ~ ., data = infr, mtry = best_mtry, ntree = 1000)

# 整理OOB误差数据
oob_error_data <- data.frame(
  Trees = rep(1:nrow(model_converge$err.rate), 3),
  Type = rep(c("OOB", "0", "1"), each = nrow(model_converge$err.rate)),
  Error = c(model_converge$err.rate[, "OOB"],
            model_converge$err.rate[, "0"],
            model_converge$err.rate[, "1"])
)

# 绘制树数量-误差图
p2 <- ggplot(oob_error_data, aes(x = Trees, y = Error, color = Type, linetype = Type)) +
  geom_line(linewidth = 1) +
  scale_color_brewer(palette = "Set2") +
  labs(title = "Error Rate by Number of Trees", x = "Number of Trees", y = "Error Rate") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "bottom")
print(p2)
ggsave("tree_convergence_plot.pdf", p2, width = 10, height = 6, dpi = 300)

# 3. 变量重要性分析（使用最终模型）
final_model <- randomForest(group ~ ., data = infr, mtry = best_mtry, ntree = 1000, importance = TRUE)
importance_df <- as.data.frame(importance(final_model))
importance_df$Variable <- rownames(importance_df)
importance_df <- importance_df[order(importance_df$MeanDecreaseGini, decreasing = TRUE), ]

# 提取Top25特征
top25_df <- importance_df[1:25, ]
top_features <- top25_df$Variable

# 绘制Top25变量重要性图
p3 <- ggplot(top25_df, aes(x = reorder(Variable, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", aes(fill = Variable)) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set3"))(25)) +
  coord_flip() +
  labs(title = "Top 25 Variable Importance Plot", x = "Gene", y = "MeanDecreaseGini") +
  theme_minimal() +
  theme(
    panel.background = element_rect(color = 'black', fill = 'transparent'),
    axis.text.x = element_text(color = 'black', size = 14),
    axis.text.y = element_text(color = 'black', size = 14),
    axis.title.x = element_text(color = 'black', size = 16),
    axis.title.y = element_text(color = 'black', size = 16),
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "none"
  )
print(p3)
ggsave("top25_randf.pdf", p3, height = 8, width = 10, dpi = 300)

# 4. 特征数量优化分析
oob_errors <- numeric(length(top_features))
for (i in 1:length(top_features)) {
  model <- randomForest(
    group ~ ., 
    data = infr[, c("group", top_features[1:i])],  # 统一使用infr数据框
    ntree = 1000,
    importance = TRUE
  )
  oob_errors[i] <- model$err.rate[nrow(model$err.rate), "OOB"]
}

# 绘制特征数量-OOB误差图
obb_df <- data.frame(feature_count = 1:length(top_features), oob_error = oob_errors)
opt_features <- which.min(obb_df$oob_error)

p4 <- ggplot(obb_df, aes(x = feature_count, y = oob_error)) +
  geom_line(color = "#2c7fb8", linewidth = 1.2) +
  geom_point(color = "#2c7fb8", size = 3) +
  geom_vline(xintercept = opt_features, color = "#e41a1c", linetype = "dashed", linewidth = 1) +
  annotate("text", x = opt_features, y = max(obb_df$oob_error),
           label = paste("Opt =", opt_features), color = "#e41a1c", hjust = -0.1, size = 6) +
  labs(title = "OOB Error vs Feature Count", x = "Number of Features", y = "OOB Error Rate") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text = element_text(color = "black", size = 8),
    axis.title = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.line = element_line(color = "black", linewidth = 0.4)
  ) +
  scale_x_continuous(breaks = seq(0, length(top_features), by = 5))
print(p4)
ggsave("obb_feature_plot.pdf", p4, width = 6, height = 4.5, dpi = 300)

