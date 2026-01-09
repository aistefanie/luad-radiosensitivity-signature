library(GSVA)
library(ggplot2)
library(reshape2)
# --- 16种免疫细胞浸润分析 --- #
df <- read.csv("LUAD_expMatrix.csv", header = TRUE)
rownames(df) <- make.names(df[, 1], unique = TRUE)
df <- df[, -1]
dim(df)

expr_matrix <- as.matrix(df)
data <- read.csv("MsigDB/16_immune_gene_sets.csv",header = T)
gene_sets <- split(data$marker, data$cell_type)

param <- GSVA::gsvaParam(exprData = expr_matrix,geneSets = gene_sets,
                         kcdf = "Gaussian", minSize = 1, maxSize = 5000)
res <- GSVA::gsva(param, verbose=TRUE, BPPARAM=SerialParam(progressbar=TRUE))

write.csv(res,"immune_cell.csv")

# Melt the results into a long format suitable for ggplot2
ssgsea_results_melted <- melt(ssgsea_results, 
                              varnames = c("Cell_Type", "Sample"), 
                              value.name = "Enrichment_Score")
# Plotting boxplot

# --- 13种免疫功能评分 --- #
data <- read.csv("MsigDB/13_Immune_function_sets.csv",header = F)
data <- t(data)
colnames(data) <- data[1,]
data <- data[-1,]
gene_sets <- split(data, colnames(data))
param <- GSVA::gsvaParam(exprData = expr_matrix,geneSets = gene_sets,
                         kcdf = "Gaussian", minSize = 1, maxSize = 5000)
res1 <- GSVA::gsva(param, verbose=TRUE, BPPARAM=SerialParam(progressbar=TRUE))

write.csv(res1,"immune_funtion.csv")


