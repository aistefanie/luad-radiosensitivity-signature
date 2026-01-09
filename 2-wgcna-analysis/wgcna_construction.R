
# wgcna筛选LUAD患者放疗与OS/vital关联的key module genes

library(WGCNA)
library(dplyr)
library(tidyr)
library(tibble)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

###### 1.数据准备 ######
exp_data <- read.csv("~/analysis/WGCNA/radiation_expression.csv", row.names = 1, check.names = FALSE) # 行基因名,列是样本ID
clinical_data <- read.csv("~/analysis/WGCNA/radiation_clinical.csv", check.names = FALSE)

# 关联临床性状
trait_data <- clinical_data %>% dplyr::select(bcr_patient_barcode, primary_therapy_outcome_success,vital_status,
    days_to_last_followup,days_to_death)   # 转换临床变量为WGCNA的数值格式
# 因子变量转换为数值型
trait_data$primary_therapy_outcome_success <- ifelse(trait_data$primary_therapy_outcome_success == "Complete Remission/Response", 1, 0)
trait_data$vital_status <- ifelse(trait_data$vital_status == "Dead", 1, 0)

# 计算总生存期OS_time (如果存在days_to_death则使用，否则使用days_to_last_followup)
trait_data$OS_time <- pmax(trait_data$days_to_death, trait_data$days_to_last_followup, na.rm = TRUE)
trait_data <- trait_data %>%
  dplyr::select(-days_to_death, -days_to_last_followup) %>% column_to_rownames("bcr_patient_barcode")
datTraits <- trait_data[substr(rownames(datExpr), 1, 12), ] # 顺序一致
rownames(datTraits) <- rownames(datExpr)
# 补全NA值
datTraits$OS_time[is.na(datTraits$OS_time)] <- median(datTraits$OS_time, na.rm = TRUE)
datTraits$vital_status[is.na(datTraits$vital_status)] <- 0

# datExpr转换
library(DESeq2)
colData <- data.frame(row.names = colnames(exp_data), condition = rep("A", ncol(exp_data)))
dds <- DESeqDataSetFromMatrix(countData = exp_data, colData = colData, design = ~ 1)
vsd <- vst(dds, blind = TRUE)
datExpr <- assay(vsd) # 方差稳定化

datExpr <- t(datExpr)
# 去除低方差基因 (选取前 5000 或 MAD top 75%)
mads <- apply(datExpr, 2, mad)
datExpr <- datExpr[, order(mads, decreasing = TRUE)[1:5000]]
# 检查缺失值
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK){
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

###### 2.WGCNA核心分析 ######
# sample clustering
sampleTree = hclust(dist(datExpr), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# 1. 确定最佳软阈值 (power)
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2)) # powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# 软阈值图
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red", cex = 0.9)
abline(h = 0.90, col = "red")
# Mean Connectivity图
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", 
     type = "n",
     main = "Mean Connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = powers, col = "red", cex = 0.9)

# 2. 构建共表达网络和模块
# 选择合适的软阈值 power，选择 R^2 > 0.9 的最小值
power <- sft$powerEstimate
power <- 13
# 构建network step
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
type = "unsigned"
corType = "pearson"
maxPOutliers = ifelse(corType=="pearson",1,0.05)

net=blockwiseModules(datExpr,power=13,
                       maxBlockSize=nGenes,TOMType=type,
                       corType=corType,maxPOutliers=maxPOutliers,
                       loadTOMs=TRUE,reassignThreshold=0,
                       minModuleSize=30,saveTOMs=TRUE,mergeCutHeight=0.25,
                       numericLabels=TRUE,pamRespectsDendro=FALSE,
                       saveTOMFileBase=paste0("expr-luad", ".tom"),verbose=3)
table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
table(moduleColors)

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(10, 8.5, 3, 3), las = 2,)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# 提取关键模块的基因列表(lightgreen、gray60)
color <- data.frame(moduleColors)
rownames(color) <- colnames(datExpr)
lightengreen_Genes <- subset(color,color$moduleColors=="lightgreen")
gray60_Genes <- subset(color,color$moduleColors=="grey60")
write.csv(lightengreen_Genes,"vital_module.csv",row.names = T)
write.csv(gray60_Genes,"OStime_module.csv",row.names = T)


####只用"vital_status", "OS_time"
target_traits <- c("vital_status", "OS_time")  # 明确需要保留的性状
moduleTraitCor_filtered <- moduleTraitCor[, target_traits, drop = FALSE]
moduleTraitPvalue_filtered <- moduleTraitPvalue[, target_traits, drop = FALSE]
textMatrix_filter <- paste(signif(moduleTraitCor_filtered, 2), 
  "\n(",signif(moduleTraitPvalue_filtered, 1),")", sep = "")
dim(textMatrix_filter) <- dim(moduleTraitCor_filtered)
datTraits_filter <- datTraits[, target_traits, drop = FALSE] 

pdf("Module-Trait Relationships.pdf", width = 6.5, height = 5.5)
par(mar = c(4, 8.5, 3, 3), las = 2,font.axis = 1)
labeledHeatmap(Matrix = moduleTraitCor_filtered,
  xLabels = names(datTraits_filter),
  yLabels = names(MEs), ySymbols = names(MEs),
  colorLabels = FALSE,colors = blueWhiteRed(50),  
  textMatrix = textMatrix_filter,setStdMargins = FALSE,
  cex.text = 0.5,zlim = c(-1, 1),
  main = "Module-Trait Relationships")
dev.off()


