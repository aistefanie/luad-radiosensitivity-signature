library(Seurat)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(pheatmap)
library(viridis)
library(CellChat)
library(patchwork)

# step 1: 数据加载和预处理 ####

setwd("/root/analysis/scRNA_data/")
# rds 文件是一个 Seurat 对象
seu <- readRDS("GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds")
seu <- CreateSeuratObject(counts = seu, project = "scRNAseq")

# anno 文件是包含 Barcode, Cell_type, Cell_subtype 的 CSV 或 TSV
anno_df <- read.delim("GSE131907_Lung_Cancer_cell_annotation.txt",header = T,sep = "\t") 

# 确保 anno_df 包含 'Barcode' 和 'Cell_type' 或 'Cell_subtype'
# 用 'Cell_type' 进行最终注释
anno_df <- anno_df %>% 
  select(Barcode = Index, Cell_type = Cell_type) # Index 列是 Barcode

# 检查 Seurat 对象和注释文件的 Barcode 是否匹配
if (all(Cells(seu) %in% anno_df$Barcode)) {
  # 将注释信息添加到 Seurat 对象的 metadata 中
  metadata_df <- seu@meta.data %>% 
    rownames_to_column("Barcode") %>%
    left_join(anno_df, by = "Barcode") %>%
    column_to_rownames("Barcode")
  
  seu <- AddMetaData(seu, metadata_df)
} else {
  stop("Barcodes in Seurat object and annotation file do not match exactly. Check your files.")
}

# step 2: 运行降维和聚类 (rds 文件中还没有 tSNE/UMAP 坐标和 Clusters) ####

# 仅在 rds 文件中缺失 tSNE 坐标和 Cluster 信息时执行
if (is.null(seu@reductions$tsne)) {
  seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
  seu <- RunPCA(seu, features = VariableFeatures(object = seu))
  seu <- FindNeighbors(seu, dims = 1:20)
  seu <- FindClusters(seu, resolution = 0.8) # Seurat cluster in metadata is 'seurat_clusters'
  seu <- RunTSNE(seu, dims = 1:20)
}
# 保存
save(seu,file = "seu.Rdata") # load(file = "seu.Rdata")

# ----------------------------------------------------------------------
# step 3: 准备绘图函数和颜色 ####

# 提取 Seurat 对象中的 tSNE 坐标和聚类/细胞类型信息
plot_data <- data.frame(
  tSNE_1 = seu@reductions$tsne@cell.embeddings[, 1],
  tSNE_2 = seu@reductions$tsne@cell.embeddings[, 2],
  Cluster = seu$seurat_clusters,           # 使用 Seurat 默认的聚类结果
  Cell_type = seu$Cell_type.y                # 使用细胞类型注释
)

plot_data_filtered <- plot_data %>% filter(!is.na(Cell_type)) # 去掉NA的细胞
# 配色
sci_colors <- c("B lymphocytes" = "#1F78B4",  "Endothelial cells" = "#A6CEE3", 
                "Epithelial cells" = "#529B52", "Fibroblasts" = "#EDC948",     
                "MAST cells" = "#FB9A99", "Myeloid cells" = "#9B529B",    
                "T lymphocytes" = "#FFA500","NK cells" = "#DB7093")
cell_type_colors <- sci_colors

# ----------------------------------------------------------------------
# Figure A: tSNE 降维图 (按 "cluster" 着色) ####

p1 <- ggplot(plot_data_filtered, aes(x = tSNE_1, y = tSNE_2, color = Cluster)) +
  geom_point(size = 0.8, alpha = 0.8) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_blank(), 
        axis.line.x = element_line(colour = "black", linewidth = 0.5), 
        axis.line.y = element_line(colour = "black", linewidth = 0.5)) +
  labs(title = "Clusters") +
  guides(color = guide_legend(override.aes = list(size = 4)))

# 添加聚类标签
cluster_centers <- plot_data_filtered %>%
  group_by(Cluster) %>%
  summarise(tSNE_1 = median(tSNE_1), tSNE_2 = median(tSNE_2))

p1 <- p1 + geom_text(data = cluster_centers, aes(label = Cluster), color = "black", size = 4)

ggsave(filename = "tSNE.pdf",
       plot = p1, device = "pdf",width = 8.5, height = 6,units = "in",
       dpi = 600, useDingbats = FALSE)

# ----------------------------------------------------------------------
# Figure B: tSNE 降维图 (按 "细胞类型" 着色) ####

p2 <- ggplot(plot_data_filtered, aes(x = tSNE_1, y = tSNE_2, color = Cell_type)) +
  geom_point(size = 0.8, alpha = 0.8) +
  scale_color_manual(values = cell_type_colors) +
  theme_bw(base_size = 18) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  plot.title = element_text(hjust = 0.5),
  panel.border = element_blank(), 
    axis.line.x = element_line(colour = "black", linewidth = 0.5), 
    axis.line.y = element_line(colour = "black", linewidth = 0.5)) +
  labs(title = "Celltype") + guides(color = guide_legend(override.aes = list(size = 4)))

ggsave(filename = "tSNE_Celltype_Distribution.pdf",
  plot = p2,device = "pdf",
  width = 8.5, height = 6,units = "in",
  dpi = 600, useDingbats = FALSE)

# ----------------------------------------------------------------------
# Figure: 细胞类型间的相互作用网络图 ####

seu_filtered <- subset(seu, subset = !is.na(Cell_type)) # 去掉 NA 分组细胞

cellchat_obj <- createCellChat(seu_filtered, group.by = "Cell_type.y")
cellchat_obj <- setIdent(cellchat_obj, ident.use = "Cell_type.y")
CellChatDB <- CellChatDB.human
cellchat_obj@DB <- CellChatDB
cellchat_obj <- subsetData(cellchat_obj)
cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
cellchat_obj <- computeCommunProb(cellchat_obj)
cellchat_obj <- aggregateNet(cellchat_obj)
# 绘制网络图
groupSize <- as.numeric(table(cellchat_obj@idents))

pdf("CellChat_Net.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), xpd = TRUE)
p3 <- netVisual_circle(cellchat_obj@net$count, vertex.weight = groupSize,
                       weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
p4 <- netVisual_circle(cellchat_obj@net$weight, vertex.weight = groupSize,
                       weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights")
dev.off()
par(mfrow = c(1, 1), xpd = FALSE)

# ----------------------------------------------------------------------
# figure D:  distribution of radsig in all cell types by counts ####

counts <- LayerData(seu_filtered, layer = "counts")[target_genes, ]
counts <- as.data.frame(t(counts))
counts$celltype <- seu_filtered@meta.data$Cell_type.y
counts$cluster <- seu_filtered@meta.data$seurat_clusters

# 分层抽样 #
seu_filtered$cellcode <- rownames(seu_filtered@meta.data)
strat_samp <- seu_filtered@meta.data %>% group_by(Cell_type.y) %>% sample_frac(size = .01)
obj_sample <- seu_filtered[,strat_samp$cellcode]
mat <- GetAssayData(obj_sample, slot = "data", assay = "RNA")[genes_to_plot, ]
nonzero_cells <- colnames(mat)[Matrix::colSums(mat) > 0]
obj_sample <- obj_sample[, nonzero_cells]
obj_sample <- obj_sample[, order(obj_sample$Cell_type.y)]

# 注释信息排序 #
celltype_counts <- table(obj_sample@meta.data$Cell_type.y) %>% 
  as.data.frame() %>%
  setNames(c("celltype", "count")) %>%
  arrange(desc(count))  # 按count从多到少排序
celltype_order <- celltype_counts$celltype # 排序后的细胞类型顺序
print(celltype_order)
annotation <- as.data.frame(obj_sample$Cell_type.y)
colnames(annotation) <- "celltype"
annotation <- annotation %>%
  mutate(celltype = factor(celltype, levels = celltype_order)) %>%
  arrange(celltype)
sorted_cell_names <- rownames(annotation)
mat4 = as.matrix(LayerData(obj_sample, layer = "counts")[target_genes, ])
mat4 <- mat4[, sorted_cell_names]  # 确保列顺序与 annotation 一致
# Z-score标准化 #
mat4 = t(scale(t(mat4)))
mat4 <- mat4[complete.cases(mat4), ] # 筛选出不包含任何NaN的行

library(ComplexHeatmap)
library(circlize)
library(grid)

ann_colors=list()
annotation_colors = c("#529B52","#9B529B","#EDC948","#FFA500","#1F78B4",
                      "#A6CEE3","#FB9A99","#DB7093")
names(annotation_colors) <- levels(factor(annotation$celltype))
ann_colors[["celltype"]] <- annotation_colors
bk=c(seq(0, 1, by = 0.01))
col_fun <- colorRamp2(bk, colorRampPalette(colors = c("#FFFFFF", "#B84D64"))(length(bk)))

col_anno <- HeatmapAnnotation(
  df = annotation,
  col = ann_colors, 
  show_legend = TRUE,
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 14, fontface = "bold"),
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  )
)

ht <- Heatmap(
  mat4,
  name = "Expression",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 12, fontface = "italic", col = "black"),
  top_annotation = col_anno,
  
  width  = unit(pmin(ncol(mat4) * 0.6, 180), "mm"),
  height = unit(pmin(nrow(mat4) * 8, 120), "mm"),
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 10, fontface = "plain"),  
    labels_gp = gpar(fontsize = 9),
    legend_width = unit(3, "mm"),
    legend_height = unit(30, "mm")
  )
)

pdf("Radsig_distribution_celltypes.pdf", width = 8, height = 5, useDingbats = FALSE)
draw(ht,heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legends = TRUE,
     legend_title_gp = gpar(fontsize = 10, fontface = "plain"),
     legend_labels_gp = gpar(fontsize = 10),
     padding = unit(c(3, 22, 3, 20), "mm"),   # 上右下左留白
     newpage = TRUE
)
dev.off()

# ----------------------------------------------------------------------
# figure E: average expression of radsig in all cell types

# 基因列表
target_genes <- c("VEGFC", "SERPINE1", "CHEK1","MEG3","SLC16A1","PRRX1")

exp_avg <- AverageExpression(seu_filtered, group.by = "Cell_type.y", assays = "RNA")[["RNA"]]
seu_filtered$cellcode <- rownames(seu_filtered@meta.data)
strat_samp <- seu_filtered@meta.data %>% group_by(Cell_type.y) %>% sample_frac(size = .03)
obj_sample <- seu_filtered[,strat_samp$cellcode]
obj_sample <- obj_sample[,order(obj_sample@meta.data$Cell_type.y)]
annotation <- as.data.frame(obj_sample$Cell_type.y)
colnames(annotation) <- "celltype"
ann_colors=list()
annotation_colors = c("#1F78B4","#A6CEE3","#529B52","#EDC948","#FB9A99","#9B529B","#DB7093","#FFA500")
names(annotation_colors) <- levels(factor(annotation$celltype))
ann_colors[["celltype"]] <- annotation_colors
cell_types <- colnames(exp_avg)

# 为每种细胞类型创建对应的颜色映射
name_mapping <- c("B lymphocytes" = "B lymphocytes",
  "Endothelial cells" = "Endothelial cells","Epithelial cells" = "Epithelial cells",
  "Fibroblasts" = "Fibroblasts","MAST cells" = "MAST cells","Myeloid cells" = "Myeloid cells",
  "T lymphocytes" = "T lymphocytes","NK cells" = "NK cells")
# 重构 annotation 数据框，行名与 exp_avg 的列名一致
annotation <- data.frame(
  celltype = factor(name_mapping[cell_types], levels = name_mapping[cell_types]),
  row.names = cell_types
)

# 从ann_colors 中提取需要的颜色
ann_colors <- list(celltype = annotation_colors[name_mapping[cell_types]])
p5 <- pheatmap::pheatmap(
  exp_avg[target_genes, ], scale = "row",
  cellwidth = 16, cellheight = 16,
  fontsize = 10, fontsize_row = 12, fontsize_col = 12, 
  annotation = annotation,
  annotation_colors = ann_colors,
  cluster_rows = F, cluster_cols = F,
  color = colorRampPalette(c("#4292c6", "white", "#e2413f"))(50),
  annotation_legend = T
)
ggsave("Avg_Expression_Heatmap.pdf",p5,height = 6,width = 5)

### tsne featureplot ###
p6 <- FeaturePlot(seu_filtered, features = c("VEGFC", "SERPINE1","CHEK1","MEG3","SLC16A1","PRRX1"),
                       reduction = "tsne")
pdf("radsig_tsne.pdf", width = 8, height = 9, useDingbats = FALSE)
print(p6)
dev.off()

### DotPlot ###
p6_dotplot <- DotPlot(seu_filtered, features = target_genes,
        group.by = "Cell_type.y", cols = c("#1F78B4","#FF7F00"), dot.scale = 6)+
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf("DotPlot.pdf", width = 5, height = 4, useDingbats = FALSE)
print(p6_dotplot)
dev.off()


