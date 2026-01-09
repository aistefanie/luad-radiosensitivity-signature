
# ========== 基因和参数 ==========
genes_to_plot <- c("CXCL12", "CXCR4", "VEGFA", "KDR", "FN1", "ITGB1", "SPP1", "CD44")
genes_valid <- genes_to_plot[genes_to_plot %in% rownames(seu_filtered)]
celltypes <- seu_filtered$Cell_type.y
ct_levels <- sort(unique(celltypes))

# log-normalized 表达矩阵
expr_log <- GetAssayData(seu_filtered, layer = "data")[genes_valid, ]

# ========== 创建空星号矩阵 ==========
stars_matrix <- matrix("", nrow = length(genes_valid), ncol = length(ct_levels))
rownames(stars_matrix) <- genes_valid
colnames(stars_matrix) <- ct_levels

# ========== 显著性 & fold-change 筛选 ==========
for (g in genes_valid) {
  for (ct in ct_levels) {
    expr_in  <- expr_log[g, celltypes == ct]
    expr_out <- expr_log[g, celltypes != ct]
    p <- tryCatch(wilcox.test(expr_in, expr_out)$p.value, error = function(e) 1)
    expr_in_raw  <- expm1(expr_in)
    expr_out_raw <- expm1(expr_out)
    avg_in  <- mean(expr_in_raw)
    avg_out <- mean(expr_out_raw)
    logfc   <- log2((avg_in + 1e-6) / (avg_out + 1e-6)) 
    if (p < 0.001 & logfc > 1.5) {
      stars_matrix[g, ct] <- "***"
    } else {
      stars_matrix[g, ct] <- ""
    }
  }
}
# ========== 平均表达热图矩阵 ==========
avg_expr <- AverageExpression(seu_filtered, features = genes_valid, group.by = "Cell_type.y", layer = "data")$RNA
pdf("LigandReceptor_Heatmap.pdf", width = 8.5, height = 4)
pheatmap(avg_expr, scale = "row",
         color = colorRampPalette(c("steelblue3", "white", "firebrick3"))(100),
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = TRUE, show_colnames = TRUE,
         fontsize_row = 14, fontsize_col = 14,
         angle_col = 45, border_color = NA,
         labels_col = colnames(avg_expr), cellwidth = 55,
         display_numbers = stars_matrix, number_color = "black",
         main = "Expression of Ligand–Receptor Pairs Across Cell Types")
dev.off()
