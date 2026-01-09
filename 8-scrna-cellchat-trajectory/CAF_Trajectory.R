library(Seurat)
library(SingleCellExperiment)
library(slingshot)

##### Fig12A #####
# 提取counts、meta、tSNE
counts   <- GetAssayData(fibro_seu, assay = "RNA", slot = "counts")
meta     <- fibro_seu@meta.data
tsne <- fibro_seu@reductions$tsne@cell.embeddings
sce <- SingleCellExperiment(assays  = list(counts = counts), colData = meta)
reducedDims(sce)$TSNE <- tsne
sce$cluster <- sce$seurat_clusters
tab <- table(sce$cluster)
tiny <- names(tab[tab < 5])
sce_filt <- sce[, !(sce$cluster %in% tiny)]
tiny <- names(tab[tab < 3])
sce_filt <- sce[, !(sce$cluster %in% tiny)]
sce_filt <- slingshot(sce_filt, clusterLabels = "cluster", reducedDim = "TSNE")
# 可视化
pdf("Fig12A_CAF_Trajectory.pdf", width = 6, height = 5)
cols <- colorRampPalette(c("navy","skyblue","gold","red3"))(100)
rad_cols <- cols[cut(fibro_seu$RadSigScore1[colnames(sce_filt)], breaks = 100)]
plot(reducedDims(sce_filt)$TSNE, col = rad_cols,
     pch = 16,cex = 0.6, xlab = "tSNE_1", ylab = "tSNE_2",
     main = "RadSigScore Along CAF Trajectory")
lines(SlingshotDataSet(sce_filt), lwd = 2, col = "darkgreen")
dev.off()

##### Fig12B #####
set.seed(123)
k <- 3
pt <- fibro_seu$pseudotime
idx_use <- which(!is.na(pt)) 
pt_use  <- pt[idx_use]
km <- kmeans(pt_use, centers = k)
fibro_seu$CAF_state <- NA
fibro_seu$CAF_state[idx_use] <- km$cluster
fibro_seu$CAF_state <- factor(
  fibro_seu$CAF_state,
  levels = 1:k,
  labels = c("Early-CAF", "Mid-CAF", "Late-CAF")
)
table(fibro_seu$CAF_state, useNA = "ifany")
df <- data.frame(RadSigGroup = fibro_seu$RadSigGroup, CAF_state   = fibro_seu$CAF_state)

df_plot <- df %>%
  filter(!is.na(CAF_state)) %>%
  group_by(RadSigGroup, CAF_state) %>%
  summarise(n = n()) %>%
  group_by(RadSigGroup) %>%
  mutate(Percent = n / sum(n) * 100)

state_cols <- c(
  "Early-CAF"="#1f78b4",
  "Mid-CAF"  ="#33a02c",
  "Late-CAF" ="#e31a1c"
)
pdf("Fig12B_CAFstate_barplot.pdf", width=5, height=4)

ggplot(df_plot, aes(x = RadSigGroup, y = Percent, fill = CAF_state)) +
  geom_bar(stat="identity", width=0.7, color="black") +
  scale_fill_manual(values = state_cols) +
  labs(x="", y="Percentage of CAF states (%)",
       title="CAF State Composition in RadSig-high vs RadSig-low") +
  theme_classic(base_size=15) +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.title = element_blank(),
        legend.text = element_text(size=12)
  )
dev.off()

##### Fig12C #####
library(pheatmap)
radsig_genes <- c("CXCL12", "SPP1", "MIF", "IL6", "FN1", "COL1A1", "SERPINE1", "VEGFC")
pt <- slingPseudotime(sce_filt)[, 1]
pseudotime <- rep(NA, length(colnames(fibro_seu)))
names(pseudotime) <- colnames(fibro_seu)
pseudotime[names(pt)] <- pt
expr <- GetAssayData(fibro_seu, slot = "data")[radsig_genes, ]
ord <- order(pseudotime, na.last = NA)
expr_ord <- expr[, ord]
pt_ord   <- pseudotime[ord]

pdf("Fig12C_RadSig_Pseudotime_Heatmap.pdf", width = 6, height = 4.8)
pheatmap(expr_ord,cluster_rows = F, cluster_cols = F,
  color = colorRampPalette(c("steelblue3","white","firebrick3"))(100),
  border_color = NA, show_colnames = FALSE, fontsize_row = 14,
  main = "RadSig Genes Along CAF Pseudotime", labels_row = radsig_genes)
dev.off()

#####  Fig12D #####
pathways.show <- c("MIF", "SPP1", "CXCL", "FN1")
pdf("Fig12D_CellChat_Bubble.pdf", width = 4.5, height = 3.8)
netVisual_bubble(cellchat_obj, sources.use = "Fibroblasts",
                 targets.use = c("T lymphocytes", "Myeloid cells", "NK cells"),
                 signaling = pathways.show, font.size = 12,
                 remove.isolate = TRUE, angle.x = 45)
dev.off()
