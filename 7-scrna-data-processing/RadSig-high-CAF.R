# =====================================================
# Figure 10: RadSig-high CAF Subcluster Analysis
# =====================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

#------------------------------------------------------
# 1) 选择 CAF（Fibroblasts）并分组 RadSig-high/low
#------------------------------------------------------

fibro_seu <- subset(seu_filtered, subset = Cell_type.y == "Fibroblasts")

median_rad <- median(fibro_seu$RadSigScore1, na.rm = TRUE)

fibro_seu$RadSigGroup <- ifelse(fibro_seu$RadSigScore1 >= median_rad,
                                "RadSig-high", "RadSig-low")
fibro_seu$RadSigGroup <- factor(fibro_seu$RadSigGroup,
                                levels = c("RadSig-low", "RadSig-high"))

#------------------------------------------------------
# 2) CAF 子集降维（tSNE）
#------------------------------------------------------

fibro_seu <- NormalizeData(fibro_seu)
fibro_seu <- FindVariableFeatures(fibro_seu)
fibro_seu <- ScaleData(fibro_seu)
fibro_seu <- RunPCA(fibro_seu)
fibro_seu <- RunTSNE(fibro_seu, dims = 1:20)

#------------------------------------------------------
# 3) Fig10A: tSNE (RadSig-high / low)
#------------------------------------------------------

pA <- DimPlot(fibro_seu, reduction = "tsne", group.by = "RadSigGroup",
               pt.size = 0.5) +
  scale_color_manual(values = c("RadSig-low"  = "#1F77B4",
                                "RadSig-high" = "#FF7F0E")) +
  labs(title = "RadSig-high CAF vs RadSig-low CAF",
       x = "tSNE_1", y = "tSNE_2") +
  theme_classic(base_size = 13) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "right",
    axis.line = element_line(color = "black", linewidth = 0.6),
    panel.border = element_blank(),
    panel.grid = element_blank())

pdf("CAF_tSNE_RadSig.pdf", width = 5.5, height = 4)
print(p9A)
dev.off()

#------------------------------------------------------
# 4) Fig10B: 差异基因分析--火山图
#------------------------------------------------------

Idents(fibro_seu) <- "RadSigGroup"

deg <- FindMarkers(fibro_seu,
                   ident.1 = "RadSig-high",
                   ident.2 = "RadSig-low",
                   logfc.threshold = 0.25,
                   min.pct = 0.1)

deg$gene <- rownames(deg)
lfc <- if ("avg_log2FC" %in% colnames(deg)) "avg_log2FC" else "avg_logFC"
deg$logFC <- deg[[lfc]]
deg$negLog10P <- -log10(deg$p_val_adj + 1e-300)

deg$change <- "NS"
deg$change[deg$logFC > 0.5 & deg$p_val_adj < 0.05] <- "Up"
deg$change[deg$logFC < -0.5 & deg$p_val_adj < 0.05] <- "Down"

top_genes <- deg %>% arrange(p_val_adj) %>% head(15) %>% pull(gene)

pB <- ggplot(deg, aes(logFC, negLog10P, color = change)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("Up" = "#FF7F0E",
                                "Down" = "#1F77B4",
                                "NS" = "grey70")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 2, color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey50") +
  ggrepel::geom_text_repel(data = subset(deg, gene %in% top_genes),
                           aes(label = gene), size = 3) +
  labs(title = "DEGs: RadSig-high vs RadSig-low CAF",
       x = "log2 Fold Change", y = "-log10 Adjusted P-value",color = "") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )
pdf("DEG_Volcano_RadSig.pdf", width = 6.5, height = 5)
print(p9B)
dev.off()

#------------------------------------------------------
# 5) Fig10C: DotPlot (CAF markers)
#------------------------------------------------------
marker_genes <- c("CXCL12", "SPP1", "MIF", "IL6",
                  "FN1", "COL1A1", "POSTN", "TGFBI")
pC <- DotPlot(fibro_seu, features = marker_genes, group.by = "RadSigGroup") +
  scale_color_gradientn(colors = c("#2166AC", "white", "#B2182B")) +
  labs(title = "CAF Marker Expression",x = "", y = "Group") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"), 
        axis.title.y = element_text(color = "black") 
  )

pdf("CAF_Marker_DotPlot.pdf", width = 6.5, height = 5)
print(p9C)
dev.off()

#------------------------------------------------------
# 6) Fig10D: GO/KEGG 富集分析（RadSig-high CAF 上调基因）
#------------------------------------------------------
up_genes <- deg %>% filter(change == "Up") %>% arrange(p_val_adj) %>% pull(gene)
entrez <- bitr(up_genes, fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)$ENTREZID
ego <- enrichGO(gene = entrez,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",
                readable = TRUE)
ego_top <- ego %>% as.data.frame() %>% head(10)

pD <- ego_top %>%
  mutate(Description = factor(Description, levels = rev(Description))) %>%
  ggplot(aes(x = Description, y = Count)) +
  geom_bar(stat = "identity", fill = "#D37A6E", width = 0.7) +
  coord_flip() + labs(x = NULL, y = "Gene Count") + 
  theme_bw(base_size = 14) + theme(
    axis.text = element_text(color = "black", size = 17),
    axis.line = element_line(color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    plot.margin = margin(10, 10, 10, 10))
pdf("GO_Enrichment_RadSigHigh.pdf", width = 8, height = 4.5)
print(p9D)
dev.off()


# ============================================================
# CAF Subtype Annotation (iCAF / myCAF / ECM-CAF)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
})

# -------------------------------------------------------------
# 1) 定义 CAF subtype marker sets
# -------------------------------------------------------------
iCAF_markers  <- c("IL6", "CXCL12", "CCL2")
myCAF_markers <- c("ACTA2", "TAGLN", "MYL9")
ECM_markers   <- c("FN1", "COL1A1", "COL3A1", "TNC")

marker_list <- list(
  iCAF = iCAF_markers,
  myCAF = myCAF_markers,
  ECM_CAF = ECM_markers
)
marker_all <- unique(unlist(marker_list))
marker_all <- marker_all[marker_all %in% rownames(fibro_seu)]

# -------------------------------------------------------------
# 2) CAF subtype scoring（AddModuleScore）
# -------------------------------------------------------------
fibro_seu <- AddModuleScore(
  fibro_seu,
  features = marker_list,
  name = "CAFsubtype"
)

colnames(fibro_seu@meta.data)[ grepl("CAFsubtype", colnames(fibro_seu@meta.data)) ] <-
  c("iCAF_score", "myCAF_score", "ECM_score")

# 按最大得分指定 CAF subtype
fibro_seu$CAF_subtype <- apply(
  fibro_seu@meta.data[, c("iCAF_score", "myCAF_score", "ECM_score")],
  1, function(x) names(which.max(x))
)
fibro_seu$CAF_subtype <- factor(
  fibro_seu$CAF_subtype,
  levels = c("iCAF_score", "myCAF_score", "ECM_score"),
  labels = c("iCAF", "myCAF", "ECM-CAF")
)

# 统一颜色
CAF_colors <- c("iCAF"="#1f77b4", "myCAF"="#2ca02c", "ECM-CAF"="#ff7f0e")

# -------------------------------------------------------------
# 3) Fig10E：CAF subtype marker 热图
# -------------------------------------------------------------

pdf("CAF_subtype_heatmap.pdf",width = 6, height = 5, useDingbats = FALSE, family = "Helvetica")
DoHeatmap(object = fibro_seu, features = marker_all, 
  group.by = "CAF_subtype",label = FALSE ) + scale_fill_gradientn(colors = c("#2166AC","white","firebrick3"),
    name = "Expression") + scale_colour_manual(values = CAF_colors,guide = "none"
  ) + labs(title="CAF Subtype Marker Expression") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(size = 16,color = "black"),
    axis.text.x.top = element_blank(), 
    axis.ticks.x.top = element_blank(),
    axis.title = element_blank(),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13))
dev.off()

# -------------------------------------------------------------
# 4) Panel F：CAF subtype 分布 (tSNE/UMAP)
# -------------------------------------------------------------

pdf("CAF_subtype_tSNE.pdf", width=6, height=4.5)
DimPlot(fibro_seu, reduction="tsne", group.by="CAF_subtype") +
  scale_color_manual(values=CAF_colors) +
  labs(title="CAF Subtypes on tSNE", x="tSNE_1", y="tSNE_2") +
  theme_classic(base_size=13) +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        axis.line = element_line(color="black"),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13, face = "bold")
  )
dev.off()

# -------------------------------------------------------------
# 5) Fig10G：各 subtype 的 RadSigScore 分布（箱线图）
# -------------------------------------------------------------
pdf("RadSigScore_by_subtype.pdf", width=6, height=4.5)
ggplot(fibro_seu@meta.data, aes(x = CAF_subtype, y = RadSigScore1, fill = CAF_subtype)) +
  geom_boxplot(outlier.alpha=0.1) +
  scale_fill_manual(values = CAF_colors) +
  labs(title="RadSigScore across CAF Subtypes",
       x="CAF subtype", y="RadSigScore") +
  theme_classic(base_size=13) +
  theme(axis.text.x = element_text(size = 16, color = "black"),
        legend.position="none",
        axis.line = element_line(color="black"))
dev.off()

# -------------------------------------------------------------
# 6) Fig10H：各 subtype 中 RadSig-high CAF 的比例
# -------------------------------------------------------------

fibro_seu$RadSigGroup <- ifelse(
  fibro_seu$RadSigScore1 >= median(fibro_seu$RadSigScore1),
  "RadSig-high", "RadSig-low"
)

df_bar <- fibro_seu@meta.data %>%
  group_by(CAF_subtype, RadSigGroup) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

pdf("RadSigHigh_fraction.pdf", width=6, height=4.5)
ggplot(df_bar, aes(x = CAF_subtype, y = freq, fill = RadSigGroup)) +
  geom_bar(stat="identity", position="fill", color="black") +
  scale_fill_manual(values=c("RadSig-high"="#ff7f0e", "RadSig-low"="#1f77b4")) +
  labs(x = NULL, y="Proportion") +
  theme_classic(base_size=13) +
  theme(axis.text.x = element_text(size = 16),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        plot.title = element_text(hjust=0.5, face="bold"),
        axis.line = element_line(color="black"))
dev.off()

# ===========================================
# Fig10I: CAF Functional Programs vs RadSigScore
# ===========================================
library(reshape2)
library(ggplot2)

# 定义 CAF 功能基因集
signatures <- list(
  iCAF      = c("IL6", "CXCL12", "CCL2"),
  myCAF     = c("ACTA2", "TAGLN", "MYL9"),
  ECM_CAF   = c("FN1", "COL1A1", "COL3A1"),
  Chemokine = c("CXCL12", "CCL2", "CCL5"),
  TGFb      = c("TGFB1", "TGFBI", "COL1A1"),
  IFN       = c("STAT1", "ISG15", "IFIT2")
)
signatures <- lapply(signatures, function(x) x[x %in% rownames(fibro_seu)])
for (i in names(signatures)) {
  fibro_seu <- AddModuleScore(
    fibro_seu,
    features = list(signatures[[i]]),
    name = paste0("Sig_", i)
  )
}  # 逐一 AddModuleScore
sig_cols <- grep("^Sig_", colnames(fibro_seu@meta.data), value = TRUE)
clean_names <- sig_cols %>%
  gsub("^Sig_", "", .) %>%
  gsub("1$", "", .)
sig_mat <- fibro_seu@meta.data[, sig_cols]
colnames(sig_mat) <- clean_names
sig_mat$RadSigGroup <- fibro_seu$RadSigGroup
df_radar <- sig_mat %>%
  group_by(RadSigGroup) %>%
  summarise(across(all_of(clean_names), mean)) %>%
  as.data.frame()

df_plot <- data.frame(
  Group = c("RadSig-high", "RadSig-low"),
  iCAF = df_radar[, "iCAF"],
  myCAF = df_radar[, "myCAF"],
  ECM_CAF = df_radar[, "ECM_CAF"],
  Chemokine = df_radar[, "Chemokine"],
  TGFb = df_radar[, "TGFb"],
  IFN = df_radar[, "IFN"]
)
df_line <- melt(df_plot, id.vars="Group",
                variable.name="Program",
                value.name="Score")

pdf("Fig10I_CAF_FunctionalPrograms_RadSig.pdf", width=7, height=5)
ggplot(df_line, aes(x = Program, y = Score, group = Group, color = Group)) +
  geom_line(size = 1.4) + geom_point(size = 3) +
  scale_color_manual(values = c("RadSig-high" = "#FF7F0E", "RadSig-low"  = "#1F77B4"
  )) + theme_classic(base_size = 14) + theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 16),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) + labs(x = "CAF Functional Program", y = "Signature Score", color = "")
dev.off()

