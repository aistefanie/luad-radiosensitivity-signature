###############################################################
# 0. 环境准备
###############################################################
library(maftools)
library(dplyr)
library(ggplot2)
library(reshape2)
library(TCGAbiolinks)
library(stringr)
library(tidyr)

# 1. 下载MAF文件
query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Simple Nucleotide Variation",
                  data.type = "Masked Somatic Mutation",
                  workflow.type = "MuTect2 Variant Aggregation and Masking")
GDCdownload(query)
maf <- GDCprepare(query)
maf <- read.maf(maf)

# 2. 加载 RadSig 高低分组
group <- read.csv("RiskGroups.csv", header=TRUE, stringsAsFactors=FALSE)性
group$SampleID <- toupper(group$SampleID)
maf_ids <- unique(maf@data$Tumor_Sample_Barcode)
group$PatientID <- str_extract(group$SampleID, "^TCGA-\\w\\w-\\w\\w\\w\\w")
convert_to_maf_id <- function(pid){
  idx <- grepl(pid, maf_ids)
  maf_ids[idx][1]
}
group$MAF_ID <- sapply(group$PatientID, convert_to_maf_id)
group$SampleID <- group$MAF_ID
group$MAF_ID <- NULL
group$PatientID <- NULL
write.csv(group, "sample_group_FIXED.csv", row.names = FALSE)

# 按组提取样本
high_samples <- group$SampleID[group$RadSigGroup == "High"]
low_samples  <- group$SampleID[group$RadSigGroup == "Low"]

# 3. 分组 MAF
maf_high <- subsetMaf(maf, tsb = high_samples)
maf_low  <- subsetMaf(maf, tsb = low_samples)

# 4. 绘制 Oncoplot（高低分组各 20 个高频突变基因）
pdf("Oncoplot_RadSig_High.pdf", width=8, height=6)
oncoplot(maf = maf_high,top = 20,draw_titv = TRUE,removeNonMutated = TRUE,
         titleText = "RadSig-High Mutation Landscape")
dev.off()

pdf("Oncoplot_RadSig_Low.pdf", width=8, height=6)
oncoplot(maf = maf_low,top = 20,draw_titv = TRUE,removeNonMutated = TRUE,
         titleText = "RadSig-Low Mutation Landscape")
dev.off()

# 5. 全队列 Oncoplot
pdf("Oncoplot_All_TCGA_LUAD.pdf", width=10, height=6)
oncoplot(maf = maf, top = 20, draw_titv = TRUE, removeNonMutated = TRUE,
  titleText = "TCGA-LUAD Mutation Landscape")
dev.off()

# 6. TMB 计算并比较 RadSig 高低组
tmb_res <- tmb(maf)
tmb_res$Group <- group$RadSigGroup[match(tmb_res$Tumor_Sample_Barcode, group$SampleID)]
tmb_res2 <- tmb_res[!is.na(tmb_res$Group), ]
table(tmb_res2$Group)
pdf("TMB_Comparison.pdf", width=6, height=5)
ggplot(tmb_res2, aes(x = Group, y = total_perMB, fill = Group)) +
  geom_boxplot() + geom_jitter(width=0.2, alpha=0.5) +
  scale_fill_manual(values=c("High"="#E64B35", "Low"="#4DBBD5")) + theme_bw()
dev.off()

t_test_result <- wilcox.test(total_perMB ~ Group, data=tmb_res) # 统计检验

# 7. RadSig 六个基因的突变情况
RadSig_genes <- c("VEGFC", "SERPINE1", "CHEK1", "MEG3", "SLC16A1", "PRRX1")
rad_mut <- subsetMaf(maf, genes=RadSig_genes, includeSyn = TRUE)
write.csv(rad_mut@data, "RadSig_gene_mutations.csv", row.names=FALSE)

####  FigC.Violin Plot & FigD.TMB Density Distribution ####
library(ggplot2)
library(ggpubr)

p1 <- ggplot(tmb_res2, aes(x = Group, y = total_perMB, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.1, color = "black", alpha = 0.7) +
  geom_jitter(width = 0.12, alpha = 0.4, size=1) +
  scale_fill_manual(values = c("High" = "#E64B35", "Low" = "#4DBBD5")) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", size = 1),
    axis.text.x  = element_text(size = 16, color = "black"),
    axis.text.y  = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16, color = "black")
  ) +
  labs(y = "TMB (mut/Mb)",
       title = "Violin Plot of TMB between RadSig Groups") +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",
                     label.y = max(tmb_res2$total_perMB) * 1.1)

p2 <- ggplot(tmb_res2, aes(x = total_perMB, color = Group, fill = Group)) +
  geom_density(alpha = 0.3, size=1.2) +
  scale_color_manual(values = c("High"="#E64B35","Low"="#4DBBD5")) +
  scale_fill_manual(values = c("High"="#E64B35","Low"="#4DBBD5")) +
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", size = 1),
    axis.text.x  = element_text(size = 16, color = "black"),
    axis.text.y  = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16, color = "black")
  ) +
  labs(x = "TMB (mut/Mb)",
       title = "TMB Density Distribution between RadSig Groups")

ggsave("TMB_Violin_Nature.pdf", p1, width=5, height=5)
ggsave("TMB_Density.pdf", p2, width=5, height=5)


####  FigE.Forest Plot + Mutation Count Barplot ####

comparison <- mafCompare(m1 = maf_high, m2 = maf_low, m1Name="High", m2Name="Low")
df2 <- df %>% select(Gene) %>%
  left_join(comparison$results, by = c("Gene" = "Hugo_Symbol")) %>%
  select(Gene, High, Low)
df2_long <- df2 %>% pivot_longer(cols = c("High", "Low"),
                                 names_to = "Group",
                                 values_to = "Count")
p_bar <- ggplot(df2_long, aes(x = Count, y = reorder(Gene, Count), fill = Group)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("High"="#E64B35", "Low"="#4DBBD5")) +
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", size = 1),
    axis.title.y = element_blank(),
    legend.title = element_blank(),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.text.x  = element_text(size = 14, color = "black"),
    axis.text.y  = element_text(size = 14, color = "black")
  ) +
  labs(x = "Mutation Count")
library(patchwork)
p_forest <- ggplot(df, aes(x = logOR, y = reorder(Gene, logOR))) +
  geom_point(size = 3, color = "#E64B35FF") +
  geom_errorbarh(aes(xmin = logOR - 0.1, xmax = logOR + 0.1),
                 height = 0.25, color = "gray40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  theme_bw(base_size = 14) +
  labs(x = "log2(Odds Ratio)") +
  theme(
    axis.title.y = element_blank(),
    panel.border = element_rect(color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.text.x  = element_text(size = 14, color = "black"),
    axis.text.y  = element_text(size = 14, color = "black")
  )
combined_plot <- p_forest | p_bar 
pdf("Driver_Mutation_Enrichment.pdf", width=9, height=5)
combined_plot
dev.off()
