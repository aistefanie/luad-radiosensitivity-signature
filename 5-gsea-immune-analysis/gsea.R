################################################################################
##  GSEA : 对每个目标基因，找到其所在的、最显著的6 个 GO 通路和 6 个 KEGG 通路并绘制
##  依赖文件: LUAD_expMatrix.csv, RiskGroups_by_median.csv
################################################################################

library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(tidyverse)
library(RColorBrewer)

OrgDb_pkg <- org.Hs.eg.db

# 2. 定义目标基因和分析设置 ----------------------------------------------------
TARGET_GENES <- c("VEGFC", "SERPINE1", "CHEK1", "MEG3", "SLC16A1", "PRRX1")
GO_ONTOLOGY <- "BP" # GO 富集本体: Biological Process
ONTOLOGIES <- c("GO", "KEGG") 

# 3. 数据读取与整合 ----------------------------------------------------------
# 3.1 读取表达 counts 矩阵
counts_matrix_raw <- read_csv("LUAD_expMatrix.csv") %>% 
  rename(Symbol = 1) 

# 处理重复基因名
counts_matrix_clean <- counts_matrix_raw %>%
  group_by(Symbol) %>%
  slice(1) %>% 
  ungroup() %>%
  column_to_rownames("Symbol") %>%
  as.matrix()

counts_matrix <- counts_matrix_clean

# 过滤低表达基因并取整
counts_matrix <- counts_matrix[rowSums(counts_matrix) > ncol(counts_matrix), ]
counts_matrix <- round(counts_matrix)

# 3.2 读取风险分组信息 
risk_group_data <- read_csv("RiskGroups_by_median.csv") %>% 
  rename(Sample = 1, RiskGroup_Col = 2) %>% 
  mutate(RiskGroup = ifelse(grepl("high", RiskGroup_Col, ignore.case = TRUE), "high_risk", 
                            ifelse(grepl("low", RiskGroup_Col, ignore.case = TRUE), "low_risk", NA))) %>%
  filter(!is.na(RiskGroup)) %>%
  select(Sample, RiskGroup) %>%
  column_to_rownames("Sample")

# 3.3 整合数据并确保样本匹配
common_samples <- intersect(colnames(counts_matrix), rownames(risk_group_data))
counts_matrix_filtered <- counts_matrix[, common_samples]
colData <- risk_group_data[common_samples, , drop = FALSE]

# 3.4 准备 DESeq2 所需的分组数据
colData$RiskGroup <- factor(colData$RiskGroup, levels = c("low_risk", "high_risk")) 
rownames(colData) <- common_samples

# 4. 差异表达分析 (DESeq2) - 获取 GSEA 输入 ----------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts_matrix_filtered,
                              colData = colData,
                              design = ~ RiskGroup)
dds <- DESeq(dds)
res_dge <- results(dds, contrast = c("RiskGroup", "high_risk", "low_risk")) %>% 
  as.data.frame() %>% 
  rownames_to_column("Symbol")

# 5. 准备 GSEA 输入的预排序基因列表 -------------------------------------------
# 基因名转换为 Entrez ID
gene_df <- bitr(res_dge$Symbol, 
                fromType = "SYMBOL", 
                toType = "ENTREZID", 
                OrgDb = OrgDb_pkg,
                drop = TRUE) %>% 
  rename(Symbol = SYMBOL) 

# 合并 logFC 和 Entrez ID
dge_gsea_input <- res_dge %>% 
  inner_join(gene_df, by = "Symbol") %>% 
  filter(!is.na(log2FoldChange))

# 创建 GSEA 所需的命名向量: 按 log2FC 降序排列
geneList <- dge_gsea_input$log2FoldChange
names(geneList) <- dge_gsea_input$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

# 6. GSEA 富集分析 (GO 和 KEGG) ----------------------------------------------------
# 6.1 GO-BP GSEA
gsea_go <- gseGO(geneList, 
                 ont = GO_ONTOLOGY, 
                 OrgDb = OrgDb_pkg, 
                 keyType = "ENTREZID",
                 pvalueCutoff = 1, 
                 pAdjustMethod = "BH",
                 verbose = FALSE)

# 6.2 KEGG GSEA
gsea_kegg <- gseKEGG(geneList,
                     organism = "hsa", 
                     keyType = "ncbi-geneid",
                     pvalueCutoff = 1, 
                     pAdjustMethod = "BH",
                     verbose = FALSE)

# 7. 循环绘制每个基因的最显著通路 GSEA 图 ------------------------------
gsea_results_list <- list("GO" = gsea_go, "KEGG" = gsea_kegg)
target_entrez_df <- bitr(TARGET_GENES, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb_pkg)
target_map <- setNames(target_entrez_df$ENTREZID, target_entrez_df$SYMBOL)

plot_info <- list()

# 外部循环：遍历 GO 和 KEGG
N_PATHWAYS <- 6
for (ontology in ONTOLOGIES) {
  gsea_res <- gsea_results_list[[ontology]]
  N_pathways <- N_PATHWAYS # GO 和 KEGG 都是 6
  
  if (is.null(gsea_res) || nrow(gsea_res@result) == 0) {
    cat(paste("Warning: No GSEA results found for", ontology, ".\n"))
    next
  }
  
  # 提取 GSEA 基因-通路映射表 (修复：替换 GSEA.pathway.table)
  gsea_gene_map_list <- gsea_res@geneSets 
  gsea_gene_map <- data.frame(
    ID = rep(names(gsea_gene_map_list), lengths(gsea_gene_map_list)),
    ENTREZID = unlist(gsea_gene_map_list),
    stringsAsFactors = FALSE
  ) %>% mutate(ENTREZID = as.character(ENTREZID)) 
  
  # 内部循环：遍历每个目标基因
  for (gene_symbol in TARGET_GENES) {
    gene_entrez <- target_map[[gene_symbol]]
    
    if (is.null(gene_entrez) || is.na(gene_entrez)) {
      plot_info[[length(plot_info) + 1]] <- data.frame(Gene = gene_symbol, Ontology = ontology, N_Plotted = 0, Top_Pathway_ID = "N/A", Top_Pathway_Desc = "No Entrez ID", Top_Pathway_NES = NA, Top_Pathway_FDR = NA, File = "N/A")
      next
    }
    
    # 筛选包含当前基因的 Top N 通路
    top_N_pathways <- gsea_gene_map %>%
      filter(ENTREZID == gene_entrez) %>%
      inner_join(gsea_res@result, by = "ID") %>%
      # 按 p.adjust 升序和 |NES| 降序排序，选取 Top N
      arrange(p.adjust, desc(abs(NES))) %>% 
      head(N_pathways)
    
    if (nrow(top_N_pathways) > 0) {
      # 绘制通路集合
      plot_term_ids <- top_N_pathways$ID
      
      p <- gseaplot2(gsea_res, 
                     geneSet = plot_term_ids, 
                     title = paste0(gene_symbol, " ", ontology, " Pathways"),
                     # 使用 Set1 调色板，确保颜色数量足够
                     color = RColorBrewer::brewer.pal(max(3, nrow(top_N_pathways)), "Set2"),
                     base_size = 18,
                     subplots = 1:2)
      
      output_file <- paste0(gene_symbol, "_GSEA_Top", N_pathways, "_", ontology, "_Plot.pdf")
      ggsave(filename = output_file, plot = p, width = 12, height = 8,device = "pdf")
      
      # 记录摘要信息 (只记录最显著的一条)
      best_pathway <- top_N_pathways[1,]
      plot_info[[length(plot_info) + 1]] <- data.frame(
        Gene = gene_symbol,
        Ontology = ontology,
        N_Plotted = nrow(top_N_pathways),
        Top_Pathway_ID = best_pathway$ID,
        Top_Pathway_Desc = best_pathway$Description,
        Top_Pathway_NES = round(best_pathway$NES, 2),
        Top_Pathway_FDR = round(best_pathway$p.adjust, 3),
        File = output_file
      )
    } else {
      plot_info[[length(plot_info) + 1]] <- data.frame(Gene = gene_symbol, Ontology = ontology, N_Plotted = 0, Top_Pathway_ID = "N/A", Top_Pathway_Desc = "No related pathway found", Top_Pathway_NES = NA, Top_Pathway_FDR = NA, File = "N/A")
    }
  }
}

# 8. 总结绘图信息并保存到 CSV 文件 ------------------------------------------------
plot_summary_df <- bind_rows(plot_info)
write_csv(plot_summary_df, "GSEA_Top6x2_Plots_Summary.csv")
print(plot_summary_df %>% select(-File))