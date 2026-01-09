# download TCGA-LUAD 数据
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

query.exp <- GDCquery(project = "TCGA-LUAD",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
GDCdownload(query.exp)

local_dir <- "./Gene_Expression_Quantification/"
local_uuids <- list.dirs(local_dir, full.names = FALSE, recursive = FALSE)
length(local_uuids)
res <- getResults(query.exp)  # 获取原始 600 个样本的结果
res_subset <- res[res$file_id %in% local_uuids, ]

# 对应的 TCGA barcode
subset_barcodes <- unique(res_subset$cases)
length(subset_barcodes)

query_exp_236 <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR-Counts",
  barcode = subset_barcodes
)

exp_data <- GDCprepare(query = query_exp_236, directory = "/root/analysis/Raw_data/") # GDCprepare构建表达矩阵
exp_matrix <- assay(exp_data)
write.csv(exp_matrix, "LUAD_expression.csv")

clinical_query <- GDCquery(project = "TCGA-LUAD", data.category = "Clinical",
                  data.type = "Clinical Supplement", data.format = "xml")

all_dirs <- list.dirs("/root/analysis/Raw_data/GDCdata/TCGA-LUAD/Clinical/Clinical_Supplement/",
                      recursive = FALSE, full.names = TRUE)
xml_files <- list.files("/root/analysis/Raw_data/GDCdata/TCGA-LUAD/Clinical/Clinical_Supplement/",
                        pattern = "\\.xml$", recursive = TRUE, full.names = TRUE)
length(all_dirs)
length(xml_files)

dirs_missing_xml <- all_dirs[!sapply(all_dirs, function(d) {
  length(list.files(d, pattern = "\\.xml$")) > 0
})]
dirs_missing_xml

clinical_data <- GDCprepare_clinic(clinical_query, clinical.info = "patient",directory = "/root/analysis/Raw_data/GDCdata/")
expr_patient_barcodes <- substr(sample_barcodes, 1, 12)
clinical_236 <- clinical_data %>%
  filter(bcr_patient_barcode %in% expr_patient_barcodes)
write.csv(clinical_236, "LUAD_clinical.csv")

# 下载TCGA_LUAD_radiation_data
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

# 检索TCGA-LUAD的临床数据
query_clin <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "xml"
)
GDCdownload(query_clin)

## query_clin对象中移除缺失 XML 文件的目录对应的记录 ##
all_dirs <- list.dirs("/root/analysis/WGCNA/GDCdata/TCGA-LUAD/Clinical/Clinical_Supplement/",
                      recursive = FALSE, full.names = TRUE)
xml_files <- list.files("/root/analysis/WGCNA/GDCdata/TCGA-LUAD/Clinical/Clinical_Supplement/",
                        pattern = "\\.xml$", recursive = TRUE, full.names = TRUE)
dirs_missing_xml <- all_dirs[!sapply(all_dirs, function(d) {
  length(list.files(d, pattern = "\\.xml$", recursive = TRUE)) > 0
})]
file_info <- query_clin[[1]][[1]]
missing_file_ids <- basename(dirs_missing_xml)
filtered_file_info <- file_info[!file_info$file_id %in% missing_file_ids, ]
query_clin[[1]][[1]] <- filtered_file_info

## filter后可正常读入xml进行后续处理
clinical <- GDCprepare_clinic(query_clin, clinical.info = "patient",
                               directory = "/root/analysis/WGCNA/GDCdata/")
							   
# 筛选接受放疗的患者（radiation_therapy = "Yes"）
radiation_patients <- clinical %>% dplyr::filter(radiation_therapy == "YES")

# 提取所需临床信息
key_clin_info <- radiation_patients %>%
  dplyr::select(bcr_patient_barcode,
    age_at_initial_pathologic_diagnosis,gender,
    stage_event_pathologic_stage,  # 病理分期
    radiation_therapy,  # 放疗
    vital_status,days_to_death,days_to_last_followup)  # 随访时间
head(key_clin_info)
write.csv(radiation_patients, "radiation_clinical.csv",row.names = F)

# 检索对应的表达量数据（仅针对筛选出的患者）
query_exp <- GDCquery(
    project = "TCGA-LUAD",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    barcode = radiation_patients$bcr_patient_barcode)
GDCdownload(query_exp)
exp_data <- GDCprepare(query_exp)
exp_matrix <- assay(exp_data)
# 去除所有表达量都为零的行
non_zero_rows <- rowSums(exp_matrix != 0) > 0
filtered <- exp_matrix[non_zero_rows, ]
write.csv(result_filtered, "radiation_expression.csv",row.names = T)

# GSE32863 dataset processing
gpl<- read.delim("../validation/GPL6884-11607.txt", sep = "\t", header = T, 
           check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")

clin <- read.csv("../validation/GSE32863_clinical.CSV",check.names = FALSE)
gse <- read.delim("../validation/GSE32863_series_matrix.txt",sep = "\t", header = T, 
                  check.names = FALSE, stringsAsFactors = FALSE)
common_samples <- intersect(clin$Accession, colnames(gse))
cat("匹配到的共同样本数量：", length(common_samples), "\n")

gse1 <- gse[, common_samples, drop = FALSE]
gse1 <- gse1[, match(clin$Accession, colnames(gse1)), drop = FALSE]
gse1$ID_REF <- gse$ID_REF

merged_data <- merge(gse1,gpl[, c("ID", "Symbol")],,by.x="ID_REF",by.y="ID",all.x=T)
library(dplyr)
dat <- merged_data[trimws(merged_data$Symbol) != "", ] %>% group_by(Symbol) %>% 
  summarise(across(where(is.numeric), mean)) %>% as.data.frame()
rownames(dat) <- dat$Symbol
expr <- dat[, -1]
write.csv(expr,"../validation/gse32863.expMtrix.csv",row.names = T)





