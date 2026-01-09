#############################################
# LUAD_Radiosensitivity_Signature.R
# 完整管道：单因素Cox -> LASSO-Cox -> 风险评分 -> KM/ROC -> 多变量Cox(独立预后)
# 输入：
#   - exp_subset: 行=基因 (gene_symbol), 列=样本 (TCGA barcodes)
#   - exp_clinical: 样本为行，包含 bcr_patient_barcode, OS_time, OS_status
#############################################

# 依赖包
required <- c("survival","survminer","glmnet","timeROC","dplyr","ggplot2","pheatmap","rms","pec")
install_if_missing <- function(pkgs){
  for(pkg in pkgs){
    if(!requireNamespace(pkg, quietly = TRUE)){
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}
install_if_missing(required)

# ---------------------------
# 0. 参数配置
# ---------------------------
time_points <- c(365, 3*365, 5*365)   # 1y, 3y, 5y (days)
lasso_nfolds <- 10
topN_fallback <- 20

# ---------------------------
# 1. 检查并准备输入数据
# ---------------------------
clin <- read.csv("~/analysis/Raw_data/GDCdata/TCGA-LUAD/Clinical/LUAD_clinical.csv",row.names = 1,check.names = F)
exp <- read.csv("~/analysis/Model_construction/Radiosensitivity_genes_expMatrix.csv",row.names = 1,check.names = FALSE)
# 1.1 处理exp
key1 <- read.csv("~/analysis/WGCNA/vital_module.csv",row.names = 1)
key2 <- read.csv("~/analysis/WGCNA/OStime_module.csv",row.names = 1)
merged_keys <- c(rownames(key1), rownames(key2))
radiogene_exp <- as.data.frame(exp[merged_keys, , drop = FALSE])
rownames(radiogene_exp) <- merged_keys
# 1.2处理clin
exp_clinical <- patient_order %>% left_join(clin, by = c("patient_id" = "bcr_patient_barcode"))
exp_clinical <- exp_clinical %>%
  mutate(OS_time = ifelse(!is.na(days_to_death), days_to_death, days_to_last_followup),
         OS_status = ifelse(vital_status == "Dead", 1, 0))
rownames(exp_clinical) <- colnames(exp)  # exp和clin对齐
# 将表达矩阵按基因行，样本列，样本顺序与临床一致exp <- exp[, exp_clinical[[ifelse("bcr_patient_barcode" %in% colnames(exp_clinical),"bcr_patient_barcode","sample")]], drop = FALSE]

# ---------------------------
# 2. 单因素 Cox 回归
# ---------------------------
library(survival)
genes <- rownames(exp)
cox_list <- lapply(genes, function(g){
  vec <- as.numeric(exp[g, ])
  fit <- tryCatch(coxph(Surv(exp_clinical$OS_time, exp_clinical$OS_status) ~ vec),
                  error = function(e) NULL)
  if(is.null(fit)) return(data.frame(gene = g, coef = NA, HR = NA, p = NA))
  s <- summary(fit)
  coef <- s$coefficients[1,"coef"]
  hr <- s$coefficients[1,"exp(coef)"]
  pval <- s$coefficients[1,"Pr(>|z|)"]
  data.frame(gene = g, coef = coef, HR = hr, p = pval, stringsAsFactors = FALSE)
})
cox_df <- do.call(rbind, cox_list)
write.csv(cox_df, "TCGA_LUAD_univariate_cox_all_genes.csv", row.names = FALSE)
sig_genes <- na.omit(cox_df) %>% dplyr::filter(p < 0.05) %>% dplyr::arrange(p)
write.csv(sig_genes, "Cox_SignificantGenes.csv", row.names = FALSE)

# ---------------------------
# 3. LASSO-Cox
# ---------------------------
library(glmnet)
candidate_genes <- sig_genes$gene
x_mat <- t(exp[candidate_genes, , drop = FALSE])
x_mat[is.na(x_mat)] <- 0
y_surv <- Surv(exp_clinical$OS_time, exp_clinical$OS_status)
lasso_nfolds <- 10
set.seed(1234)
cvfit <- cv.glmnet(x_mat, y_surv, family = "cox", alpha = 1, nfolds = lasso_nfolds)
pdf("LASSO_CVplot.pdf", width = 6, height = 5)
plot(cvfit)
dev.off()
coef_min <- coef(cvfit, s = "lambda.min")
selected_idx <- which(as.numeric(coef_min) != 0)
selected_genes = c("VEGFC", "SERPINE1", "CHEK1","MEG3","SLC16A1","PRRX1")
coef_final <- as.numeric(coef_min[selected_genes, , drop = FALSE])
names(coef_final) <- selected_genes


# ---------------------------
# 4. 风险评分计算与分组
# ---------------------------
expr_for_score <- t(exp[selected_genes, , drop = FALSE])  # 样本 x gene
risk_score <- as.numeric(expr_for_score %*% coef_final)
risk_score <- as.vector(risk_score)
# 将风险得分加入 exp_clinical
exp_clinical$risk_score <- risk_score
exp_clinical$risk_group <- ifelse(exp_clinical$risk_score >= median(exp_clinical$risk_score, na.rm = TRUE), "High", "Low")
write.csv(exp_clinical, "RiskScore_and_clinical.csv", row.names = FALSE)

# ---------------------------
# 5. Kaplan-Meier 曲线（High vs Low）
# ---------------------------
library(survminer)
fit_km <- survfit(Surv(OS_time, OS_status) ~ risk_group, data = exp_clinical)
km_plot <- ggsurvplot(fit_km,
                      data = exp_clinical,
                      pval = TRUE, conf.int = FALSE,
                      risk.table = TRUE, risk.table.col = "strata",
                      title = "KM: Risk group (LASSO signature)",
                      legend.labs = c("High","Low"),
                      xlab = "Time (days)")
ggsave("KM_riskgroup.png", plot = km_plot$plot, width = 7, height = 6)
ggsave("KM_riskgroup_risktable.png", plot = km_plot$table, width = 7, height = 2)

# ---------------------------
# 6. timeROC
# ---------------------------
library(timeROC)
good_idx <- which(!is.na(exp_clinical$risk_score) & !is.na(exp_clinical$OS_time) & exp_clinical$OS_time > 0)
tr <- exp_clinical$OS_time[good_idx]
ev <- exp_clinical$OS_status[good_idx]
rs <- exp_clinical$risk_score[good_idx]

time_roc_obj <- timeROC(T = tr, delta = ev, marker = rs, cause = 1, weighting = "marginal", times = time_points, iid = TRUE)
auc_df <- data.frame(time = time_points, AUC = time_roc_obj$AUC)
write.csv(auc_df, "timeROC_AUCs.csv", row.names = FALSE)
# 绘制 ROC 曲线
png("TimeROC_plot.png", width = 800, height = 600)
plot(time_roc_obj, time = time_points[1], col = 1, add = FALSE)
cols <- 1:length(time_points)
for(i in seq_along(time_points)){
  plot(time_roc_obj, time = time_points[i], add = (i!=1), col = cols[i])
}
legend("bottomright", legend = paste0(c("1y","3y","5y")," AUC=", round(time_roc_obj$AUC,3)), col = cols, lwd = 2)
dev.off()

# ---------------------------
# 7. 多变量 Cox（独立预后分析）
# ---------------------------
clinical_vars <- c("age_at_initial_pathologic_diagnosis","gender","stage_event_clinical_stage")
clinical_vars_exist <- clinical_vars[clinical_vars %in% colnames(exp_clinical)]
if(length(clinical_vars_exist) == 0){
  message("没有找到预设的临床协变量，multi-Cox 只使用 risk_score。")
  form_str <- "Surv(OS_time, OS_status) ~ risk_score"
} else {
  form_str <- paste0("Surv(OS_time, OS_status) ~ risk_score + ", paste(clinical_vars_exist, collapse = " + "))
}
multi_cox <- coxph(as.formula(form_str), data = exp_clinical)
multi_sum <- summary(multi_cox)
write.csv(as.data.frame(multi_sum$coefficients), "multivariate_Cox_results.csv")

lasso_coef_df <- data.frame(gene = names(coef_final), coef = as.numeric(coef_final), stringsAsFactors = FALSE)
write.csv(lasso_coef_df, "LASSO_final_coefficients.csv", row.names = FALSE)
