# ==============================================================================
# STEP 5: COMPREHENSIVE STATISTICAL REPORT GENERATOR
#   Produces publication-ready CSV tables + Statistical_Report_Summary.txt
#   covering cohorts, 22-cell LM22 tests, 5-gene univariate/multivariable stats,
#   pairwise DeLong ROC comparisons, calibration (Brier), Cox PH, and software.
# ==============================================================================

suppressPackageStartupMessages({
  library(pROC)
  library(survival)
})

dir.create("revision_outputs/stat_report", showWarnings = FALSE, recursive = TRUE)

# Load all verified analysis outputs
load("Phase3_Immune_Phenotypes.RData")                      # harmonized_matrix, pheno_immune, immune_fractions
load("Phase5_Final_Hub_Genes.RData")                        # hub_genes
load("revision_outputs/Step2_LeakageFree_Results.RData")    # Nested CV, End-to-End CV, LOCO
load("revision_outputs/Step3_TCGA_External_Results.RData")  # TCGA External ROC, Biomarkers, Cox
load("revision_outputs/Step4_Sensitivity_Results.RData")    # Sensitivity Analysis

y_geo_fac <- factor(pheno_immune$Immune_Phenotype, levels = c("Immune-Cold", "Immune-Hot"))
y_geo_bin <- as.integer(y_geo_fac == "Immune-Hot")

# Helper: Rank-biserial correlation effect size for Wilcoxon rank-sum test
rank_biserial <- function(x_hot, x_cold) {
  n1 <- length(x_hot); n0 <- length(x_cold)
  u  <- unname(wilcox.test(x_hot, x_cold, exact = FALSE)$statistic)
  (2 * u) / (n1 * n0) - 1  # positive => higher in Hot
}

# Helper: Manual VIF from design matrix (no extra package needed)
calc_vif <- function(X_df) {
  sapply(colnames(X_df), function(col) {
    f <- as.formula(paste(col, "~ ."))
    r2 <- summary(lm(f, data = X_df))$r.squared
    1 / (1 - r2)
  })
}

# Helper: Brier Score
brier_score <- function(y, p) mean((p - y)^2)

# ==============================================================================
# 1. COHORT CHARACTERISTICS & BATCH INDEPENDENCE
# ==============================================================================
tab_geo  <- table(pheno_immune$Batch, y_geo_fac)
chi_geo  <- chisq.test(tab_geo)
cramer_v <- sqrt(unname(chi_geo$statistic) / (sum(tab_geo) * (min(dim(tab_geo)) - 1)))

cohort_df <- data.frame(
  Cohort      = c(rownames(tab_geo), "GEO Discovery (Pooled)", "TCGA-LUAD (External Classification)", "TCGA-LUAD (Complete Survival + Stage)"),
  Platform    = c(rep("Affymetrix GPL570", 5), rep("Illumina RNA-seq (GDC STAR-TPM)", 2)),
  N_Total     = c(rowSums(tab_geo), sum(tab_geo), length(tcga_pheno), nrow(cox_data)),
  N_Hot       = c(tab_geo[, "Immune-Hot"],  sum(tab_geo[, "Immune-Hot"]),  sum(tcga_pheno == "Immune-Hot"),  sum(cox_data$Risk_Group == "Immune-Hot-like (Low Risk)")),
  N_Cold      = c(tab_geo[, "Immune-Cold"], sum(tab_geo[, "Immune-Cold"]), sum(tcga_pheno == "Immune-Cold"), sum(cox_data$Risk_Group == "Immune-Cold-like (High Risk)")),
  stringsAsFactors = FALSE
)
cohort_df$Pct_Hot  <- round(100 * cohort_df$N_Hot  / cohort_df$N_Total, 1)
cohort_df$Pct_Cold <- round(100 * cohort_df$N_Cold / cohort_df$N_Total, 1)

write.csv(cohort_df, "revision_outputs/stat_report/Table_Step5_1_Cohort_Characteristics.csv", row.names = FALSE)

# ==============================================================================
# 2. FULL 22-CELL LM22 DIFFERENTIAL INFILTRATION REPORT (GEO & TCGA)
# ==============================================================================
summarize_lm22 <- function(frac_mat, pheno_fac, cohort_lbl) {
  hot_idx  <- which(pheno_fac == "Immune-Hot")
  cold_idx <- which(pheno_fac == "Immune-Cold")
  rows <- lapply(colnames(frac_mat), function(cell) {
    xh <- frac_mat[hot_idx, cell]
    xc <- frac_mat[cold_idx, cell]
    wt <- wilcox.test(xh, xc, exact = FALSE)
    data.frame(
      Cohort         = cohort_lbl,
      Cell_Type      = cell,
      Hot_Mean_SD    = sprintf("%.4f ± %.4f", mean(xh), sd(xh)),
      Cold_Mean_SD   = sprintf("%.4f ± %.4f", mean(xc), sd(xc)),
      Hot_Median_IQR = sprintf("%.4f [%.4f-%.4f]", median(xh), quantile(xh, 0.25), quantile(xh, 0.75)),
      Cold_Median_IQR= sprintf("%.4f [%.4f-%.4f]", median(xc), quantile(xc, 0.25), quantile(xc, 0.75)),
      Diff_Mean      = round(mean(xh) - mean(xc), 4),
      Rank_Biserial_r= round(rank_biserial(xh, xc), 4),
      Wilcoxon_P     = signif(wt$p.value, 3),
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  df$BH_FDR <- signif(p.adjust(df$Wilcoxon_P, method = "BH"), 3)
  df[order(df$BH_FDR, -abs(df$Rank_Biserial_r)), ]
}

lm22_geo_df  <- summarize_lm22(immune_fractions, y_geo_fac,  "GEO_Discovery_n536")
lm22_tcga_df <- summarize_lm22(tcga_lm22,        tcga_pheno, "TCGA_External_n516")
lm22_all_df  <- rbind(lm22_geo_df, lm22_tcga_df)
write.csv(lm22_all_df, "revision_outputs/stat_report/Table_Step5_2_LM22_Differential_GEO_TCGA.csv", row.names = FALSE)

# ==============================================================================
# 3. FIVE-GENE UNIVARIATE, MULTIVARIABLE GLM, VIF & CALIBRATION REPORT
# ==============================================================================
X_geo_raw <- t(harmonized_matrix[hub_genes, , drop = FALSE])
X_geo_z   <- scale(X_geo_raw)

# Fit standardized multivariable GLM predicting Immune-Hot = 1 (matches locked_glm)
df_glm_std <- data.frame(X_geo_z, y = y_geo_bin, check.names = TRUE)
glm_std    <- glm(y ~ ., data = df_glm_std, family = binomial)
glm_sum    <- summary(glm_std)$coefficients
glm_ci     <- suppressMessages(confint.default(glm_std))
vifs       <- calc_vif(as.data.frame(X_geo_z, check.names = TRUE))

gene_rows <- lapply(hub_genes, function(g) {
  xh <- X_geo_raw[y_geo_bin == 1, g]
  xc <- X_geo_raw[y_geo_bin == 0, g]
  wt <- wilcox.test(xh, xc, exact = FALSE)
  ro <- roc(y_geo_bin, X_geo_raw[, g], levels = c(0, 1), direction = "<", ci = TRUE, quiet = TRUE)
  
  beta <- glm_sum[g, "Estimate"]
  se   <- glm_sum[g, "Std. Error"]
  zval <- glm_sum[g, "z value"]
  pval <- glm_sum[g, "Pr(>|z|)"]
  or   <- exp(beta)
  or_lo<- exp(glm_ci[g, 1])
  or_hi<- exp(glm_ci[g, 2])
  
  data.frame(
    Gene            = g,
    Hot_Mean_SD     = sprintf("%.3f ± %.3f", mean(xh), sd(xh)),
    Cold_Mean_SD    = sprintf("%.3f ± %.3f", mean(xc), sd(xc)),
    Delta_Log2      = round(mean(xh) - mean(xc), 3),
    Rank_Biserial_r = round(rank_biserial(xh, xc), 4),
    Univar_Wilcox_P = signif(wt$p.value, 3),
    SingleGene_AUC  = sprintf("%.4f [%.4f-%.4f]", auc(ro), ro$ci[1], ro$ci[3]),
    Std_Beta        = round(beta, 4),
    Std_SE          = round(se, 4),
    OR_per_SD       = sprintf("%.3f [%.3f-%.3f]", or, or_lo, or_hi),
    Wald_Z          = round(zval, 3),
    Multivar_P      = signif(pval, 3),
    VIF             = round(unname(vifs[g]), 3),
    stringsAsFactors = FALSE
  )
})
gene_stat_df <- do.call(rbind, gene_rows)
gene_stat_df$Univar_BH_FDR <- signif(p.adjust(gene_stat_df$Univar_Wilcox_P, method = "BH"), 3)
write.csv(gene_stat_df, "revision_outputs/stat_report/Table_Step5_3_FiveGene_Stats_and_GLM.csv", row.names = FALSE)

# Calibration (Brier score) & Overall Model Goodness-of-Fit
null_brier_geo  <- brier_score(y_geo_bin, rep(mean(y_geo_bin), length(y_geo_bin)))
brier_geo_oof   <- c(
  GLM = brier_score(y_geo_bin, oof_glm),
  RF  = brier_score(y_geo_bin, oof_rf),
  SVM = brier_score(y_geo_bin, oof_svm),
  ANN = brier_score(y_geo_bin, oof_ann)
)
brier_loco <- c(
  GLM = brier_score(y_geo_bin, loco_oof_glm),
  RF  = brier_score(y_geo_bin, loco_oof_rf),
  SVM = brier_score(y_geo_bin, loco_oof_svm),
  ANN = brier_score(y_geo_bin, loco_oof_ann)
)
mcfadden_r2 <- 1 - (glm_std$deviance / glm_std$null.deviance)

# ==============================================================================
# 4. MASTER BENCHMARKING TABLE + PAIRWISE DELONG ROC TESTS
# ==============================================================================
master_bench <- rbind(partA_table, partB_table, loco_pooled, partA_tcga_table)
write.csv(master_bench, "revision_outputs/stat_report/Table_Step5_4_Master_Model_Benchmark.csv", row.names = FALSE)

run_delong <- function(r1, r2, name1, name2, context) {
  dt <- roc.test(r1, r2, method = "delong", paired = TRUE)
  data.frame(
    Context   = context,
    Model_1   = name1,
    Model_2   = name2,
    AUC_1     = round(as.numeric(auc(r1)), 4),
    AUC_2     = round(as.numeric(auc(r2)), 4),
    Delta_AUC = round(as.numeric(auc(r1)) - as.numeric(auc(r2)), 4),
    Z_Stat    = round(as.numeric(dt$statistic), 3),
    DeLong_P  = signif(dt$p.value, 3),
    stringsAsFactors = FALSE
  )
}

delong_df <- rbind(
  # Nested 10-Fold OOF pairwise comparisons (GEO)
  run_delong(res_glm_app$roc, res_glm_oof$roc, "GLM (Apparent)", "GLM (Nested OOF)", "GEO Apparent vs OOF"),
  run_delong(res_glm_oof$roc, res_svm_oof$roc, "GLM", "SVM-Radial",  "GEO Nested 10-Fold OOF"),
  run_delong(res_glm_oof$roc, res_rf_oof$roc,  "GLM", "RF",          "GEO Nested 10-Fold OOF"),
  run_delong(res_glm_oof$roc, res_ann_oof$roc, "GLM", "ANN",         "GEO Nested 10-Fold OOF"),
  run_delong(res_svm_oof$roc, res_rf_oof$roc,  "SVM-Radial", "RF",   "GEO Nested 10-Fold OOF"),
  run_delong(res_svm_oof$roc, res_ann_oof$roc, "SVM-Radial", "ANN",  "GEO Nested 10-Fold OOF"),
  # External TCGA-LUAD pairwise comparisons
  run_delong(ext_glm$roc, ext_svm$roc, "GLM", "SVM-Radial", "TCGA-LUAD External"),
  run_delong(ext_glm$roc, ext_rf$roc,  "GLM", "RF",         "TCGA-LUAD External"),
  run_delong(ext_glm$roc, ext_ann$roc, "GLM", "ANN",        "TCGA-LUAD External")
)
write.csv(delong_df, "revision_outputs/stat_report/Table_Step5_4b_Pairwise_DeLong_Tests.csv", row.names = FALSE)

# ==============================================================================
# 5. SOFTWARE & PACKAGE REPRODUCIBILITY RECORD
# ==============================================================================
pkg_list <- c("GEOquery", "affy", "frma", "hgu133plus2.db", "sva", "WGCNA",
              "glmnet", "randomForest", "e1071", "keras", "reticulate",
              "pROC", "survival", "survminer", "cluster", "ggplot2", "patchwork",
              "pheatmap", "data.table")
ver_df <- data.frame(
  Package = c("R", pkg_list),
  Version = c(
    paste0(R.version$major, ".", R.version$minor),
    sapply(pkg_list, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "Not Installed")
  ),
  stringsAsFactors = FALSE
)
write.csv(ver_df, "revision_outputs/stat_report/Table_Step5_5_Software_Versions.csv", row.names = FALSE)

# ==============================================================================
# 6. WRITE HUMAN-READABLE MASTER STATISTICAL REPORT
# ==============================================================================
rep_file <- "revision_outputs/stat_report/Statistical_Report_Summary.txt"
sink(rep_file)
cat("==============================================================================\n")
cat("COMPREHENSIVE STATISTICAL REPORT — LUAD 5-GENE IMMUNE PHENOTYPING STUDY\n")
cat("Locked Consensus Signature: ", paste(hub_genes, collapse = ", "), "\n", sep = "")
cat("==============================================================================\n\n")

cat("1. COHORT CHARACTERISTICS & BATCH INDEPENDENCE\n")
print(cohort_df, row.names = FALSE)
cat(sprintf("\nGEO Cohort x Phenotype Chi-square: X2 = %.3f, df = %d, p = %.4f, Cramer's V = %.4f\n\n",
            chi_geo$statistic, chi_geo$parameter, chi_geo$p.value, cramer_v))

cat("2. TOP 10 DIFFERENTIALLY INFILTRATED LM22 CELL TYPES (GEO DISCOVERY, n = 536)\n")
print(head(lm22_geo_df[, c("Cell_Type", "Hot_Median_IQR", "Cold_Median_IQR", "Rank_Biserial_r", "BH_FDR")], 10), row.names = FALSE)

cat("\n3. FIVE-GENE DIFFERENTIAL EXPRESSION, MULTIVARIABLE GLM & VIF (GEO, n = 536)\n")
print(gene_stat_df[, c("Gene", "Delta_Log2", "SingleGene_AUC", "OR_per_SD", "Multivar_P", "VIF", "Univar_BH_FDR")], row.names = FALSE)
cat(sprintf("\nMultivariable GLM McFadden Pseudo-R2 = %.4f | Likelihood Ratio p = %.3e\n",
            mcfadden_r2, pchisq(glm_std$null.deviance - glm_std$deviance, df = 5, lower.tail = FALSE)))
cat(sprintf("GEO Calibration Brier Scores (Null = %.4f): GLM=%.4f, RF=%.4f, SVM=%.4f, ANN=%.4f\n",
            null_brier_geo, brier_geo_oof["GLM"], brier_geo_oof["RF"], brier_geo_oof["SVM"], brier_geo_oof["ANN"]))
cat(sprintf("LOCO Pre-ComBat Brier Scores:              GLM=%.4f, RF=%.4f, SVM=%.4f, ANN=%.4f\n\n",
            brier_loco["GLM"], brier_loco["RF"], brier_loco["SVM"], brier_loco["ANN"]))

cat("4. MASTER CLASSIFICATION BENCHMARKING TABLE\n")
print(master_bench[, c("Model", "Evaluation", "AUC", "CI_Lower", "CI_Upper", "Balanced_Accuracy", "Sensitivity", "Specificity", "F1_Score")], row.names = FALSE)

cat("\n5. PAIRWISE PAIRED DELONG ROC TESTS\n")
print(delong_df, row.names = FALSE)

cat("\n6. TCGA-LUAD SURVIVAL & STAGE STATISTICS (FULL 5-GENE PANEL)\n")
print(cox_summary_table, row.names = FALSE)
cat(sprintf("\nMultivariate Cox Concordance (C-index) = %.4f (SE = %.4f)\n",
            summary(cox_multi)$concordance[1], summary(cox_multi)$concordance[2]))
cat(sprintf("Schoenfeld PH Global Test: chi2 = %.3f, df = %d, p = %.4f\n",
            zph_multi$table["GLOBAL", "chisq"], zph_multi$table["GLOBAL", "df"], zph_multi$table["GLOBAL", "p"]))
cat(sprintf("Pathologic Stage I-IV Kruskal-Wallis (n = 508): chi2 = %.3f, p = %.4g\n\n",
            kw_stage$statistic, kw_stage$p.value))

cat("7. SOFTWARE & PACKAGE VERSIONS\n")
print(ver_df, row.names = FALSE)
sink()

# Print concise output to console
cat(readLines(rep_file), sep = "\n")
cat("\nAll tables saved in revision_outputs/stat_report/\n")
