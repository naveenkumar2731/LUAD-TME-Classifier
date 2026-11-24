# ==============================================================================
# 10_generate_publication_figures.R
# Unified 600 DPI publication-grade figure generator:
#   1. Figure_Rev1_PCA_Batch_Correction.png
#   2. Figure_Rev2_LM22_Immune_Landscape.png
#   3. Figure_Rev3_WGCNA_Network.png
#   4. Figure_Rev4_TriML_Feature_Selection.png
#   5. Figure_Step2A_NestedCV_ROC.png
#   6. Figure_Step2B_ANN_Training_Curves.png
#   7. Figure_Step2C_EndToEnd_NestedCV_Stability.png
#   8. Figure_Step2D_LOCO_Validation.png
# (Figures 3A-3C are produced by 07_tcga_external_validation.R and Figure 4
#  is produced by 08_threshold_sensitivity.R.)
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(WGCNA)
  library(glmnet)
  library(randomForest)
  library(e1071)
  library(pROC)
})

dir.create("revision_outputs/figures", showWarnings = FALSE, recursive = TRUE)

theme_pub <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title        = element_text(face = "bold", size = base_size + 1, hjust = 0.5),
      plot.subtitle     = element_text(size = base_size - 1, hjust = 0.5, color = "gray30"),
      axis.title        = element_text(face = "bold", size = base_size),
      axis.text         = element_text(size = base_size - 1, color = "black"),
      legend.title      = element_text(face = "bold", size = base_size - 1),
      legend.text       = element_text(size = base_size - 2),
      legend.background = element_rect(fill = ggplot2::alpha("white", 0.90), color = "gray80"),
      plot.margin       = ggplot2::margin(10, 12, 10, 10)
    )
}

roc_to_df <- function(roc_obj, label) {
  df <- data.frame(
    FPR   = 1 - roc_obj$specificities,
    TPR   = roc_obj$sensitivities,
    Model = label,
    stringsAsFactors = FALSE
  )
  df[order(df$FPR, df$TPR), ]
}

load("Phase3_Immune_Phenotypes.RData")
load("Phase4_WGCNA_Input.RData")
load("Phase4_WGCNA_Modules_Corrected.RData")
load("Phase5_Final_Hub_Genes.RData")
load("revision_outputs/Step2_LeakageFree_Results.RData")

y_bin       <- as.integer(pheno_immune$Immune_Phenotype == "Immune-Hot")
cohort_cols <- c("GSE30219" = "#E64B35FF", "GSE31210" = "#4DBBD5FF",
                 "GSE37745" = "#00A087FF", "GSE50081" = "#3C5488FF")

# ==============================================================================
# 1. PCA BEFORE vs. AFTER COMBAT HARMONIZATION
# ==============================================================================
load("Phase2_Mapped_Data_and_Metadata.RData")
raw_sub <- mapped_matrix[rownames(harmonized_matrix), pheno_immune$Sample, drop = FALSE]
rm(mapped_matrix, meta_31210, meta_50081, meta_37745, meta_30219)
gc()

rv_raw  <- apply(raw_sub, 1, var, na.rm = TRUE)
raw_sub <- raw_sub[is.finite(rv_raw) & rv_raw > 0, , drop = FALSE]

pca_bef <- prcomp(t(raw_sub), scale. = TRUE)
var_bef <- round(100 * summary(pca_bef)$importance[2, 1:2], 1)
df_bef  <- data.frame(PC1 = pca_bef$x[, 1], PC2 = pca_bef$x[, 2], Batch = pheno_immune$Batch)

pca_aft <- prcomp(t(harmonized_matrix), scale. = TRUE)
var_aft <- round(100 * summary(pca_aft)$importance[2, 1:2], 1)
df_aft  <- data.frame(PC1 = pca_aft$x[, 1], PC2 = pca_aft$x[, 2], Batch = pheno_immune$Batch)

p_pca1 <- ggplot(df_bef, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 2.2, alpha = 0.72) +
  scale_color_manual(values = cohort_cols) +
  labs(title = "A. Before ComBat (Pre-Harmonization fRMA)",
       x = sprintf("PC1 (%.1f%% variance)", var_bef[1]),
       y = sprintf("PC2 (%.1f%% variance)", var_bef[2]),
       color = "GEO Cohort") +
  theme_pub(14)

p_pca2 <- ggplot(df_aft, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 2.2, alpha = 0.72) +
  scale_color_manual(values = cohort_cols) +
  labs(title = "B. After ComBat (Harmonized LUAD Cohort, n = 536)",
       x = sprintf("PC1 (%.1f%% variance)", var_aft[1]),
       y = sprintf("PC2 (%.1f%% variance)", var_aft[2]),
       color = "GEO Cohort") +
  theme_pub(14)

p_rev1 <- (p_pca1 | p_pca2) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.text = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"))

ggsave("revision_outputs/figures/Figure_Rev1_PCA_Batch_Correction.png",
       plot = p_rev1, width = 12.8, height = 5.6, dpi = 600)

# ==============================================================================
# 2. LM22 IMMUNE INFILTRATION LANDSCAPE HEATMAP (Z-SCORE CAPPED ±3 SD)
# ==============================================================================
mat_z     <- t(scale(immune_fractions))
mat_z[!is.finite(mat_z)] <- 0
mat_z_cap <- pmax(pmin(mat_z, 3), -3)

hc_cols <- hclust(dist(scale(immune_fractions), method = "euclidean"), method = "ward.D2")
hc_rows <- hclust(dist(mat_z, method = "euclidean"), method = "ward.D2")

ann_col <- data.frame(
  Phenotype = factor(pheno_immune$Immune_Phenotype, levels = c("Immune-Hot", "Immune-Cold")),
  Cohort    = factor(pheno_immune$Batch)
)
rownames(ann_col) <- pheno_immune$Sample

ann_colors <- list(
  Phenotype = c("Immune-Hot" = "#E64B35FF", "Immune-Cold" = "#3C5488FF"),
  Cohort    = cohort_cols
)

png("revision_outputs/figures/Figure_Rev2_LM22_Immune_Landscape.png",
    width = 13.0, height = 8.0, units = "in", res = 600)
pheatmap(
  mat_z_cap,
  cluster_cols      = hc_cols,
  cluster_rows      = hc_rows,
  annotation_col    = ann_col,
  annotation_colors = ann_colors,
  show_colnames     = FALSE,
  fontsize          = 13,
  fontsize_row      = 12.5,
  treeheight_col    = 45,
  treeheight_row    = 35,
  breaks            = seq(-3, 3, length.out = 101),
  color             = colorRampPalette(c("#3C5488FF", "#F7F7F7", "#E64B35FF"))(100),
  main              = "LM22 Immune Infiltration Landscape (n = 536, Z-Score ±3 SD)"
)
dev.off()

# ==============================================================================
# 3. SIGNED WGCNA SCALE-FREE TOPOLOGY & MODULE-TRAIT CORRELATIONS
# ==============================================================================
sft_df <- data.frame(
  Power     = sft$fitIndices[, 1],
  Signed_R2 = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  Mean_Conn = sft$fitIndices[, 5]
)

p_w1 <- ggplot(sft_df, aes(x = Power, y = Signed_R2)) +
  geom_hline(yintercept = 0.85, linetype = "dashed", color = "#3C5488FF", linewidth = 0.9) +
  geom_point(aes(color = Power == 9), size = 3.2, show.legend = FALSE) +
  geom_text(aes(label = Power, color = Power == 9), vjust = -0.85, size = 4.2, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c("FALSE" = "gray30", "TRUE" = "#E64B35FF")) +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(title = "A. Scale-Free Topology Fit",
       subtitle = "Optimal soft threshold β = 9 (R² ≥ 0.85)",
       x = "Soft-Thresholding Power (β)", y = "Signed R²") +
  theme_pub(13)

p_w2 <- ggplot(sft_df, aes(x = Power, y = Mean_Conn)) +
  geom_line(color = "gray65", linewidth = 0.8) +
  geom_point(aes(color = Power == 9), size = 3.2, show.legend = FALSE) +
  geom_text(aes(label = Power, color = Power == 9), vjust = -0.85, size = 4.0, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c("FALSE" = "gray30", "TRUE" = "#E64B35FF")) +
  coord_cartesian(ylim = c(0, max(sft_df$Mean_Conn) * 1.12)) +
  labs(title = "B. Mean Network Connectivity",
       subtitle = "Preserved connectivity at β = 9",
       x = "Soft-Thresholding Power (β)", y = "Mean Connectivity") +
  theme_pub(13)

n_samp_w  <- nrow(datExpr_filtered)
mod_pvals <- WGCNA::corPvalueStudent(moduleTraitCor_clean, n_samp_w)
mod_df <- data.frame(
  Module = rownames(moduleTraitCor_clean),
  Cor    = as.numeric(moduleTraitCor_clean[, 1]),
  Pval   = as.numeric(mod_pvals[, 1]),
  stringsAsFactors = FALSE
)
mod_df        <- mod_df[order(mod_df$Cor), ]
mod_df$Module <- factor(mod_df$Module, levels = mod_df$Module)
mod_df$Label  <- sprintf("r = %.2f  (p = %.1e)", mod_df$Cor, mod_df$Pval)

p_w3 <- ggplot(mod_df, aes(x = Cor, y = Module, fill = Cor)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  geom_col(width = 0.72, color = "black", linewidth = 0.35) +
  geom_text(aes(label = Label, hjust = ifelse(Cor < 0, 1.05, -0.05)), size = 3.8, fontface = "bold") +
  scale_fill_gradient2(low = "#3C5488FF", mid = "white", high = "#E64B35FF", midpoint = 0, limits = c(-0.7, 0.7)) +
  coord_cartesian(xlim = c(-1.10, 0.85)) +
  labs(title = "C. WGCNA Module-Trait Correlations (Immune-Cold = 1)",
       subtitle = "MEmagenta (134 genes) shows strongest immune association (r = -0.58, p = 9e-50)",
       x = "Pearson Correlation with Immune-Cold Phenotype", y = NULL) +
  theme_pub(13) +
  theme(legend.position = "none", axis.text.y = element_text(face = "bold", size = 12))

p_rev3 <- ((p_w1 / p_w2) | p_w3) + plot_layout(widths = c(1, 1.35))
ggsave("revision_outputs/figures/Figure_Rev3_WGCNA_Network.png",
       plot = p_rev3, width = 14.0, height = 7.6, dpi = 600)

# ==============================================================================
# 4. TRI-ML CONSENSUS FEATURE SELECTION (4-PANEL)
# ==============================================================================
X_mod  <- t(harmonized_matrix[target_genes, , drop = FALSE])
y_fac  <- factor(pheno_immune$Immune_Phenotype, levels = c("Immune-Cold", "Immune-Hot"))
y_cold <- ifelse(y_fac == "Immune-Cold", 1, 0)

set.seed(123)
cv_lasso <- cv.glmnet(x = as.matrix(X_mod), y = y_cold, family = "binomial", alpha = 1, nfolds = 10)
lasso_df <- data.frame(
  LogLam = log(cv_lasso$lambda),
  CVM    = cv_lasso$cvm,
  CVLo   = cv_lasso$cvlo,
  CVUp   = cv_lasso$cvup
)

p_m1 <- ggplot(lasso_df, aes(x = LogLam, y = CVM)) +
  geom_ribbon(aes(ymin = CVLo, ymax = CVUp), fill = "#3C5488FF", alpha = 0.16) +
  geom_line(color = "#3C5488FF", linewidth = 1.1) +
  geom_point(color = "#E64B35FF", size = 1.6) +
  geom_vline(xintercept = log(cv_lasso$lambda.min), linetype = "dashed", color = "#E64B35FF", linewidth = 0.9) +
  geom_vline(xintercept = log(cv_lasso$lambda.1se), linetype = "dotted", color = "gray40", linewidth = 0.8) +
  annotate("text", x = log(cv_lasso$lambda.min) - 0.18, y = max(lasso_df$CVM) * 0.90,
           label = sprintf("Optimal λmin\n(%d genes)", length(lasso_genes)),
           hjust = 1, size = 3.9, fontface = "bold", color = "#E64B35FF") +
  labs(title = "A. LASSO 10-Fold CV Selection", x = "Log(λ)", y = "Binomial Deviance") +
  theme_pub(13)

set.seed(123)
rf_fs  <- randomForest(x = X_mod, y = y_fac, ntree = 1000, importance = TRUE)
rf_imp <- sort(importance(rf_fs)[, "MeanDecreaseGini"], decreasing = TRUE)[1:25]
rf_df  <- data.frame(
  Gene  = factor(names(rf_imp), levels = rev(names(rf_imp))),
  Gini  = as.numeric(rf_imp),
  IsHub = names(rf_imp) %in% hub_genes
)

p_m2 <- ggplot(rf_df, aes(x = Gini, y = Gene, fill = IsHub)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.28) +
  scale_fill_manual(values = c("FALSE" = "#9BB7D4", "TRUE" = "#E64B35FF"),
                    labels = c("FALSE" = "Candidate", "TRUE" = "Consensus Hub (5/5)")) +
  labs(title = "B. Random Forest Top-25 Gini Importance",
       x = "Mean Decrease in Gini Impurity", y = NULL, fill = NULL) +
  theme_pub(13) +
  theme(legend.position = "inside", legend.position.inside = c(0.72, 0.16),
        axis.text.y = element_text(face = "bold.italic", size = 9.5))

set.seed(123)
kf_rfe  <- 5L
f_rfe   <- sample(rep(1:kf_rfe, length.out = nrow(X_mod)))
acc_rfe <- numeric(length(svm_genes))
for (i in seq_along(svm_genes)) {
  cur <- svm_genes[1:i]
  fa  <- numeric(kf_rfe)
  for (f in seq_len(kf_rfe)) {
    m <- svm(x = X_mod[f_rfe != f, cur, drop = FALSE], y = y_fac[f_rfe != f], kernel = "linear", scale = TRUE)
    fa[f] <- mean(predict(m, X_mod[f_rfe == f, cur, drop = FALSE]) == y_fac[f_rfe == f])
  }
  acc_rfe[i] <- mean(fa)
}
rfe_df <- data.frame(Features = seq_along(svm_genes), Accuracy = acc_rfe)
best_k <- which.max(acc_rfe)

p_m3 <- ggplot(rfe_df, aes(x = Features, y = Accuracy)) +
  geom_line(color = "#3C5488FF", linewidth = 1.2) +
  geom_point(color = "#E64B35FF", size = 2.4) +
  geom_vline(xintercept = best_k, linetype = "dashed", color = "gray45", linewidth = 0.8) +
  annotate("text", x = best_k - 1, y = min(acc_rfe) + 0.006,
           label = sprintf("Peak: %.1f%% (%d genes)", 100 * max(acc_rfe), best_k),
           hjust = 1, size = 3.9, fontface = "bold", color = "gray25") +
  labs(title = "C. SVM-RFE Feature Subset Curve",
       x = "Number of Top Ranked Features", y = "5-Fold CV Accuracy") +
  theme_pub(13)

circle_pts <- function(cx, cy, r = 1.05, n = 200) {
  th <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + r * cos(th), y = cy + r * sin(th))
}
circ_df <- rbind(
  transform(circle_pts(-0.58,  0.35), Set = "LASSO"),
  transform(circle_pts( 0.58,  0.35), Set = "Random Forest"),
  transform(circle_pts( 0.00, -0.52), Set = "SVM-RFE")
)

n_all <- length(Reduce(intersect, list(lasso_genes, rf_genes, svm_genes)))
n_lr  <- length(setdiff(intersect(lasso_genes, rf_genes), svm_genes))
n_ls  <- length(setdiff(intersect(lasso_genes, svm_genes), rf_genes))
n_rs  <- length(setdiff(intersect(rf_genes, svm_genes), lasso_genes))
n_l   <- length(setdiff(lasso_genes, union(rf_genes, svm_genes)))
n_r   <- length(setdiff(rf_genes, union(lasso_genes, svm_genes)))
n_s   <- length(setdiff(svm_genes, union(lasso_genes, rf_genes)))

lbl_venn <- data.frame(
  x   = c(-1.00, 1.00,  0.00, 0.00, -0.50,  0.50, 0.00),
  y   = c( 0.50, 0.50, -0.98, 0.86, -0.22, -0.22, 0.08),
  txt = as.character(c(n_l, n_r, n_s, n_lr, n_ls, n_rs, n_all)),
  sz  = c(5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 6.2)
)

p_m4 <- ggplot() +
  geom_polygon(data = circ_df, aes(x = x, y = y, fill = Set, color = Set), alpha = 0.35, linewidth = 1.1) +
  scale_fill_manual(values  = c("LASSO" = "#E64B35FF", "Random Forest" = "#00A087FF", "SVM-RFE" = "#3C5488FF")) +
  scale_color_manual(values = c("LASSO" = "#B03020FF", "Random Forest" = "#007360FF", "SVM-RFE" = "#26375AFF")) +
  geom_text(data = lbl_venn, aes(x = x, y = y, label = txt, size = sz), fontface = "bold", show.legend = FALSE) +
  scale_size_identity() +
  annotate("text", x = -1.18, y =  1.58, label = sprintf("LASSO (n=%d)", length(lasso_genes)), fontface = "bold", size = 4.2, color = "#B03020FF") +
  annotate("text", x =  1.18, y =  1.58, label = sprintf("Random Forest (n=%d)", length(rf_genes)), fontface = "bold", size = 4.2, color = "#007360FF") +
  annotate("text", x =  0.00, y = -1.78, label = sprintf("SVM-RFE (n=%d)", length(svm_genes)), fontface = "bold", size = 4.2, color = "#26375AFF") +
  annotate("text", x =  0.00, y = -2.26,
           label = paste0("Locked 5-Gene Consensus:\n", paste(hub_genes, collapse = ", ")),
           fontface = "bold.italic", size = 4.0, lineheight = 1.15, color = "black") +
  coord_equal(xlim = c(-2.15, 2.15), ylim = c(-2.58, 1.78), clip = "off") +
  labs(title = "D. Tri-ML Consensus Intersection") +
  theme_void(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.margin = ggplot2::margin(10, 14, 12, 14),
        legend.position = "none")

p_rev4 <- (p_m1 | p_m2) / (p_m3 | p_m4)
ggsave("revision_outputs/figures/Figure_Rev4_TriML_Feature_Selection.png",
       plot = p_rev4, width = 13.8, height = 10.6, dpi = 600)

# ==============================================================================
# 5. NESTED 10-FOLD CV ROC COMPARISON (LOCKED 5-GENE SIGNATURE)
# ==============================================================================
lbl_glm_app <- sprintf("GLM Apparent (AUC = %.4f [%.4f-%.4f])",
                       as.numeric(pROC::auc(res_glm_app$roc)), res_glm_app$roc$ci[1], res_glm_app$roc$ci[3])
lbl_svm_oof <- sprintf("SVM-Radial OOF (AUC = %.4f [%.4f-%.4f])",
                       as.numeric(pROC::auc(res_svm_oof$roc)), res_svm_oof$roc$ci[1], res_svm_oof$roc$ci[3])
lbl_glm_oof <- sprintf("GLM OOF (AUC = %.4f [%.4f-%.4f])",
                       as.numeric(pROC::auc(res_glm_oof$roc)), res_glm_oof$roc$ci[1], res_glm_oof$roc$ci[3])
lbl_rf_oof  <- sprintf("Random Forest OOF (AUC = %.4f [%.4f-%.4f])",
                       as.numeric(pROC::auc(res_rf_oof$roc)),  res_rf_oof$roc$ci[1],  res_rf_oof$roc$ci[3])
lbl_ann_oof <- sprintf("ANN OOF (AUC = %.4f [%.4f-%.4f])",
                       as.numeric(pROC::auc(res_ann_oof$roc)), res_ann_oof$roc$ci[1], res_ann_oof$roc$ci[3])

df_roc_2a <- rbind(
  roc_to_df(res_glm_app$roc, lbl_glm_app),
  roc_to_df(res_svm_oof$roc, lbl_svm_oof),
  roc_to_df(res_glm_oof$roc, lbl_glm_oof),
  roc_to_df(res_rf_oof$roc,  lbl_rf_oof),
  roc_to_df(res_ann_oof$roc, lbl_ann_oof)
)
df_roc_2a$Model <- factor(df_roc_2a$Model,
                          levels = c(lbl_glm_app, lbl_svm_oof, lbl_glm_oof, lbl_rf_oof, lbl_ann_oof))

p_2a <- ggplot(df_roc_2a, aes(x = FPR, y = TPR, color = Model, linetype = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray55", linewidth = 0.8) +
  geom_line(linewidth = 1.25) +
  scale_color_manual(values = setNames(c("gray35", "#3C5488FF", "#E64B35FF", "#00A087FF", "#7E6148FF"),
                                       levels(df_roc_2a$Model))) +
  scale_linetype_manual(values = setNames(c("dashed", "solid", "solid", "solid", "solid"),
                                          levels(df_roc_2a$Model))) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Nested 10-Fold CV Benchmarks (Locked 5-Gene Signature)",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)",
       color = NULL, linetype = NULL) +
  theme_pub(14) +
  theme(legend.position = c(0.61, 0.20), legend.key.width = unit(1.4, "cm"))

ggsave("revision_outputs/figures/Figure_Step2A_NestedCV_ROC.png",
       plot = p_2a, width = 7.5, height = 7.0, dpi = 600)

# ==============================================================================
# 6. AGGREGATED 10-FOLD SHUFFLED ANN TRAINING CURVES
# ==============================================================================
fold_dfs <- lapply(seq_along(ann_fold_histories), function(f) {
  h         <- ann_fold_histories[[f]]
  acc_k     <- intersect(c("accuracy", "acc"), names(h))[1]
  val_acc_k <- intersect(c("val_accuracy", "val_acc"), names(h))[1]
  n_ep      <- length(h$loss)
  data.frame(
    Fold       = f,
    Epoch      = seq_len(n_ep),
    Train_Loss = as.numeric(h$loss),
    Val_Loss   = as.numeric(h$val_loss),
    Train_Acc  = as.numeric(h[[acc_k]]),
    Val_Acc    = as.numeric(h[[val_acc_k]])
  )
})
ann_all     <- do.call(rbind, fold_dfs)
stop_epochs <- sapply(fold_dfs, nrow)
ep_counts   <- table(ann_all$Epoch)
ann_sub     <- ann_all[ann_all$Epoch %in% as.integer(names(ep_counts)[ep_counts >= 5]), ]

agg_curve <- function(df, col_tr, col_va, metric_label) {
  eps <- sort(unique(df$Epoch))
  mk <- function(col, split_lbl) {
    m <- sapply(eps, function(e) mean(df[[col]][df$Epoch == e], na.rm = TRUE))
    s <- sapply(eps, function(e) sd(df[[col]][df$Epoch == e], na.rm = TRUE))
    data.frame(Epoch = eps, Mean = m, Lower = m - s, Upper = m + s, Split = split_lbl, Metric = metric_label)
  }
  rbind(mk(col_tr, "Training"), mk(col_va, "Validation"))
}

loss_sum <- agg_curve(ann_sub, "Train_Loss", "Val_Loss", "Loss")
acc_sum  <- agg_curve(ann_sub, "Train_Acc",  "Val_Acc",  "Accuracy")

p_loss <- ggplot(loss_sum, aes(x = Epoch, y = Mean, color = Split, fill = Split, linetype = Split)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.18, color = NA) +
  geom_line(linewidth = 1.25) +
  scale_color_manual(values = c("Training" = "#3C5488FF", "Validation" = "#E64B35FF")) +
  scale_fill_manual(values  = c("Training" = "#3C5488FF", "Validation" = "#E64B35FF")) +
  scale_linetype_manual(values = c("Training" = "solid", "Validation" = "dashed")) +
  labs(title = "A. ANN Training vs. Validation Loss",
       subtitle = sprintf("Mean ± SD across 10 folds (median stop: %d epochs)", round(median(stop_epochs))),
       x = "Epoch", y = "Binary Cross-Entropy Loss", color = NULL, fill = NULL, linetype = NULL) +
  theme_pub(14) +
  theme(legend.position = c(0.76, 0.84))

p_acc <- ggplot(acc_sum, aes(x = Epoch, y = Mean, color = Split, fill = Split, linetype = Split)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.18, color = NA) +
  geom_line(linewidth = 1.25) +
  scale_color_manual(values = c("Training" = "#3C5488FF", "Validation" = "#E64B35FF")) +
  scale_fill_manual(values  = c("Training" = "#3C5488FF", "Validation" = "#E64B35FF")) +
  scale_linetype_manual(values = c("Training" = "solid", "Validation" = "dashed")) +
  labs(title = "B. ANN Training vs. Validation Accuracy",
       subtitle = "Shuffled in-fold splits with dropout & early stopping",
       x = "Epoch", y = "Classification Accuracy", color = NULL, fill = NULL, linetype = NULL) +
  theme_pub(14) +
  theme(legend.position = c(0.76, 0.20))

ggsave("revision_outputs/figures/Figure_Step2B_ANN_Training_Curves.png",
       plot = p_loss | p_acc, width = 12.0, height = 5.2, dpi = 600)

# ==============================================================================
# 7. END-TO-END 5-FOLD NESTED CV ROC & FEATURE SELECTION STABILITY
# ==============================================================================
roc_e2e_glm <- roc(y_bin, e2e_oof_glm, levels = c(0, 1), direction = "<", ci = TRUE, quiet = TRUE)
roc_e2e_rf  <- roc(y_bin, e2e_oof_rf,  levels = c(0, 1), direction = "<", ci = TRUE, quiet = TRUE)
roc_e2e_svm <- roc(y_bin, e2e_oof_svm, levels = c(0, 1), direction = "<", ci = TRUE, quiet = TRUE)

l_eg <- sprintf("End-to-End GLM (AUC = %.4f [%.4f-%.4f])", as.numeric(pROC::auc(roc_e2e_glm)), roc_e2e_glm$ci[1], roc_e2e_glm$ci[3])
l_er <- sprintf("End-to-End RF (AUC = %.4f [%.4f-%.4f])",  as.numeric(pROC::auc(roc_e2e_rf)),  roc_e2e_rf$ci[1],  roc_e2e_rf$ci[3])
l_es <- sprintf("End-to-End SVM (AUC = %.4f [%.4f-%.4f])", as.numeric(pROC::auc(roc_e2e_svm)), roc_e2e_svm$ci[1], roc_e2e_svm$ci[3])

df_roc_2c <- rbind(
  roc_to_df(roc_e2e_glm, l_eg),
  roc_to_df(roc_e2e_rf,  l_er),
  roc_to_df(roc_e2e_svm, l_es)
)
df_roc_2c$Model <- factor(df_roc_2c$Model, levels = c(l_eg, l_er, l_es))

p_2c_left <- ggplot(df_roc_2c, aes(x = FPR, y = TPR, color = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray55", linewidth = 0.8) +
  geom_line(linewidth = 1.25) +
  scale_color_manual(values = setNames(c("#E64B35FF", "#00A087FF", "#3C5488FF"), levels(df_roc_2c$Model))) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "A. End-to-End Leakage-Free 5-Fold CV",
       subtitle = "WGCNA + Tri-ML repeated strictly inside training folds",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)",
       color = NULL) +
  theme_pub(13) +
  theme(legend.position = c(0.57, 0.16))

stab_long <- rbind(
  data.frame(Gene = stability_df$Gene, Method = "LASSO",             Count = stability_df$LASSO_Freq),
  data.frame(Gene = stability_df$Gene, Method = "Random Forest",     Count = stability_df$RF_Freq),
  data.frame(Gene = stability_df$Gene, Method = "SVM-RFE",           Count = stability_df$SVM_RFE_Freq),
  data.frame(Gene = stability_df$Gene, Method = "Tri-ML\nConsensus", Count = stability_df$Consensus_Freq)
)
stab_long$Method <- factor(stab_long$Method,
                           levels = c("LASSO", "Random Forest", "SVM-RFE", "Tri-ML\nConsensus"))
stab_long$Gene   <- factor(stab_long$Gene, levels = rev(c("CD38", "IGHM", "FAM30A", "P2RX5", "IGHG3")))
stab_long$Label  <- paste0(stab_long$Count, "/5")

p_2c_right <- ggplot(stab_long, aes(x = Method, y = Gene, fill = Count)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = Label), size = 5.2, fontface = "bold",
            color = ifelse(stab_long$Count >= 4, "white", "black")) +
  scale_fill_gradient(low = "#EDF2F8", high = "#3C5488FF", limits = c(0, 5), breaks = 0:5) +
  labs(title = "B. Feature Selection Stability (5 Folds)",
       subtitle = "Re-selection frequency in independent training folds",
       x = NULL, y = NULL, fill = "Folds\n(of 5)") +
  theme_pub(13) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(face = "bold.italic", size = 13),
        axis.text.x = element_text(face = "bold", size = 12))

ggsave("revision_outputs/figures/Figure_Step2C_EndToEnd_NestedCV_Stability.png",
       plot = p_2c_left | p_2c_right, width = 13.0, height = 5.8, dpi = 600)

# ==============================================================================
# 8. PRE-COMBAT LEAVE-ONE-COHORT-OUT (LOCO) VALIDATION
# ==============================================================================
loco_long <- rbind(
  data.frame(Cohort = loco_per_cohort$HeldOut_Cohort, Model = "GLM",        AUC = loco_per_cohort$GLM_AUC),
  data.frame(Cohort = loco_per_cohort$HeldOut_Cohort, Model = "RF",         AUC = loco_per_cohort$RF_AUC),
  data.frame(Cohort = loco_per_cohort$HeldOut_Cohort, Model = "SVM-Radial", AUC = loco_per_cohort$SVM_AUC),
  data.frame(Cohort = loco_per_cohort$HeldOut_Cohort, Model = "ANN",        AUC = loco_per_cohort$ANN_AUC)
)
loco_long$Model      <- factor(loco_long$Model, levels = c("GLM", "RF", "SVM-Radial", "ANN"))
loco_long$Cohort_Lbl <- paste0(
  loco_long$Cohort, "\n(n=",
  loco_per_cohort$N_Total[match(loco_long$Cohort, loco_per_cohort$HeldOut_Cohort)], ")"
)

p_2d_left <- ggplot(loco_long, aes(x = Cohort_Lbl, y = AUC, fill = Model)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.70, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", AUC)),
            position = position_dodge(width = 0.78), vjust = -0.4, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("GLM" = "#E64B35FF", "RF" = "#00A087FF",
                               "SVM-Radial" = "#3C5488FF", "ANN" = "#7E6148FF")) +
  coord_cartesian(ylim = c(0.70, 0.96)) +
  labs(title = "A. Per-Cohort LOCO AUC (Pre-ComBat fRMA)",
       subtitle = "Each cohort held out & Z-scaled independently",
       x = "Held-Out GEO Cohort", y = "ROC AUC", fill = "Classifier") +
  theme_pub(13) +
  theme(legend.position = "top")

roc_l_glm <- roc(y_bin, loco_oof_glm, levels = c(0, 1), direction = "<", ci = TRUE, quiet = TRUE)
roc_l_rf  <- roc(y_bin, loco_oof_rf,  levels = c(0, 1), direction = "<", ci = TRUE, quiet = TRUE)
roc_l_svm <- roc(y_bin, loco_oof_svm, levels = c(0, 1), direction = "<", ci = TRUE, quiet = TRUE)
roc_l_ann <- roc(y_bin, loco_oof_ann, levels = c(0, 1), direction = "<", ci = TRUE, quiet = TRUE)

ll_ann <- sprintf("ANN (AUC = %.4f [%.4f-%.4f])",        as.numeric(pROC::auc(roc_l_ann)), roc_l_ann$ci[1], roc_l_ann$ci[3])
ll_svm <- sprintf("SVM-Radial (AUC = %.4f [%.4f-%.4f])", as.numeric(pROC::auc(roc_l_svm)), roc_l_svm$ci[1], roc_l_svm$ci[3])
ll_glm <- sprintf("GLM (AUC = %.4f [%.4f-%.4f])",        as.numeric(pROC::auc(roc_l_glm)), roc_l_glm$ci[1], roc_l_glm$ci[3])
ll_rf  <- sprintf("RF (AUC = %.4f [%.4f-%.4f])",         as.numeric(pROC::auc(roc_l_rf)),  roc_l_rf$ci[1],  roc_l_rf$ci[3])

df_roc_2d <- rbind(
  roc_to_df(roc_l_ann, ll_ann),
  roc_to_df(roc_l_svm, ll_svm),
  roc_to_df(roc_l_glm, ll_glm),
  roc_to_df(roc_l_rf,  ll_rf)
)
df_roc_2d$Model <- factor(df_roc_2d$Model, levels = c(ll_ann, ll_svm, ll_glm, ll_rf))

p_2d_right <- ggplot(df_roc_2d, aes(x = FPR, y = TPR, color = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray55", linewidth = 0.8) +
  geom_line(linewidth = 1.25) +
  scale_color_manual(values = setNames(c("#7E6148FF", "#3C5488FF", "#E64B35FF", "#00A087FF"), levels(df_roc_2d$Model))) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "B. Pooled LOCO ROC Curves (n = 536)",
       subtitle = "Cross-cohort generalization without ComBat",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)",
       color = NULL) +
  theme_pub(13) +
  theme(legend.position = c(0.60, 0.18))

ggsave("revision_outputs/figures/Figure_Step2D_LOCO_Validation.png",
       plot = p_2d_left | p_2d_right, width = 13.2, height = 5.8, dpi = 600)

message("Generated all publication figures in revision_outputs/figures/")
