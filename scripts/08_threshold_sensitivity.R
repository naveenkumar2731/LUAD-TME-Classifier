# ==============================================================================
# 08_threshold_sensitivity.R
# Sensitivity analysis of immune phenotype definitions and decision cutoffs:
#   Part A: Clustering algorithm perturbation (k = 2) -> ARI, Kappa, 10-fold OOF AUC
#   Part B: Direct quantile cutoffs (CD8+ only vs. multi-lineage effector score)
#   Part C: Three-tier Ward.D2 stratification (k = 3) monotonicity check
#   Part D: Classifier probability decision-threshold sweep (GEO OOF & TCGA)
#   Part E: High-resolution (600 DPI) summary figure
# ==============================================================================

suppressPackageStartupMessages({
  library(cluster)
  library(pROC)
  library(e1071)
  library(ggplot2)
  library(patchwork)
})

dir.create("revision_outputs/figures", showWarnings = FALSE, recursive = TRUE)

load("Phase3_Immune_Phenotypes.RData")
load("Phase5_Final_Hub_Genes.RData")
load("revision_outputs/Step2_LeakageFree_Results.RData")
load("revision_outputs/Step3_TCGA_External_Results.RData")

stopifnot(all(rownames(immune_fractions) == pheno_immune$Sample))
stopifnot(all(colnames(harmonized_matrix) == pheno_immune$Sample))

X_5    <- t(harmonized_matrix[hub_genes, , drop = FALSE])
Z_lm22 <- scale(immune_fractions)
Z_lm22[!is.finite(Z_lm22)] <- 0
y_base <- as.integer(pheno_immune$Immune_Phenotype == "Immune-Hot")

calc_ari <- function(a, b) {
  tab   <- table(a, b)
  n     <- sum(tab)
  c2    <- function(x) sum(x * (x - 1) / 2)
  idx   <- c2(tab)
  ri    <- c2(rowSums(tab))
  ci    <- c2(colSums(tab))
  tot   <- c2(n)
  exp_i <- (ri * ci) / tot
  max_i <- 0.5 * (ri + ci)
  if (max_i == exp_i) return(0)
  (idx - exp_i) / (max_i - exp_i)
}

calc_kappa <- function(a, b) {
  tab <- table(factor(a, levels = c(0, 1)), factor(b, levels = c(0, 1)))
  n   <- sum(tab)
  po  <- sum(diag(tab)) / n
  pe  <- sum(rowSums(tab) * colSums(tab)) / (n^2)
  if (pe == 1) return(1)
  (po - pe) / (1 - pe)
}

label_hot_by_cd8 <- function(cl_vec, cd8_z) {
  means  <- tapply(cd8_z, cl_vec, mean)
  hot_cl <- names(which.max(means))
  as.integer(as.character(cl_vec) == hot_cl)
}

eval_10fold_oof <- function(X_mat, y_bin, seed = 42L) {
  set.seed(seed)
  k     <- 10L
  folds <- integer(length(y_bin))
  for (cls in unique(y_bin)) {
    idx <- which(y_bin == cls)
    folds[idx] <- sample(rep(1:k, length.out = length(idx)))
  }
  p_glm <- numeric(length(y_bin))
  p_svm <- numeric(length(y_bin))
  y_fac <- factor(ifelse(y_bin == 1, "Hot", "Cold"), levels = c("Cold", "Hot"))

  for (f in seq_len(k)) {
    tr <- which(folds != f)
    te <- which(folds == f)
    mu <- colMeans(X_mat[tr, , drop = FALSE])
    sd_tr <- apply(X_mat[tr, , drop = FALSE], 2, sd)
    sd_tr[sd_tr == 0] <- 1
    Xtr <- scale(X_mat[tr, , drop = FALSE], center = mu, scale = sd_tr)
    Xte <- scale(X_mat[te, , drop = FALSE], center = mu, scale = sd_tr)

    m_g <- glm(y ~ ., data = data.frame(Xtr, y = y_bin[tr], check.names = TRUE), family = binomial)
    p_glm[te] <- predict(m_g, newdata = data.frame(Xte, check.names = TRUE), type = "response")

    m_s <- svm(x = Xtr, y = y_fac[tr], kernel = "radial", scale = FALSE, probability = TRUE)
    p_svm[te] <- attr(predict(m_s, Xte, probability = TRUE), "probabilities")[, "Hot"]
  }
  r_g <- roc(y_bin, p_glm, levels = c(0, 1), direction = "<", ci = TRUE, quiet = TRUE)
  r_s <- roc(y_bin, p_svm, levels = c(0, 1), direction = "<", ci = TRUE, quiet = TRUE)
  list(
    glm_auc = as.numeric(pROC::auc(r_g)), glm_lo = as.numeric(r_g$ci[1]), glm_hi = as.numeric(r_g$ci[3]),
    svm_auc = as.numeric(pROC::auc(r_s)), svm_lo = as.numeric(r_s$ci[1]), svm_hi = as.numeric(r_s$ci[3]),
    oof_glm = p_glm
  )
}

# ==============================================================================
# PART A: CLUSTERING ALGORITHM PERTURBATION (k = 2)
# ==============================================================================
cd8_z <- Z_lm22[, "T cells CD8"]

y_wardD2 <- y_base
y_wardD  <- label_hot_by_cd8(cutree(hclust(dist(Z_lm22, "euclidean"), method = "ward.D"), k = 2), cd8_z)
y_manh   <- label_hot_by_cd8(cutree(hclust(dist(Z_lm22, "manhattan"), method = "ward.D2"), k = 2), cd8_z)
cor_dist <- as.dist(1 - cor(t(Z_lm22)))
y_cord   <- label_hot_by_cd8(cutree(hclust(cor_dist, method = "ward.D2"), k = 2), cd8_z)

set.seed(42)
y_pam <- label_hot_by_cd8(pam(Z_lm22, k = 2)$clustering, cd8_z)

set.seed(42)
y_kmeans <- label_hot_by_cd8(kmeans(Z_lm22, centers = 2, nstart = 50)$cluster, cd8_z)

cluster_defs <- list(
  "Ward.D2 (Euclidean, Baseline)" = y_wardD2,
  "Ward.D (Euclidean)"            = y_wardD,
  "Ward.D2 (Manhattan)"           = y_manh,
  "Ward.D2 (Correlation 1-r)"     = y_cord,
  "PAM (k-Medoids)"               = y_pam,
  "K-Means (nstart = 50)"         = y_kmeans
)

partA_df <- do.call(rbind, lapply(names(cluster_defs), function(nm) {
  yt <- cluster_defs[[nm]]
  ev <- eval_10fold_oof(X_5, yt)
  data.frame(
    Definition       = nm,
    Category         = "Clustering (k=2)",
    N_Hot            = sum(yt == 1),
    N_Cold           = sum(yt == 0),
    Pct_Hot          = round(100 * mean(yt == 1), 1),
    ARI_vs_Base      = round(calc_ari(y_base, yt), 4),
    Kappa_vs_Base    = round(calc_kappa(y_base, yt), 4),
    GLM_OOF_AUC      = round(ev$glm_auc, 4),
    GLM_CI_Lo        = round(ev$glm_lo, 4),
    GLM_CI_Hi        = round(ev$glm_hi, 4),
    SVM_OOF_AUC      = round(ev$svm_auc, 4),
    stringsAsFactors = FALSE
  )
}))

# ==============================================================================
# PART B: DIRECT QUANTILE CUTOFFS (CD8+ ONLY vs. MULTI-LINEAGE EFFECTOR SCORE)
# ==============================================================================
cd8_raw   <- immune_fractions[, "T cells CD8"]
eff_cols  <- c("T cells CD8", "T cells CD4 memory activated", "Macrophages M1",
               "T cells follicular helper", "Plasma cells")
eff_score <- rowMeans(Z_lm22[, eff_cols, drop = FALSE])

cut_probs <- c(0.75, 0.67, 0.625, 0.50)
cut_names <- c("Top 25%", "Top 33%", "Top 37.5%", "Top 50% (Median)")

eval_quantiles <- function(score_vec, prefix, cat_lbl) {
  do.call(rbind, lapply(seq_along(cut_probs), function(i) {
    thr <- quantile(score_vec, probs = cut_probs[i], names = FALSE)
    yt  <- as.integer(score_vec >= thr)
    ev  <- eval_10fold_oof(X_5, yt)
    data.frame(
      Definition       = paste0(prefix, " (", cut_names[i], ")"),
      Category         = cat_lbl,
      N_Hot            = sum(yt == 1),
      N_Cold           = sum(yt == 0),
      Pct_Hot          = round(100 * mean(yt == 1), 1),
      ARI_vs_Base      = round(calc_ari(y_base, yt), 4),
      Kappa_vs_Base    = round(calc_kappa(y_base, yt), 4),
      GLM_OOF_AUC      = round(ev$glm_auc, 4),
      GLM_CI_Lo        = round(ev$glm_lo, 4),
      GLM_CI_Hi        = round(ev$glm_hi, 4),
      SVM_OOF_AUC      = round(ev$svm_auc, 4),
      stringsAsFactors = FALSE
    )
  }))
}

partB_df <- rbind(
  eval_quantiles(cd8_raw,   "CD8+ T-cell",            "CD8+ Quantile"),
  eval_quantiles(eff_score, "Multi-Lineage Effector", "Multi-Lineage Quantile")
)

sensitivity_master <- rbind(partA_df, partB_df)
write.csv(sensitivity_master, "revision_outputs/Table_Step4A_Phenotype_Sensitivity.csv", row.names = FALSE)

# ==============================================================================
# PART C: THREE-TIER STRATIFICATION (k = 3)
# ==============================================================================
hc_base <- hclust(dist(Z_lm22, method = "euclidean"), method = "ward.D2")
cl3     <- cutree(hc_base, k = 3)
eff_3   <- tapply(eff_score, cl3, mean)
ord3    <- names(sort(eff_3, decreasing = FALSE))

tier3 <- factor(
  ifelse(cl3 == ord3[1], "Immune-Cold",
         ifelse(cl3 == ord3[2], "Immune-Intermediate", "Immune-Hot")),
  levels = c("Immune-Cold", "Immune-Intermediate", "Immune-Hot")
)

oof_clip      <- pmin(pmax(oof_glm, 1e-4), 1 - 1e-4)
hot_score_oof <- log(oof_clip / (1 - oof_clip))

tier_df <- data.frame(
  Sample    = pheno_immune$Sample,
  Tier3     = tier3,
  Tier_Num  = as.integer(tier3),
  Hot_Score = hot_score_oof,
  OOF_Prob  = oof_glm,
  CD8_Frac  = cd8_raw,
  Eff_Score = eff_score
)

tier_summary <- data.frame(
  Tier             = levels(tier3),
  N                = as.integer(table(tier3)),
  Median_Hot_Score = round(tapply(tier_df$Hot_Score, tier_df$Tier3, median), 3),
  IQR_Hot_Score    = round(tapply(tier_df$Hot_Score, tier_df$Tier3, IQR), 3),
  Mean_CD8_Frac    = round(tapply(tier_df$CD8_Frac,  tier_df$Tier3, mean), 4),
  Mean_Eff_Z       = round(tapply(tier_df$Eff_Score, tier_df$Tier3, mean), 3)
)

kw_3     <- kruskal.test(Hot_Score ~ Tier3, data = tier_df)
trend_3  <- suppressWarnings(cor.test(tier_df$Tier_Num, tier_df$Hot_Score, method = "spearman"))
pw_3     <- pairwise.wilcox.test(tier_df$Hot_Score, tier_df$Tier3, p.adjust.method = "BH")
ext_idx  <- tier_df$Tier3 %in% c("Immune-Cold", "Immune-Hot")
roc_ext3 <- roc(as.integer(tier_df$Tier3[ext_idx] == "Immune-Hot"),
                tier_df$OOF_Prob[ext_idx], levels = c(0, 1), direction = "<", ci = TRUE, quiet = TRUE)

write.csv(tier_summary, "revision_outputs/Table_Step4B_3Tier_Summary.csv", row.names = FALSE)

# ==============================================================================
# PART D: CLASSIFIER DECISION-THRESHOLD SWEEP
# ==============================================================================
sweep_thresholds <- function(y_true, probs, thrs = seq(0.20, 0.80, by = 0.05)) {
  do.call(rbind, lapply(thrs, function(t) {
    pr   <- as.integer(probs >= t)
    tp   <- sum(pr == 1 & y_true == 1)
    tn   <- sum(pr == 0 & y_true == 0)
    fp   <- sum(pr == 1 & y_true == 0)
    fn   <- sum(pr == 0 & y_true == 1)
    sens <- tp / (tp + fn)
    spec <- tn / (tn + fp)
    prec <- ifelse((tp + fp) > 0, tp / (tp + fp), NA)
    f1   <- ifelse(!is.na(prec) && (prec + sens) > 0, 2 * prec * sens / (prec + sens), NA)
    data.frame(
      Threshold         = round(t, 2),
      Accuracy          = round((tp + tn) / length(y_true), 4),
      Balanced_Accuracy = round((sens + spec) / 2, 4),
      Sensitivity       = round(sens, 4),
      Specificity       = round(spec, 4),
      F1_Score          = round(f1, 4),
      Youden_J          = round(sens + spec - 1, 4)
    )
  }))
}

geo_sweep   <- sweep_thresholds(y_base, oof_glm)
tcga_sweep  <- sweep_thresholds(y_tcga_bin, tcga_prob_glm)
best_geo_t  <- geo_sweep[which.max(geo_sweep$Youden_J), ]
best_tcga_t <- tcga_sweep[which.max(tcga_sweep$Youden_J), ]

write.csv(geo_sweep,  "revision_outputs/Table_Step4C_GEO_ProbThreshold_Sweep.csv",  row.names = FALSE)
write.csv(tcga_sweep, "revision_outputs/Table_Step4C_TCGA_ProbThreshold_Sweep.csv", row.names = FALSE)

save(
  partA_df, partB_df, sensitivity_master,
  tier_df, tier_summary, kw_3, trend_3, pw_3, roc_ext3,
  geo_sweep, tcga_sweep, best_geo_t, best_tcga_t,
  file = "revision_outputs/Step4_Sensitivity_Results.RData"
)

# ==============================================================================
# PART E: HIGH-RESOLUTION (600 DPI) THRESHOLD SENSITIVITY FIGURE
# ==============================================================================
theme_pub <- function(base_size = 13) {
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

plot_a_df            <- sensitivity_master
plot_a_df$Definition <- factor(plot_a_df$Definition, levels = rev(plot_a_df$Definition))
plot_a_df$Category   <- factor(plot_a_df$Category,
                               levels = c("Clustering (k=2)", "Multi-Lineage Quantile", "CD8+ Quantile"))

p_4a <- ggplot(plot_a_df, aes(x = GLM_OOF_AUC, y = Definition, color = Category)) +
  geom_vline(xintercept = partA_df$GLM_OOF_AUC[1], linetype = "dashed", color = "gray45", linewidth = 0.7) +
  geom_errorbar(aes(xmin = GLM_CI_Lo, xmax = GLM_CI_Hi), orientation = "y", width = 0.28, linewidth = 0.95) +
  geom_point(size = 3.3) +
  geom_text(aes(x = GLM_CI_Hi + 0.012, label = sprintf("%.3f", GLM_OOF_AUC)),
            hjust = 0, size = 3.6, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c("Clustering (k=2)"       = "#3C5488FF",
                                "Multi-Lineage Quantile" = "#E64B35FF",
                                "CD8+ Quantile"          = "#00A087FF")) +
  coord_cartesian(xlim = c(0.50, 0.98)) +
  labs(
    title = "A. Sensitivity to Immune Phenotype Definition",
    subtitle = "10-Fold OOF GLM AUC (95% CI) of locked 5-gene panel",
    x = "10-Fold OOF ROC AUC", y = NULL, color = NULL
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_pub(13) +
  theme(legend.position = "bottom")

p_4b <- ggplot(tier_df, aes(x = Tier3, y = Hot_Score, fill = Tier3)) +
  geom_violin(alpha = 0.35, trim = FALSE, linewidth = 0.6) +
  geom_boxplot(width = 0.22, outlier.size = 1.2, alpha = 0.85) +
  scale_fill_manual(values = c("Immune-Cold"         = "#3C5488FF",
                               "Immune-Intermediate" = "#F39B7FFF",
                               "Immune-Hot"          = "#E64B35FF")) +
  scale_x_discrete(labels = c(
    sprintf("Cold\n(n=%d)", tier_summary$N[1]),
    sprintf("Intermediate\n(n=%d)", tier_summary$N[2]),
    sprintf("Hot\n(n=%d)", tier_summary$N[3])
  )) +
  labs(
    title = "B. 3-Tier Stratification (k = 3)",
    subtitle = sprintf("KW p = %.2e | ρ = %.3f | Extreme AUC = %.3f",
                       kw_3$p.value, trend_3$estimate, as.numeric(pROC::auc(roc_ext3))),
    x = "Ward.D2 TME Tier (k = 3)",
    y = "5-Gene OOF Hot Score (Log-Odds)"
  ) +
  theme_pub(13) +
  theme(legend.position = "none")

sweep_long <- rbind(
  data.frame(Threshold = geo_sweep$Threshold, Metric = "Balanced Accuracy", Value = geo_sweep$Balanced_Accuracy),
  data.frame(Threshold = geo_sweep$Threshold, Metric = "Sensitivity",       Value = geo_sweep$Sensitivity),
  data.frame(Threshold = geo_sweep$Threshold, Metric = "Specificity",       Value = geo_sweep$Specificity),
  data.frame(Threshold = geo_sweep$Threshold, Metric = "F1 Score",          Value = geo_sweep$F1_Score)
)
sweep_long$Metric <- factor(sweep_long$Metric,
                            levels = c("Balanced Accuracy", "Sensitivity", "Specificity", "F1 Score"))

p_4c <- ggplot(sweep_long, aes(x = Threshold, y = Value, color = Metric, linetype = Metric)) +
  geom_vline(xintercept = 0.50, linetype = "dotted", color = "gray40", linewidth = 0.8) +
  geom_vline(xintercept = best_geo_t$Threshold, linetype = "dashed", color = "#E64B35FF", linewidth = 0.8) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.0) +
  scale_color_manual(values = c("Balanced Accuracy" = "#000000",
                                "Sensitivity"       = "#E64B35FF",
                                "Specificity"       = "#3C5488FF",
                                "F1 Score"          = "#00A087FF")) +
  scale_linetype_manual(values = c("Balanced Accuracy" = "solid",
                                   "Sensitivity"       = "solid",
                                   "Specificity"       = "solid",
                                   "F1 Score"          = "dashed")) +
  annotate("text", x = best_geo_t$Threshold - 0.015, y = 0.43,
           label = sprintf("Youden opt.\n(t = %.2f)", best_geo_t$Threshold),
           hjust = 1, size = 3.6, fontface = "bold", color = "#E64B35FF") +
  coord_cartesian(ylim = c(0.35, 1.0)) +
  labs(
    title = "C. Probability Cutoff Sensitivity (GEO OOF)",
    subtitle = "Stability of operating metrics across decision thresholds",
    x = "Predicted Immune-Hot Probability Cutoff", y = "Metric Value",
    color = NULL, linetype = NULL
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_pub(13) +
  theme(legend.position = "bottom")

p_4 <- (p_4a | (p_4b / p_4c)) + plot_layout(widths = c(1.38, 1))
ggsave("revision_outputs/figures/Figure_Step4_Threshold_Sensitivity.png",
       plot = p_4, width = 14.6, height = 8.8, dpi = 600)
message("Saved revision_outputs/Step4_Sensitivity_Results.RData and Figure_Step4_Threshold_Sensitivity.png")
