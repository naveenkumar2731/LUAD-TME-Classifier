# ==============================================================================
# 05_triml_feature_selection.R
# Consensus Tri-ML feature selection (LASSO, Random Forest, and Linear SVM-RFE)
# on candidate WGCNA module genes, followed by SVM-RFE subset cross-validation
# and multivariable logistic regression on the locked five-gene signature.
# ==============================================================================

suppressPackageStartupMessages({
  library(glmnet)
  library(randomForest)
  library(e1071)
  library(pROC)
})

# 1. Load Phase 3 phenotypes and Phase 4 candidate module genes
load("Phase3_Immune_Phenotypes.RData")
load("Phase4_WGCNA_Modules_Corrected.RData")

X      <- t(harmonized_matrix[target_genes, , drop = FALSE])
y      <- factor(pheno_immune$Immune_Phenotype, levels = c("Immune-Cold", "Immune-Hot"))
y_cold <- ifelse(y == "Immune-Cold", 1, 0)
y_hot  <- ifelse(y == "Immune-Hot",  1, 0)

# 2. Algorithm 1: LASSO penalized logistic regression (10-fold CV at lambda.min)
set.seed(123)
cv_lasso    <- cv.glmnet(x = as.matrix(X), y = y_cold, family = "binomial",
                         alpha = 1, nfolds = 10)
lasso_coefs <- coef(cv_lasso, s = "lambda.min")
lasso_genes <- setdiff(rownames(lasso_coefs)[which(lasso_coefs[, 1] != 0)], "(Intercept)")
message("LASSO selected ", length(lasso_genes), " genes.")

# 3. Algorithm 2: Random Forest feature importance (Top 30 by MeanDecreaseGini)
set.seed(123)
rf_model      <- randomForest(x = X, y = y, ntree = 1000, importance = TRUE)
rf_importance <- importance(rf_model)
rf_ranked     <- sort(rf_importance[, "MeanDecreaseGini"], decreasing = TRUE)
rf_genes      <- names(rf_ranked[1:min(30L, length(rf_ranked))])
message("Random Forest selected top ", length(rf_genes), " genes.")

# 4. Algorithm 3: Linear SVM-RFE (Top 30 ranked features)
svm_rfe_ranking <- function(x, y) {
  feature_list    <- colnames(x)
  ranked_features <- character(0)

  while (length(feature_list) > 0) {
    svm_fit <- svm(x = x[, feature_list, drop = FALSE], y = y,
                   kernel = "linear", scale = TRUE)
    w       <- t(svm_fit$coefs) %*% svm_fit$SV
    weights <- as.numeric(w[1, ])^2

    worst_feature   <- feature_list[which.min(weights)]
    ranked_features <- c(worst_feature, ranked_features)
    feature_list    <- setdiff(feature_list, worst_feature)
  }
  ranked_features
}

ranked_svm_genes <- svm_rfe_ranking(X, y)
svm_genes        <- ranked_svm_genes[1:min(30L, length(ranked_svm_genes))]
message("SVM-RFE selected top ", length(svm_genes), " genes.")

# 5. Strict three-way consensus intersection
hub_genes <- Reduce(intersect, list(lasso_genes, rf_genes, svm_genes))
message("Locked consensus hub genes: ", paste(hub_genes, collapse = ", "))

# 6. Apparent multivariable logistic regression on locked signature
hub_expr           <- as.data.frame(X[, hub_genes, drop = FALSE])
hub_expr$Phenotype <- y_cold

fit_hub   <- glm(Phenotype ~ ., data = hub_expr, family = binomial)
pred_prob <- predict(fit_hub, type = "response")
roc_obj   <- roc(y_cold, pred_prob, ci = TRUE, quiet = TRUE)

message(sprintf("Apparent GLM AUC = %.4f (95%% CI: %.4f-%.4f)",
                as.numeric(auc(roc_obj)), roc_obj$ci[1], roc_obj$ci[3]))

# 7. SVM-RFE 5-fold CV subset accuracy curve and top-30 OOF ROC
set.seed(123)
k_folds <- 5L
folds   <- sample(rep(1:k_folds, length.out = nrow(X)))

rfe_accuracy <- numeric(length(svm_genes))
for (i in seq_along(svm_genes)) {
  cur_feats <- svm_genes[1:i]
  fold_acc  <- numeric(k_folds)
  for (f in seq_len(k_folds)) {
    m <- svm(x = X[folds != f, cur_feats, drop = FALSE], y = y[folds != f],
             kernel = "linear", scale = TRUE)
    fold_acc[f] <- mean(predict(m, X[folds == f, cur_feats, drop = FALSE]) == y[folds == f])
  }
  rfe_accuracy[i] <- mean(fold_acc)
}

svm_rfe_oof <- numeric(nrow(X))
for (f in seq_len(k_folds)) {
  m_cv <- svm(x = X[folds != f, svm_genes, drop = FALSE], y = y[folds != f],
              kernel = "linear", scale = TRUE, probability = TRUE)
  p_cv <- predict(m_cv, X[folds == f, svm_genes, drop = FALSE], probability = TRUE)
  svm_rfe_oof[folds == f] <- attr(p_cv, "probabilities")[, "Immune-Hot"]
}
svm_rfe_roc <- roc(y_hot, svm_rfe_oof, levels = c(0, 1), direction = "<", quiet = TRUE)
message(sprintf("SVM-RFE peak 5-fold CV accuracy = %.2f%% (%d genes) | Top-30 OOF AUC = %.4f",
                100 * max(rfe_accuracy), which.max(rfe_accuracy), as.numeric(auc(svm_rfe_roc))))

# 8. Save Phase 5 outputs
save(hub_genes, lasso_genes, rf_genes, svm_genes, roc_obj, fit_hub,
     rfe_accuracy, svm_rfe_roc,
     file = "Phase5_Final_Hub_Genes.RData")
write.csv(data.frame(Hub_Gene = hub_genes),
          file = "Final_Core_Immune_Evasion_Signature.csv", row.names = FALSE)
message("Saved Phase5_Final_Hub_Genes.RData")