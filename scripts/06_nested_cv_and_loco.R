# ==============================================================================
# 06_nested_cv_and_loco.R
# Leakage-controlled internal and cross-cohort validation:
#   Part A: Strictly nested 10-fold CV on the locked five-gene signature
#   Part B: End-to-end 5-fold nested CV (WGCNA + Tri-ML inside training folds)
#   Part C: Leave-One-Cohort-Out (LOCO) validation on pre-ComBat fRMA data
# ==============================================================================

suppressPackageStartupMessages({
  library(reticulate)
  library(keras)
})

if (Sys.getenv("RETICULATE_PYTHON") == "") {
  if (exists("py_require_legacy_keras", mode = "function")) {
    py_require_legacy_keras()
  } else {
    Sys.setenv(TF_USE_LEGACY_KERAS = "1")
    reticulate::py_require("tf-keras")
  }
  reticulate::py_require(c("tensorflow", "numpy"))
}

invisible(callback_early_stopping(monitor = "val_loss", patience = as.integer(15)))

suppressPackageStartupMessages({
  library(pROC)
  library(glmnet)
  library(randomForest)
  library(e1071)
  library(WGCNA)
})

allowWGCNAThreads()
dir.create("revision_outputs", showWarnings = FALSE, recursive = TRUE)

# 1. Load verified Phase 3 and Phase 5 outputs
load("Phase3_Immune_Phenotypes.RData")
load("Phase5_Final_Hub_Genes.RData")

stopifnot(all(colnames(harmonized_matrix) == pheno_immune$Sample))

X_all   <- t(harmonized_matrix)
X_5     <- X_all[, hub_genes, drop = FALSE]
y_fac   <- factor(pheno_immune$Immune_Phenotype, levels = c("Immune-Cold", "Immune-Hot"))
y_bin   <- as.integer(y_fac == "Immune-Hot")
batches <- pheno_immune$Batch

make_stratified_folds <- function(y, k = 10L, seed = 42L) {
  set.seed(seed)
  folds <- integer(length(y))
  for (cls in unique(y)) {
    idx <- which(y == cls)
    folds[idx] <- sample(rep(1:k, length.out = length(idx)))
  }
  folds
}

compute_metrics <- function(y_true, probs, model_name, eval_type) {
  roc_o    <- roc(y_true, probs, levels = c(0, 1), direction = "<", ci = TRUE, quiet = TRUE)
  pred_cls <- as.integer(probs >= 0.5)
  tp <- sum(pred_cls == 1 & y_true == 1)
  tn <- sum(pred_cls == 0 & y_true == 0)
  fp <- sum(pred_cls == 1 & y_true == 0)
  fn <- sum(pred_cls == 0 & y_true == 1)
  sens <- tp / (tp + fn)
  spec <- tn / (tn + fp)
  prec <- ifelse((tp + fp) > 0, tp / (tp + fp), NA)
  f1   <- ifelse(!is.na(prec) && (prec + sens) > 0, 2 * prec * sens / (prec + sens), NA)
  acc  <- (tp + tn) / length(y_true)
  bacc <- (sens + spec) / 2

  list(
    roc = roc_o,
    summary = data.frame(
      Model             = model_name,
      Evaluation        = eval_type,
      AUC               = round(as.numeric(auc(roc_o)), 4),
      CI_Lower          = round(as.numeric(roc_o$ci[1]), 4),
      CI_Upper          = round(as.numeric(roc_o$ci[3]), 4),
      Accuracy          = round(acc, 4),
      Balanced_Accuracy = round(bacc, 4),
      Sensitivity       = round(sens, 4),
      Specificity       = round(spec, 4),
      Precision         = round(prec, 4),
      F1_Score          = round(f1, 4),
      stringsAsFactors  = FALSE
    )
  )
}

tune_predict_rf <- function(X_tr_sc, y_tr_fac, y_tr_bin, X_te_sc, rf_grid, seed = 42L) {
  inner_folds <- make_stratified_folds(y_tr_bin, k = 5L, seed = seed + 100L)
  df_tr_rf    <- data.frame(X_tr_sc, Class = y_tr_fac, check.names = TRUE)
  df_te_rf    <- data.frame(X_te_sc, check.names = TRUE)

  best_auc <- -1
  best_p   <- list(mtry = 2, nodesize = 5)
  for (g in seq_len(nrow(rf_grid))) {
    in_probs <- numeric(nrow(df_tr_rf))
    for (inf in 1:5) {
      itr <- which(inner_folds != inf)
      iva <- which(inner_folds == inf)
      set.seed(42)
      m_in <- randomForest(Class ~ ., data = df_tr_rf[itr, ],
                           mtry = rf_grid$mtry[g], nodesize = rf_grid$nodesize[g],
                           ntree = 300)
      in_probs[iva] <- predict(m_in, newdata = df_tr_rf[iva, ], type = "prob")[, "Immune-Hot"]
    }
    in_auc <- as.numeric(auc(roc(y_tr_bin, in_probs, levels = c(0, 1), direction = "<", quiet = TRUE)))
    if (in_auc > best_auc) {
      best_auc <- in_auc
      best_p   <- list(mtry = rf_grid$mtry[g], nodesize = rf_grid$nodesize[g])
    }
  }
  set.seed(seed)
  fit_rf <- randomForest(Class ~ ., data = df_tr_rf,
                         mtry = best_p$mtry, nodesize = best_p$nodesize, ntree = 1000)
  predict(fit_rf, newdata = df_te_rf, type = "prob")[, "Immune-Hot"]
}

tune_predict_svm <- function(X_tr_sc, y_tr_fac, y_tr_bin, X_te_sc, svm_grid, seed = 42L) {
  inner_folds <- make_stratified_folds(y_tr_bin, k = 5L, seed = seed + 200L)
  best_auc <- -1
  best_p   <- list(cost = 1, gamma = 0.1)
  for (g in seq_len(nrow(svm_grid))) {
    in_probs <- numeric(nrow(X_tr_sc))
    for (inf in 1:5) {
      itr <- which(inner_folds != inf)
      iva <- which(inner_folds == inf)
      m_in <- svm(x = X_tr_sc[itr, , drop = FALSE], y = y_tr_fac[itr],
                  kernel = "radial", cost = svm_grid$cost[g], gamma = svm_grid$gamma[g],
                  scale = FALSE, probability = TRUE)
      p_in <- predict(m_in, X_tr_sc[iva, , drop = FALSE], probability = TRUE)
      in_probs[iva] <- attr(p_in, "probabilities")[, "Immune-Hot"]
    }
    in_auc <- as.numeric(auc(roc(y_tr_bin, in_probs, levels = c(0, 1), direction = "<", quiet = TRUE)))
    if (in_auc > best_auc) {
      best_auc <- in_auc
      best_p   <- list(cost = svm_grid$cost[g], gamma = svm_grid$gamma[g])
    }
  }
  fit_svm <- svm(x = X_tr_sc, y = y_tr_fac, kernel = "radial",
                 cost = best_p$cost, gamma = best_p$gamma,
                 scale = FALSE, probability = TRUE)
  p_svm <- predict(fit_svm, X_te_sc, probability = TRUE)
  attr(p_svm, "probabilities")[, "Immune-Hot"]
}

fit_predict_ann <- function(X_tr_sc, y_tr, X_te_sc, seed = 42L) {
  early_stop <- callback_early_stopping(
    monitor = "val_loss", patience = as.integer(15), restore_best_weights = TRUE
  )

  reticulate::py_run_string(sprintf(
    paste(
      "import tensorflow as tf, gc",
      "try:",
      "    tf.keras.backend.clear_session()",
      "except Exception:",
      "    pass",
      "try:",
      "    tf.keras.utils.set_random_seed(%d)",
      "except Exception:",
      "    tf.random.set_seed(%d)",
      "gc.collect()",
      sep = "\n"
    ),
    seed, seed
  ))

  set.seed(seed)
  shuf   <- sample(seq_len(nrow(X_tr_sc)))
  X_tr_s <- unname(as.matrix(X_tr_sc[shuf, , drop = FALSE]))
  storage.mode(X_tr_s) <- "double"
  y_tr_s <- as.array(as.numeric(y_tr[shuf]))

  X_te_s <- unname(as.matrix(X_te_sc))
  storage.mode(X_te_s) <- "double"

  model <- keras_model_sequential(layers = list(
    layer_input(shape = c(as.integer(ncol(X_tr_s)))),
    layer_dense(units = 32, activation = "relu"),
    layer_dropout(rate = 0.3),
    layer_dense(units = 16, activation = "relu"),
    layer_dropout(rate = 0.2),
    layer_dense(units = 1, activation = "sigmoid")
  ))
  model$compile(
    loss = "binary_crossentropy",
    optimizer = optimizer_adam(learning_rate = 0.001),
    metrics = list("accuracy")
  )
  hist <- model$fit(
    x = X_tr_s, y = y_tr_s,
    epochs = as.integer(100), batch_size = as.integer(16),
    validation_split = 0.2, callbacks = list(early_stop), verbose = 0
  )
  preds  <- as.numeric(model$predict(X_te_s, verbose = 0))
  hist_r <- lapply(reticulate::py_to_r(hist$history), as.numeric)
  list(preds = preds, history = hist_r)
}

# ==============================================================================
# PART A: STRICTLY NESTED 10-FOLD CV ON LOCKED 5-GENE PANEL
# ==============================================================================
df_app        <- data.frame(X_5, y = y_bin, check.names = TRUE)
glm_app       <- glm(y ~ ., data = df_app, family = binomial)
glm_app_probs <- predict(glm_app, type = "response")
res_glm_app   <- compute_metrics(y_bin, glm_app_probs, "Logistic Regression", "Apparent (In-Sample)")

outer_k     <- 10L
outer_folds <- make_stratified_folds(y_bin, k = outer_k, seed = 42L)

oof_glm <- numeric(nrow(X_5))
oof_rf  <- numeric(nrow(X_5))
oof_svm <- numeric(nrow(X_5))
oof_ann <- numeric(nrow(X_5))

rf_grid  <- expand.grid(mtry = c(1, 2, 3, 4, 5), nodesize = c(1, 3, 5, 10, 15))
svm_grid <- expand.grid(cost = c(0.1, 1, 10, 100), gamma = c(0.001, 0.01, 0.1, 1))

ann_fold_histories <- vector("list", outer_k)

for (f in seq_len(outer_k)) {
  te_idx <- which(outer_folds == f)
  tr_idx <- which(outer_folds != f)

  mu    <- colMeans(X_5[tr_idx, , drop = FALSE])
  sd_tr <- apply(X_5[tr_idx, , drop = FALSE], 2, sd)
  sd_tr[sd_tr == 0] <- 1
  X_tr_sc <- scale(X_5[tr_idx, , drop = FALSE], center = mu, scale = sd_tr)
  X_te_sc <- scale(X_5[te_idx, , drop = FALSE], center = mu, scale = sd_tr)

  y_tr_bin <- y_bin[tr_idx]
  y_tr_fac <- y_fac[tr_idx]

  fit_glm <- glm(y ~ ., data = data.frame(X_tr_sc, y = y_tr_bin, check.names = TRUE), family = binomial)
  oof_glm[te_idx] <- predict(fit_glm, newdata = data.frame(X_te_sc, check.names = TRUE), type = "response")

  oof_rf[te_idx]  <- tune_predict_rf(X_tr_sc, y_tr_fac, y_tr_bin, X_te_sc, rf_grid, seed = 42L + f)
  oof_svm[te_idx] <- tune_predict_svm(X_tr_sc, y_tr_fac, y_tr_bin, X_te_sc, svm_grid, seed = 42L + f)

  ann_out <- fit_predict_ann(X_tr_sc, y_tr_bin, X_te_sc, seed = 42L + f)
  oof_ann[te_idx] <- ann_out$preds
  ann_fold_histories[[f]] <- ann_out$history
}

res_glm_oof <- compute_metrics(y_bin, oof_glm, "Logistic Regression",  "Nested 10-Fold OOF")
res_rf_oof  <- compute_metrics(y_bin, oof_rf,  "Random Forest",        "Nested 10-Fold OOF")
res_svm_oof <- compute_metrics(y_bin, oof_svm, "SVM (Radial)",         "Nested 10-Fold OOF")
res_ann_oof <- compute_metrics(y_bin, oof_ann, "ANN (Early Stopping)", "Nested 10-Fold OOF")

partA_table <- rbind(
  res_glm_app$summary,
  res_glm_oof$summary,
  res_rf_oof$summary,
  res_svm_oof$summary,
  res_ann_oof$summary
)
write.csv(partA_table, "revision_outputs/Table_Step2A_NestedCV_5Gene.csv", row.names = FALSE)

# ==============================================================================
# PART B: END-TO-END LEAKAGE-FREE 5-FOLD NESTED CV
# ==============================================================================
svm_rfe_fast <- function(x, y) {
  feature_list    <- colnames(x)
  ranked_features <- character(0)
  while (length(feature_list) > 0) {
    svm_fit <- svm(x = x[, feature_list, drop = FALSE], y = y, kernel = "linear", scale = TRUE)
    w       <- t(svm_fit$coefs) %*% svm_fit$SV
    weights <- as.numeric(w[1, ])^2
    names(weights) <- feature_list
    n_drop <- if (length(feature_list) > 80) max(1L, floor(length(feature_list) * 0.1)) else 1L
    worst  <- names(sort(weights, decreasing = FALSE))[1:n_drop]
    ranked_features <- c(worst, ranked_features)
    feature_list    <- setdiff(feature_list, worst)
  }
  ranked_features
}

e2e_k     <- 5L
e2e_folds <- make_stratified_folds(y_bin, k = e2e_k, seed = 2026L)

e2e_oof_glm <- numeric(nrow(X_all))
e2e_oof_rf  <- numeric(nrow(X_all))
e2e_oof_svm <- numeric(nrow(X_all))

fold_selected_genes <- vector("list", e2e_k)
fold_lasso_genes    <- vector("list", e2e_k)
fold_rf_genes       <- vector("list", e2e_k)
fold_svm_genes      <- vector("list", e2e_k)

for (f in seq_len(e2e_k)) {
  te_idx <- which(e2e_folds == f)
  tr_idx <- which(e2e_folds != f)

  X_tr_full <- X_all[tr_idx, , drop = FALSE]
  y_tr_fac  <- y_fac[tr_idx]
  y_tr_bin  <- y_bin[tr_idx]
  y_tr_cold <- ifelse(y_tr_fac == "Immune-Cold", 1, 0)

  g_vars <- apply(X_tr_full, 2, var, na.rm = TRUE)
  top5k  <- names(sort(g_vars, decreasing = TRUE))[1:5000]
  dat_tr <- X_tr_full[, top5k, drop = FALSE]

  cor <- WGCNA::cor
  bnet <- blockwiseModules(
    dat_tr, power = 9, networkType = "signed", TOMType = "signed",
    minModuleSize = 30, deepSplit = 2, pamRespectsDendro = FALSE,
    mergeCutHeight = 0.25, maxBlockSize = 5000, numericLabels = FALSE, verbose = 0
  )
  cor <- stats::cor

  MEs_tr     <- bnet$MEs[, colnames(bnet$MEs) != "MEgrey", drop = FALSE]
  mod_cor    <- cor(MEs_tr, y_tr_cold, use = "p")
  top_mod    <- substring(rownames(mod_cor)[which.max(abs(mod_cor[, 1]))], 3)
  cand_genes <- colnames(dat_tr)[bnet$colors == top_mod]

  X_tr_cand <- X_tr_full[, cand_genes, drop = FALSE]

  set.seed(123L + f)
  cv_l    <- cv.glmnet(x = as.matrix(X_tr_cand), y = y_tr_cold, family = "binomial", alpha = 1, nfolds = 10)
  lc      <- coef(cv_l, s = "lambda.min")
  g_lasso <- setdiff(rownames(lc)[which(lc[, 1] != 0)], "(Intercept)")

  set.seed(123L + f)
  rf_fs <- randomForest(x = X_tr_cand, y = y_tr_fac, ntree = 1000, importance = TRUE)
  g_rf  <- names(sort(importance(rf_fs)[, "MeanDecreaseGini"], decreasing = TRUE))[1:min(30L, ncol(X_tr_cand))]

  g_svm <- svm_rfe_fast(X_tr_cand, y_tr_fac)[1:min(30L, ncol(X_tr_cand))]

  g_inter <- Reduce(intersect, list(g_lasso, g_rf, g_svm))
  if (length(g_inter) == 0) {
    all_counts <- table(c(g_lasso, g_rf, g_svm))
    g_inter    <- names(all_counts[all_counts >= 2])
  }

  fold_lasso_genes[[f]]    <- g_lasso
  fold_rf_genes[[f]]       <- g_rf
  fold_svm_genes[[f]]      <- g_svm
  fold_selected_genes[[f]] <- g_inter

  mu    <- colMeans(X_all[tr_idx, g_inter, drop = FALSE])
  sd_tr <- apply(X_all[tr_idx, g_inter, drop = FALSE], 2, sd)
  sd_tr[sd_tr == 0] <- 1
  X_tr_sel <- scale(X_all[tr_idx, g_inter, drop = FALSE], center = mu, scale = sd_tr)
  X_te_sel <- scale(X_all[te_idx, g_inter, drop = FALSE], center = mu, scale = sd_tr)

  df_tr_b   <- data.frame(X_tr_sel, check.names = TRUE)
  df_tr_b$y <- y_tr_bin
  df_te_b   <- data.frame(X_te_sel, check.names = TRUE)
  f_glm     <- glm(y ~ ., data = df_tr_b, family = binomial)
  e2e_oof_glm[te_idx] <- predict(f_glm, newdata = df_te_b, type = "response")

  set.seed(42L + f)
  f_rf <- randomForest(x = X_tr_sel, y = y_tr_fac, ntree = 1000)
  e2e_oof_rf[te_idx] <- predict(f_rf, newdata = X_te_sel, type = "prob")[, "Immune-Hot"]

  f_svm <- svm(x = X_tr_sel, y = y_tr_fac, kernel = "radial", scale = FALSE, probability = TRUE)
  p_svm <- predict(f_svm, X_te_sel, probability = TRUE)
  e2e_oof_svm[te_idx] <- attr(p_svm, "probabilities")[, "Immune-Hot"]
}

res_e2e_glm <- compute_metrics(y_bin, e2e_oof_glm, "End-to-End GLM", "5-Fold Full Nested CV")
res_e2e_rf  <- compute_metrics(y_bin, e2e_oof_rf,  "End-to-End RF",  "5-Fold Full Nested CV")
res_e2e_svm <- compute_metrics(y_bin, e2e_oof_svm, "End-to-End SVM", "5-Fold Full Nested CV")

partB_table <- rbind(res_e2e_glm$summary, res_e2e_rf$summary, res_e2e_svm$summary)
write.csv(partB_table, "revision_outputs/Table_Step2B_EndToEnd_NestedCV.csv", row.names = FALSE)

stability_df <- data.frame(
  Gene           = hub_genes,
  LASSO_Freq     = sapply(hub_genes, function(g) sum(sapply(fold_lasso_genes,    function(v) g %in% v))),
  RF_Freq        = sapply(hub_genes, function(g) sum(sapply(fold_rf_genes,       function(v) g %in% v))),
  SVM_RFE_Freq   = sapply(hub_genes, function(g) sum(sapply(fold_svm_genes,      function(v) g %in% v))),
  Consensus_Freq = sapply(hub_genes, function(g) sum(sapply(fold_selected_genes, function(v) g %in% v))),
  Total_Folds    = e2e_k
)
write.csv(stability_df, "revision_outputs/Table_Step2B_Feature_Stability.csv", row.names = FALSE)

# ==============================================================================
# PART C: LEAVE-ONE-COHORT-OUT (LOCO) VALIDATION (PRE-COMBAT fRMA)
# ==============================================================================
load("Phase2_Mapped_Data_and_Metadata.RData")
stopifnot(all(pheno_immune$Sample %in% colnames(mapped_matrix)))

X_raw_5 <- t(mapped_matrix[hub_genes, pheno_immune$Sample, drop = FALSE])
rm(mapped_matrix, meta_31210, meta_50081, meta_37745, meta_30219)
gc()

X_raw_z <- X_raw_5
for (b in unique(batches)) {
  idx_b <- which(batches == b)
  X_raw_z[idx_b, ] <- scale(X_raw_5[idx_b, , drop = FALSE])
}

unique_cohorts <- sort(unique(batches))
loco_rows      <- list()
loco_oof_glm   <- numeric(nrow(X_raw_z))
loco_oof_rf    <- numeric(nrow(X_raw_z))
loco_oof_svm   <- numeric(nrow(X_raw_z))
loco_oof_ann   <- numeric(nrow(X_raw_z))

for (i in seq_along(unique_cohorts)) {
  coh    <- unique_cohorts[i]
  te_idx <- which(batches == coh)
  tr_idx <- which(batches != coh)

  X_tr   <- X_raw_z[tr_idx, , drop = FALSE]
  X_te   <- X_raw_z[te_idx, , drop = FALSE]
  y_tr_b <- y_bin[tr_idx]
  y_te_b <- y_bin[te_idx]
  y_tr_f <- y_fac[tr_idx]

  m_glm <- glm(y ~ ., data = data.frame(X_tr, y = y_tr_b, check.names = TRUE), family = binomial)
  p_glm <- predict(m_glm, newdata = data.frame(X_te, check.names = TRUE), type = "response")
  loco_oof_glm[te_idx] <- p_glm

  p_rf  <- tune_predict_rf(X_tr, y_tr_f, y_tr_b, X_te, rf_grid, seed = 300L + i)
  loco_oof_rf[te_idx] <- p_rf

  p_svm <- tune_predict_svm(X_tr, y_tr_f, y_tr_b, X_te, svm_grid, seed = 300L + i)
  loco_oof_svm[te_idx] <- p_svm

  p_ann <- fit_predict_ann(X_tr, y_tr_b, X_te, seed = 300L + i)$preds
  loco_oof_ann[te_idx] <- p_ann

  get_auc <- function(p) round(as.numeric(auc(roc(y_te_b, p, levels = c(0, 1), direction = "<", quiet = TRUE))), 4)
  loco_rows[[coh]] <- data.frame(
    HeldOut_Cohort   = coh,
    N_Total          = length(te_idx),
    N_Hot            = sum(y_te_b == 1),
    N_Cold           = sum(y_te_b == 0),
    GLM_AUC          = get_auc(p_glm),
    RF_AUC           = get_auc(p_rf),
    SVM_AUC          = get_auc(p_svm),
    ANN_AUC          = get_auc(p_ann),
    stringsAsFactors = FALSE
  )
}

loco_per_cohort <- do.call(rbind, loco_rows)
loco_pooled <- rbind(
  compute_metrics(y_bin, loco_oof_glm, "Logistic Regression",  "LOCO Pooled (Pre-ComBat)")$summary,
  compute_metrics(y_bin, loco_oof_rf,  "Random Forest",        "LOCO Pooled (Pre-ComBat)")$summary,
  compute_metrics(y_bin, loco_oof_svm, "SVM (Radial)",         "LOCO Pooled (Pre-ComBat)")$summary,
  compute_metrics(y_bin, loco_oof_ann, "ANN (Early Stopping)", "LOCO Pooled (Pre-ComBat)")$summary
)

write.csv(loco_per_cohort, "revision_outputs/Table_Step2C_LOCO_PerCohort.csv", row.names = FALSE)
write.csv(loco_pooled,     "revision_outputs/Table_Step2C_LOCO_Pooled.csv",    row.names = FALSE)

save(
  hub_genes, partA_table, res_glm_app, res_glm_oof, res_rf_oof, res_svm_oof, res_ann_oof,
  oof_glm, oof_rf, oof_svm, oof_ann, ann_fold_histories,
  partB_table, stability_df, fold_selected_genes, e2e_oof_glm, e2e_oof_rf, e2e_oof_svm,
  loco_per_cohort, loco_pooled, loco_oof_glm, loco_oof_rf, loco_oof_svm, loco_oof_ann,
  file = "revision_outputs/Step2_LeakageFree_Results.RData"
)
message("Saved revision_outputs/Step2_LeakageFree_Results.RData")