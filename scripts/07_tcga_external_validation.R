# ==============================================================================
# 07_tcga_external_validation.R
# External cross-platform validation in GDC TCGA-LUAD RNA-seq (STAR-TPM):
#   Part A: Independent LM22 nu-SVR deconvolution & zero-shot classifier ROC
#   Part B: Orthogonal TME & immunotherapy biomarker correlations (GEP-18, TLS, CYT)
#   Part C: Complete five-gene survival & pathologic stage Cox regression
#   Part D: High-resolution (600 DPI) external validation figures (3A, 3B, 3C)
# ==============================================================================

options(timeout = 600)

suppressPackageStartupMessages({
  library(data.table)
  library(e1071)
  library(randomForest)
  library(pROC)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(patchwork)
  library(parallel)
})

dir.create("revision_outputs/figures", showWarnings = FALSE, recursive = TRUE)

download_first_working <- function(urls, destfile) {
  if (file.exists(destfile) && file.info(destfile)$size > 10000) {
    return(destfile)
  }
  for (u in urls) {
    ok <- tryCatch({
      status <- download.file(u, destfile = destfile, method = "curl",
                              extra = "-L -f --retry 2 --connect-timeout 30", quiet = TRUE)
      status == 0 && file.exists(destfile) && file.info(destfile)$size > 10000
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (ok) return(destfile)
    if (file.exists(destfile)) unlink(destfile)
  }
  NULL
}

# 1. Load locked GEO discovery data and retrieve TCGA-LUAD expression + clinical data
load("Phase3_Immune_Phenotypes.RData")
load("Phase5_Final_Hub_Genes.RData")

clin_cache <- "revision_outputs/TCGA_LUAD_Clinical.RData"
if (file.exists(clin_cache)) {
  load(clin_cache)
} else if (file.exists("Phase6_TCGA_Validation_Cohort.RData")) {
  load("Phase6_TCGA_Validation_Cohort.RData")
  save(tcga_clinical, file = clin_cache)
} else {
  clin_file <- download_first_working(
    urls = c(
      "https://tcga.xenahubs.net/download/TCGA.LUAD.sampleMap/LUAD_clinicalMatrix",
      "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap/LUAD_clinicalMatrix"
    ),
    destfile = "revision_outputs/LUAD_clinicalMatrix"
  )
  stopifnot("Failed to download TCGA-LUAD clinical matrix." = !is.null(clin_file))
  tcga_clinical <- as.data.frame(fread(clin_file))
  save(tcga_clinical, file = clin_cache)
}

tcga_cache <- "revision_outputs/TCGA_LUAD_GDC_STAR_TPM.RData"
if (!file.exists(tcga_cache)) {
  tpm_file <- download_first_working(
    urls = c(
      "https://gdc.xenahubs.net/download/TCGA-LUAD.star_tpm.tsv.gz",
      "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.star_tpm.tsv.gz"
    ),
    destfile = "revision_outputs/TCGA-LUAD.star_tpm.tsv.gz"
  )
  stopifnot("Failed to download GDC TCGA-LUAD STAR-TPM matrix." = !is.null(tpm_file))

  dt_tpm <- fread(tpm_file)

  pm_file <- download_first_working(
    urls = c(
      "https://gdc.xenahubs.net/download/gencode.v36.annotation.gtf.gene.probemap",
      "https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap",
      "https://gdc.xenahubs.net/download/gencode.v22.annotation.gtf.gene.probemap"
    ),
    destfile = "revision_outputs/gencode.annotation.gtf.gene.probemap"
  )

  # Retain primary solid tumors ("01") and deduplicate to one sample per patient
  sample_cols <- colnames(dt_tpm)[-1]
  prim_cols   <- sample_cols[substr(sample_cols, 14, 15) == "01"]
  pat_ids     <- substr(prim_cols, 1, 12)
  keep_idx    <- !duplicated(pat_ids)
  prim_cols   <- prim_cols[keep_idx]
  pat_ids     <- pat_ids[keep_idx]

  ens_ids <- dt_tpm[[1]]
  if (!is.null(pm_file)) {
    pm      <- fread(pm_file)
    sym_map <- pm$gene[match(ens_ids, pm$id)]
    na_m    <- is.na(sym_map)
    if (any(na_m)) {
      sym_map[na_m] <- pm$gene[match(sub("\\..*$", "", ens_ids[na_m]), sub("\\..*$", "", pm$id))]
    }
  } else {
    suppressPackageStartupMessages({
      library(AnnotationDbi)
      library(org.Hs.eg.db)
    })
    sym_map <- unname(mapIds(org.Hs.eg.db, keys = sub("\\..*$", "", ens_ids),
                             column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
  }

  dt_sub <- dt_tpm[, c(colnames(dt_tpm)[1], prim_cols), with = FALSE]
  dt_sub[, Gene := sym_map]
  dt_sub <- dt_sub[!is.na(Gene) & nzchar(Gene)]
  dt_sub[, (colnames(dt_tpm)[1]) := NULL]

  dt_col        <- dt_sub[, lapply(.SD, max), by = Gene]
  tcga_gdc_log2 <- as.matrix(dt_col[, -1, with = FALSE])
  rownames(tcga_gdc_log2) <- dt_col$Gene
  colnames(tcga_gdc_log2) <- pat_ids
  tcga_barcodes <- setNames(prim_cols, pat_ids)

  save(tcga_gdc_log2, tcga_barcodes, file = tcga_cache)
  rm(dt_tpm, dt_sub, dt_col)
  gc()
} else {
  load(tcga_cache)
}

stopifnot(all(hub_genes %in% rownames(tcga_gdc_log2)))

# 2. Independent TCGA-LUAD deconvolution (local CIBERSORT LM22 nu-SVR)
ciber_cache <- "revision_outputs/TCGA_LUAD_LM22_CIBERSORT.RData"
if (!file.exists(ciber_cache)) {
  lm22        <- as.matrix(read.delim("LM22.txt", row.names = 1, check.names = FALSE))
  shared_lm22 <- intersect(rownames(lm22), rownames(tcga_gdc_log2))

  is_log   <- max(tcga_gdc_log2, na.rm = TRUE) < 40
  Y_lin    <- if (is_log) pmax(2^tcga_gdc_log2[shared_lm22, , drop = FALSE] - 1, 0) else tcga_gdc_log2[shared_lm22, , drop = FALSE]
  X_sig    <- lm22[shared_lm22, , drop = FALSE]
  X_sig_sc <- (X_sig - mean(X_sig)) / sd(as.numeric(X_sig))

  cibersort_one <- function(y_vec, X_mat) {
    y_sc <- as.numeric(scale(y_vec))
    y_sc[!is.finite(y_sc)] <- 0
    best_rmse <- Inf
    best_w    <- rep(0, ncol(X_mat))
    for (nu_val in c(0.25, 0.5, 0.75)) {
      fit <- e1071::svm(x = X_mat, y = y_sc, type = "nu-regression",
                        kernel = "linear", nu = nu_val, scale = FALSE)
      w <- as.numeric(t(fit$coefs) %*% fit$SV)
      w[w < 0] <- 0
      if (sum(w) > 0) w <- w / sum(w)
      rmse <- sqrt(mean((as.numeric(X_mat %*% w) - y_sc)^2))
      if (rmse < best_rmse) {
        best_rmse <- rmse
        best_w    <- w
      }
    }
    best_w
  }

  n_cores <- max(1L, min(10L, detectCores() - 1L))
  cl <- makeCluster(n_cores)
  clusterExport(cl, varlist = c("X_sig_sc", "cibersort_one"), envir = environment())
  clusterEvalQ(cl, library(e1071))

  sample_list <- lapply(seq_len(ncol(Y_lin)), function(j) Y_lin[, j])
  w_list      <- parLapply(cl, sample_list, function(yv) cibersort_one(yv, X_sig_sc))
  stopCluster(cl)

  tcga_lm22 <- do.call(rbind, w_list)
  rownames(tcga_lm22) <- colnames(Y_lin)
  colnames(tcga_lm22) <- colnames(X_sig)

  save(tcga_lm22, shared_lm22, file = ciber_cache)
} else {
  load(ciber_cache)
}

tcga_lm22_z <- scale(tcga_lm22)
tcga_lm22_z[!is.finite(tcga_lm22_z)] <- 0
hc_tcga        <- hclust(dist(tcga_lm22_z, method = "euclidean"), method = "ward.D2")
cl_tcga        <- cutree(hc_tcga, k = 2)
cd8_means_tcga <- tapply(tcga_lm22_z[, "T cells CD8"], cl_tcga, mean)
hot_id_tcga    <- as.integer(names(which.max(cd8_means_tcga)))

tcga_pheno <- factor(ifelse(cl_tcga == hot_id_tcga, "Immune-Hot", "Immune-Cold"),
                     levels = c("Immune-Cold", "Immune-Hot"))
names(tcga_pheno) <- rownames(tcga_lm22)
y_tcga_bin        <- as.integer(tcga_pheno == "Immune-Hot")

# ==============================================================================
# PART A: ZERO-SHOT EXTERNAL VALIDATION OF LOCKED GEO CLASSIFIERS
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

X_geo_5  <- scale(t(harmonized_matrix[hub_genes, , drop = FALSE]))
y_geo_f  <- factor(pheno_immune$Immune_Phenotype, levels = c("Immune-Cold", "Immune-Hot"))
y_geo_b  <- as.integer(y_geo_f == "Immune-Hot")
X_tcga_5 <- scale(t(tcga_gdc_log2[hub_genes, rownames(tcga_lm22), drop = FALSE]))

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
      Model = model_name, Evaluation = eval_type,
      AUC = round(as.numeric(auc(roc_o)), 4),
      CI_Lower = round(as.numeric(roc_o$ci[1]), 4),
      CI_Upper = round(as.numeric(roc_o$ci[3]), 4),
      Accuracy = round(acc, 4), Balanced_Accuracy = round(bacc, 4),
      Sensitivity = round(sens, 4), Specificity = round(spec, 4),
      Precision = round(prec, 4), F1_Score = round(f1, 4),
      stringsAsFactors = FALSE
    )
  )
}

make_stratified_folds <- function(y, k = 5L, seed = 42L) {
  set.seed(seed)
  folds <- integer(length(y))
  for (cls in unique(y)) {
    idx <- which(y == cls)
    folds[idx] <- sample(rep(1:k, length.out = length(idx)))
  }
  folds
}
geo_folds <- make_stratified_folds(y_geo_b, k = 5L, seed = 42L)

# 1. Locked GEO GLM
df_geo_glm     <- data.frame(X_geo_5, y = y_geo_b, check.names = TRUE)
df_tcga_glm    <- data.frame(X_tcga_5, check.names = TRUE)
locked_glm     <- glm(y ~ ., data = df_geo_glm, family = binomial)
tcga_prob_glm  <- predict(locked_glm, newdata = df_tcga_glm, type = "response")
tcga_score_hot <- as.numeric(predict(locked_glm, newdata = df_tcga_glm, type = "link"))

# 2. Locked GEO Random Forest (tuned strictly on GEO via 5-fold CV)
rf_grid     <- expand.grid(mtry = c(1, 2, 3, 4, 5), nodesize = c(1, 3, 5, 10, 15))
df_geo_rf   <- data.frame(X_geo_5, Class = y_geo_f, check.names = TRUE)
best_rf_auc <- -1
best_rf_p   <- list(mtry = 2, nodesize = 5)
for (g in seq_len(nrow(rf_grid))) {
  in_p <- numeric(nrow(df_geo_rf))
  for (f in 1:5) {
    set.seed(42)
    m_in <- randomForest(Class ~ ., data = df_geo_rf[geo_folds != f, ],
                         mtry = rf_grid$mtry[g], nodesize = rf_grid$nodesize[g], ntree = 300)
    in_p[geo_folds == f] <- predict(m_in, newdata = df_geo_rf[geo_folds == f, ], type = "prob")[, "Immune-Hot"]
  }
  a <- as.numeric(auc(roc(y_geo_b, in_p, levels = c(0, 1), direction = "<", quiet = TRUE)))
  if (a > best_rf_auc) {
    best_rf_auc <- a
    best_rf_p   <- list(mtry = rf_grid$mtry[g], nodesize = rf_grid$nodesize[g])
  }
}
set.seed(42)
locked_rf    <- randomForest(Class ~ ., data = df_geo_rf,
                             mtry = best_rf_p$mtry, nodesize = best_rf_p$nodesize, ntree = 1000)
tcga_prob_rf <- predict(locked_rf, newdata = df_tcga_glm, type = "prob")[, "Immune-Hot"]

# 3. Locked GEO SVM-Radial (tuned strictly on GEO via 5-fold CV)
svm_grid     <- expand.grid(cost = c(0.1, 1, 10, 100), gamma = c(0.001, 0.01, 0.1, 1))
best_svm_auc <- -1
best_svm_p   <- list(cost = 1, gamma = 0.1)
for (g in seq_len(nrow(svm_grid))) {
  in_p <- numeric(nrow(X_geo_5))
  for (f in 1:5) {
    m_in <- svm(x = X_geo_5[geo_folds != f, , drop = FALSE], y = y_geo_f[geo_folds != f],
                kernel = "radial", cost = svm_grid$cost[g], gamma = svm_grid$gamma[g],
                scale = FALSE, probability = TRUE)
    p_in <- predict(m_in, X_geo_5[geo_folds == f, , drop = FALSE], probability = TRUE)
    in_p[geo_folds == f] <- attr(p_in, "probabilities")[, "Immune-Hot"]
  }
  a <- as.numeric(auc(roc(y_geo_b, in_p, levels = c(0, 1), direction = "<", quiet = TRUE)))
  if (a > best_svm_auc) {
    best_svm_auc <- a
    best_svm_p   <- list(cost = svm_grid$cost[g], gamma = svm_grid$gamma[g])
  }
}
locked_svm    <- svm(x = X_geo_5, y = y_geo_f, kernel = "radial",
                     cost = best_svm_p$cost, gamma = best_svm_p$gamma,
                     scale = FALSE, probability = TRUE)
tcga_prob_svm <- attr(predict(locked_svm, X_tcga_5, probability = TRUE), "probabilities")[, "Immune-Hot"]

# 4. Locked GEO ANN (5-seed ensemble trained strictly on shuffled GEO data)
ann_preds_mat <- matrix(0, nrow = nrow(X_tcga_5), ncol = 5)
for (s in 1:5) {
  seed_i     <- 40L + s
  early_stop <- callback_early_stopping(monitor = "val_loss", patience = as.integer(15), restore_best_weights = TRUE)
  reticulate::py_run_string(sprintf(
    "import tensorflow as tf, gc\ntry: tf.keras.backend.clear_session()\nexcept Exception: pass\ntry: tf.keras.utils.set_random_seed(%d)\nexcept Exception: tf.random.set_seed(%d)\ngc.collect()",
    seed_i, seed_i
  ))
  set.seed(seed_i)
  shuf   <- sample(seq_len(nrow(X_geo_5)))
  X_tr_a <- unname(as.matrix(X_geo_5[shuf, , drop = FALSE])); storage.mode(X_tr_a) <- "double"
  y_tr_a <- as.array(as.numeric(y_geo_b[shuf]))
  X_te_a <- unname(as.matrix(X_tcga_5)); storage.mode(X_te_a) <- "double"

  m_ann <- keras_model_sequential(layers = list(
    layer_input(shape = c(as.integer(ncol(X_tr_a)))),
    layer_dense(units = 32, activation = "relu"),
    layer_dropout(rate = 0.3),
    layer_dense(units = 16, activation = "relu"),
    layer_dropout(rate = 0.2),
    layer_dense(units = 1, activation = "sigmoid")
  ))
  m_ann$compile(loss = "binary_crossentropy", optimizer = optimizer_adam(learning_rate = 0.001), metrics = list("accuracy"))
  m_ann$fit(x = X_tr_a, y = y_tr_a, epochs = as.integer(100), batch_size = as.integer(16),
            validation_split = 0.2, callbacks = list(early_stop), verbose = 0)
  ann_preds_mat[, s] <- as.numeric(m_ann$predict(X_te_a, verbose = 0))
}
tcga_prob_ann <- rowMeans(ann_preds_mat)

ext_glm <- compute_metrics(y_tcga_bin, tcga_prob_glm, "Logistic Regression",  "TCGA-LUAD External")
ext_rf  <- compute_metrics(y_tcga_bin, tcga_prob_rf,  "Random Forest",        "TCGA-LUAD External")
ext_svm <- compute_metrics(y_tcga_bin, tcga_prob_svm, "SVM (Radial)",         "TCGA-LUAD External")
ext_ann <- compute_metrics(y_tcga_bin, tcga_prob_ann, "ANN (Early Stopping)", "TCGA-LUAD External")

partA_tcga_table <- rbind(ext_glm$summary, ext_rf$summary, ext_svm$summary, ext_ann$summary)
write.csv(partA_tcga_table, "revision_outputs/Table_Step3A_TCGA_External_Classification.csv", row.names = FALSE)

# ==============================================================================
# PART B: ORTHOGONAL IMMUNOLOGICAL & CLINICAL BIOMARKER VALIDATION
# ==============================================================================
tcga_expr <- tcga_gdc_log2[, rownames(tcga_lm22), drop = FALSE]

mean_z_sig <- function(expr_mat, genes) {
  g <- intersect(genes, rownames(expr_mat))
  rowMeans(scale(t(expr_mat[g, , drop = FALSE])))
}

cyt_score   <- colMeans(tcga_expr[intersect(c("GZMA", "PRF1"), rownames(tcga_expr)), , drop = FALSE])
gep18_genes <- c("CD3D", "IDO1", "CIITA", "CD3E", "CCL5", "GZMK", "CD2", "HLA-DRA", "CXCL13",
                 "IL2RG", "NKG7", "HLA-E", "CXCR6", "LAG3", "TAGAP", "CXCL10", "STAT1", "GZMB")
gep18_score <- mean_z_sig(tcga_expr, gep18_genes)
tls_genes   <- c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21",
                 "CXCL9", "CXCL10", "CXCL11", "CXCL13", "MS4A1", "LAMP3")
tls_score   <- mean_z_sig(tcga_expr, tls_genes)
cd8_eff     <- mean_z_sig(tcga_expr, c("CD8A", "CD8B", "IFNG"))

gep_high     <- as.integer(gep18_score >= median(gep18_score))
cyt_high     <- as.integer(cyt_score   >= median(cyt_score))
roc_gep_orth <- roc(gep_high, tcga_prob_glm, levels = c(0, 1), direction = "<", ci = TRUE, quiet = TRUE)
roc_cyt_orth <- roc(cyt_high, tcga_prob_glm, levels = c(0, 1), direction = "<", ci = TRUE, quiet = TRUE)

markers_list <- list(
  "T-cell Inflamed GEP (18-gene)" = gep18_score,
  "Mature TLS Signature"          = tls_score,
  "Cytolytic Activity (CYT)"      = cyt_score,
  "CD8+ Effector Score"           = cd8_eff,
  "PDCD1 (PD-1)"                  = as.numeric(tcga_expr["PDCD1", ]),
  "CTLA4"                         = as.numeric(tcga_expr["CTLA4", ]),
  "TIGIT"                         = as.numeric(tcga_expr["TIGIT", ]),
  "LAG3"                          = as.numeric(tcga_expr["LAG3", ]),
  "CD274 (PD-L1)"                 = as.numeric(tcga_expr["CD274", ]),
  "CD8A Expression"               = as.numeric(tcga_expr["CD8A", ]),
  "LM22 CD8+ T-cell Fraction"     = as.numeric(tcga_lm22[, "T cells CD8"]),
  "LM22 Plasma Cell Fraction"     = as.numeric(tcga_lm22[, "Plasma cells"])
)

pred_hot_lbl <- factor(ifelse(tcga_prob_glm >= 0.5, "Predicted Hot", "Predicted Cold"),
                       levels = c("Predicted Cold", "Predicted Hot"))

cor_table <- do.call(rbind, lapply(names(markers_list), function(nm) {
  v  <- markers_list[[nm]]
  ct <- suppressWarnings(cor.test(tcga_score_hot, v, method = "spearman"))
  wt <- wilcox.test(v ~ pred_hot_lbl)
  data.frame(
    Biomarker        = nm,
    Spearman_Rho     = round(as.numeric(ct$estimate), 4),
    Spearman_P       = signif(ct$p.value, 3),
    Median_Pred_Cold = round(median(v[pred_hot_lbl == "Predicted Cold"], na.rm = TRUE), 3),
    Median_Pred_Hot  = round(median(v[pred_hot_lbl == "Predicted Hot"],  na.rm = TRUE), 3),
    Wilcoxon_P       = signif(wt$p.value, 3),
    stringsAsFactors = FALSE
  )
}))
cor_table$Spearman_FDR <- signif(p.adjust(cor_table$Spearman_P, method = "BH"), 3)
cor_table$Wilcoxon_FDR <- signif(p.adjust(cor_table$Wilcoxon_P, method = "BH"), 3)
write.csv(cor_table, "revision_outputs/Table_Step3B_TCGA_Immune_Correlations.csv", row.names = FALSE)

# ==============================================================================
# PART C: FULL FIVE-GENE SURVIVAL & PATHOLOGIC STAGE ANALYSIS
# ==============================================================================
pat_tcga  <- rownames(tcga_lm22)
clin_prim <- tcga_clinical[substr(tcga_clinical$sampleID, 14, 15) == "01", ]
idx_m     <- match(pat_tcga, substr(clin_prim$sampleID, 1, 12))
clin_m    <- clin_prim[idx_m, ]
miss_c    <- is.na(idx_m)
if (any(miss_c)) {
  clin_m[miss_c, ] <- tcga_clinical[match(pat_tcga[miss_c], substr(tcga_clinical$sampleID, 1, 12)), ]
}

os_days   <- ifelse(!is.na(clin_m$days_to_death),
                    as.numeric(clin_m$days_to_death),
                    as.numeric(clin_m$days_to_last_followup))
os_stat   <- ifelse(clin_m$vital_status %in% c("DECEASED", "Dead", "dead", 1), 1, 0)
raw_stage <- trimws(as.character(clin_m$pathologic_stage))
raw_stage[raw_stage %in% c("", "[Discrepancy]", "Not Available", "NA", NA)] <- NA
stage_clean <- factor(sub("[AB]$", "", raw_stage),
                      levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))

surv_master <- data.frame(
  Patient_ID = pat_tcga,
  OS_Time    = os_days / 30.44,
  OS_Status  = os_stat,
  Hot_Score  = tcga_score_hot,
  Cold_Risk  = -tcga_score_hot,
  Stage      = stage_clean,
  stringsAsFactors = FALSE
)

stage_df <- surv_master[!is.na(surv_master$Stage), ]
kw_stage <- kruskal.test(Hot_Score ~ Stage, data = stage_df)

cox_data <- surv_master[!is.na(surv_master$OS_Time) & surv_master$OS_Time > 0 &
                        !is.na(surv_master$Cold_Risk) & !is.na(surv_master$Stage), ]

med_cut <- median(cox_data$Cold_Risk)
cox_data$Risk_Group <- factor(
  ifelse(cox_data$Cold_Risk >= med_cut, "Immune-Cold-like (High Risk)", "Immune-Hot-like (Low Risk)"),
  levels = c("Immune-Hot-like (Low Risk)", "Immune-Cold-like (High Risk)")
)

cox_uni        <- coxph(Surv(OS_Time, OS_Status) ~ Risk_Group, data = cox_data)
cox_multi      <- coxph(Surv(OS_Time, OS_Status) ~ Risk_Group + Stage, data = cox_data)
cox_multi_cont <- coxph(Surv(OS_Time, OS_Status) ~ Cold_Risk + Stage, data = cox_data)
zph_multi      <- cox.zph(cox_multi)

format_cox <- function(fit, label) {
  s <- summary(fit)
  data.frame(
    Model    = label,
    Term     = rownames(s$coefficients),
    HR       = round(s$conf.int[, "exp(coef)"], 3),
    CI_Lower = round(s$conf.int[, "lower .95"], 3),
    CI_Upper = round(s$conf.int[, "upper .95"], 3),
    Z_Stat   = round(s$coefficients[, "z"], 3),
    P_Value  = signif(s$coefficients[, "Pr(>|z|)"], 3),
    stringsAsFactors = FALSE
  )
}

cox_summary_table <- rbind(
  format_cox(cox_uni,        "Univariate (Median Split)"),
  format_cox(cox_multi,      "Multivariate (Median Split + Stage)"),
  format_cox(cox_multi_cont, "Multivariate (Continuous Score + Stage)")
)
write.csv(cox_summary_table, "revision_outputs/Table_Step3C_TCGA_Cox_Summary.csv", row.names = FALSE)

save(
  hub_genes, tcga_lm22, tcga_pheno, y_tcga_bin,
  locked_glm, locked_rf, locked_svm,
  tcga_prob_glm, tcga_prob_rf, tcga_prob_svm, tcga_prob_ann, tcga_score_hot,
  ext_glm, ext_rf, ext_svm, ext_ann, partA_tcga_table,
  roc_gep_orth, roc_cyt_orth, markers_list, cor_table,
  kw_stage, cox_data, cox_uni, cox_multi, cox_multi_cont, zph_multi, cox_summary_table,
  file = "revision_outputs/Step3_TCGA_External_Results.RData"
)

# ==============================================================================
# PART D: PUBLICATION-GRADE 600 DPI FIGURES (3A, 3B, 3C)
# ==============================================================================
theme_pub <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title        = element_text(face = "bold", size = base_size + 1, hjust = 0.5),
      plot.subtitle     = element_text(size = base_size - 1, hjust = 0.5, color = "gray30"),
      axis.title        = element_text(face = "bold", size = base_size),
      axis.text         = element_text(size = base_size - 1, color = "black"),
      legend.title      = element_text(face = "bold", size = base_size - 1),
      legend.text       = element_text(size = base_size - 2),
      legend.background = element_rect(fill = ggplot2::alpha("white", 0.85), color = "gray80"),
      plot.margin       = ggplot2::margin(10, 12, 10, 10)
    )
}

roc_to_df <- function(roc_obj, label) {
  df <- data.frame(FPR = 1 - roc_obj$specificities, TPR = roc_obj$sensitivities, Model = label)
  df[order(df$FPR, df$TPR), ]
}

# Figure 3A: External Classification ROC
l_g <- sprintf("GLM (AUC = %.4f [%.4f-%.4f])",           auc(ext_glm$roc), ext_glm$roc$ci[1], ext_glm$roc$ci[3])
l_s <- sprintf("SVM-Radial (AUC = %.4f [%.4f-%.4f])",    auc(ext_svm$roc), ext_svm$roc$ci[1], ext_svm$roc$ci[3])
l_a <- sprintf("ANN (AUC = %.4f [%.4f-%.4f])",           auc(ext_ann$roc), ext_ann$roc$ci[1], ext_ann$roc$ci[3])
l_r <- sprintf("Random Forest (AUC = %.4f [%.4f-%.4f])", auc(ext_rf$roc),  ext_rf$roc$ci[1],  ext_rf$roc$ci[3])
l_o <- sprintf("Orthogonal GEP-18 (AUC = %.4f)",         auc(roc_gep_orth))

df_3a <- rbind(
  roc_to_df(ext_glm$roc,  l_g),
  roc_to_df(ext_svm$roc,  l_s),
  roc_to_df(ext_ann$roc,  l_a),
  roc_to_df(ext_rf$roc,   l_r),
  roc_to_df(roc_gep_orth, l_o)
)
df_3a$Model <- factor(df_3a$Model, levels = c(l_g, l_s, l_a, l_r, l_o))

p_3a <- ggplot(df_3a, aes(x = FPR, y = TPR, color = Model, linetype = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray55", linewidth = 0.8) +
  geom_line(linewidth = 1.25) +
  scale_color_manual(values = setNames(c("#E64B35FF", "#3C5488FF", "#7E6148FF", "#00A087FF", "gray30"),
                                       levels(df_3a$Model))) +
  scale_linetype_manual(values = setNames(c("solid", "solid", "solid", "solid", "dashed"),
                                          levels(df_3a$Model))) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = sprintf("External Classifier Validation in TCGA-LUAD (n = %d)", length(y_tcga_bin)),
    subtitle = "Locked GEO-trained 5-gene models applied to independent RNA-seq cohort",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = NULL, linetype = NULL
  ) +
  theme_pub(14) +
  theme(legend.position = c(0.61, 0.20), legend.key.width = unit(1.3, "cm"))

ggsave("revision_outputs/figures/Figure_Step3A_TCGA_External_ROC.png",
       plot = p_3a, width = 7.6, height = 7.0, dpi = 600)

# Figure 3B: Orthogonal Immunological Biomarker Correlations
cor_plot_df <- cor_table
cor_plot_df$Biomarker <- factor(cor_plot_df$Biomarker,
                                levels = cor_plot_df$Biomarker[order(cor_plot_df$Spearman_Rho)])
cor_plot_df$Sig_Lbl <- ifelse(cor_plot_df$Spearman_FDR < 0.001,
                              sprintf("ρ = %.2f***", cor_plot_df$Spearman_Rho),
                              ifelse(cor_plot_df$Spearman_FDR < 0.01,
                                     sprintf("ρ = %.2f**", cor_plot_df$Spearman_Rho),
                                     ifelse(cor_plot_df$Spearman_FDR < 0.05,
                                            sprintf("ρ = %.2f*", cor_plot_df$Spearman_Rho),
                                            sprintf("ρ = %.2f", cor_plot_df$Spearman_Rho))))

sub_fdr <- if (max(cor_table$Spearman_FDR, na.rm = TRUE) < 0.001) {
  "Spearman ρ with locked 5-gene Hot Score (all BH FDR < 0.001)"
} else {
  "Spearman ρ with locked 5-gene Hot Score (*FDR<0.05, **<0.01, ***<0.001)"
}

p_3b_left <- ggplot(cor_plot_df, aes(x = Spearman_Rho, y = Biomarker, fill = Spearman_Rho)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.3) +
  geom_text(aes(label = Sig_Lbl), hjust = -0.10, size = 3.8, fontface = "bold") +
  scale_fill_gradient(low = "#9BB7D4", high = "#E64B35FF", guide = "none") +
  coord_cartesian(xlim = c(min(0, min(cor_plot_df$Spearman_Rho) - 0.05),
                           max(cor_plot_df$Spearman_Rho) + 0.22)) +
  labs(
    title = "A. Correlation with Independent TME Markers",
    subtitle = sub_fdr,
    x = "Spearman Rank Correlation (ρ)", y = NULL
  ) +
  theme_pub(13)

scat_df <- data.frame(Hot_Score = tcga_score_hot, GEP18 = gep18_score,
                      TLS = tls_score, Group = tcga_pheno)

p_3b_top <- ggplot(scat_df, aes(x = Hot_Score, y = GEP18, color = Group)) +
  geom_point(alpha = 0.55, size = 1.8) +
  geom_smooth(method = "lm", color = "black", linewidth = 1.0, se = TRUE) +
  scale_color_manual(values = c("Immune-Cold" = "#3C5488FF", "Immune-Hot" = "#E64B35FF")) +
  labs(
    title = sprintf("B. 18-Gene T-cell GEP (ρ = %.3f)",
                    cor_table$Spearman_Rho[cor_table$Biomarker == "T-cell Inflamed GEP (18-gene)"]),
    x = "5-Gene Immune-Hot Score", y = "Ayers 18-Gene GEP (Z)", color = "TCGA TME"
  ) +
  theme_pub(12) + theme(legend.position = "top")

p_3b_bot <- ggplot(scat_df, aes(x = Hot_Score, y = TLS, color = Group)) +
  geom_point(alpha = 0.55, size = 1.8) +
  geom_smooth(method = "lm", color = "black", linewidth = 1.0, se = TRUE) +
  scale_color_manual(values = c("Immune-Cold" = "#3C5488FF", "Immune-Hot" = "#E64B35FF")) +
  labs(
    title = sprintf("C. Mature TLS Signature (ρ = %.3f)",
                    cor_table$Spearman_Rho[cor_table$Biomarker == "Mature TLS Signature"]),
    x = "5-Gene Immune-Hot Score", y = "TLS Signature (Z)", color = NULL
  ) +
  theme_pub(12) + theme(legend.position = "none")

p_3b <- (p_3b_left | (p_3b_top / p_3b_bot)) + plot_layout(widths = c(1.38, 1))
ggsave("revision_outputs/figures/Figure_Step3B_TCGA_Immune_Correlations.png",
       plot = p_3b, width = 13.8, height = 6.8, dpi = 600)

# Figure 3C: Full 5-Gene Kaplan-Meier Survival Curve
uni_row   <- cox_summary_table[cox_summary_table$Model == "Univariate (Median Split)", ]
adj_row   <- cox_summary_table[cox_summary_table$Model == "Multivariate (Median Split + Stage)" &
                               grepl("Risk_Group", cox_summary_table$Term), ]
fit_km_3c <- survfit(Surv(OS_Time, OS_Status) ~ Risk_Group, data = cox_data)

annot_txt <- sprintf(
  "Univariate Cox:\n  HR = %.2f (95%% CI: %.2f–%.2f), p = %.4f\nStage-Adjusted Cox:\n  HR = %.2f (95%% CI: %.2f–%.2f), p = %.4f",
  uni_row$HR, uni_row$CI_Lower, uni_row$CI_Upper, uni_row$P_Value,
  adj_row$HR, adj_row$CI_Lower, adj_row$CI_Upper, adj_row$P_Value
)

png("revision_outputs/figures/Figure_Step3C_TCGA_Survival_Full5Gene.png",
    width = 8.8, height = 7.4, units = "in", res = 600)
p_3c <- ggsurvplot(
  fit_km_3c,
  data = cox_data,
  conf.int = TRUE,
  palette = c("#E64B35FF", "#3C5488FF"),
  title = sprintf("TCGA-LUAD Full 5-Gene Survival Validation (n = %d)", nrow(cox_data)),
  xlab = "Time (Months)",
  ylab = "Overall Survival Probability",
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.25,
  legend.title = "Median Split",
  legend.labs = c("Immune-Hot-like (Low Risk)", "Immune-Cold-like (High Risk)"),
  font.main = c(15, "bold"),
  font.x = c(14, "bold"),
  font.y = c(14, "bold"),
  font.tickslab = c(12.5, "plain"),
  font.legend = c(12.5, "plain"),
  risk.table.fontsize = 4.3,
  ggtheme = theme_classic(base_size = 13.5)
)
p_3c$plot <- p_3c$plot +
  annotate("label", x = 102, y = 0.87, label = annot_txt,
           hjust = 0, vjust = 0.5, size = 3.9, fontface = "bold",
           fill = ggplot2::alpha("white", 0.92), color = "black")
print(p_3c)
dev.off()

message("Saved revision_outputs/Step3_TCGA_External_Results.RData and Figures 3A-3C")