# ==============================================================================
# 02_annotation_and_combat.R
# Map GPL570 probes to HGNC gene symbols (collapsing duplicates by maximum
# expression), filter histologically confirmed LUAD tumors from GEO metadata,
# and apply empirical Bayes ComBat batch harmonization.
# ==============================================================================

suppressPackageStartupMessages({
  library(hgu133plus2.db)
  library(data.table)
  library(GEOquery)
  library(sva)
  library(ggplot2)
  library(patchwork)
  library(parallel)
})

setDTthreads(min(10L, max(1L, detectCores() - 1L)))

# 1. Load fRMA probe-level matrix and collapse probes to gene symbols
load("fRMA_Normalized_Data.RData")

map_probes_to_genes <- function(expr_mat) {
  mapped_probes <- mapIds(hgu133plus2.db,
                          keys = rownames(expr_mat),
                          column = "SYMBOL",
                          keytype = "PROBEID",
                          multiVals = "first")
  dt <- as.data.table(expr_mat, keep.rownames = "Probe")
  dt[, Gene := mapped_probes]
  dt <- dt[!is.na(Gene)]
  dt[, Probe := NULL]
  dt_collapsed <- dt[, lapply(.SD, max), by = Gene]

  collapsed_mat <- as.matrix(dt_collapsed[, -1, with = FALSE])
  rownames(collapsed_mat) <- dt_collapsed$Gene
  collapsed_mat
}

mapped_matrix <- map_probes_to_genes(combined_frma)
rm(combined_frma)
gc()

# 2. Retrieve clinical metadata from GEO
unlink(list.files(tempdir(), pattern = "GSE|GPL", full.names = TRUE))

get_geo_meta <- function(geo_id) {
  gse <- getGEO(geo_id, GSEMatrix = TRUE, getGPL = FALSE,
                destdir = tempdir(), returnType = "ExpressionSet")
  pData(gse[[1]])
}

meta_31210 <- get_geo_meta("GSE31210")
meta_50081 <- get_geo_meta("GSE50081")
meta_37745 <- get_geo_meta("GSE37745")
meta_30219 <- get_geo_meta("GSE30219")

# Retain pre-ComBat mapped matrix for independent LOCO validation
save(mapped_matrix, meta_31210, meta_50081, meta_37745, meta_30219,
     file = "Phase2_Mapped_Data_and_Metadata.RData")

# 3. Filter strictly to histologically confirmed LUAD samples
build_pheno <- function(meta, batch, filter_col, filter_pattern) {
  idx <- grepl(filter_pattern, meta[[filter_col]], ignore.case = TRUE)
  meta_sub <- meta[idx, ]

  supp_col <- grep("supplementary_file", colnames(meta_sub), value = TRUE)[1]
  supps_clean <- if (!is.na(supp_col)) {
    sub("\\.CEL.*$", "", basename(as.character(meta_sub[[supp_col]])), ignore.case = TRUE)
  } else {
    meta_sub$geo_accession
  }

  data.frame(
    GSM      = meta_sub$geo_accession,
    SuppFile = supps_clean,
    Batch    = batch,
    stringsAsFactors = FALSE
  )
}

target_pheno <- rbind(
  build_pheno(meta_31210, "GSE31210", "gene alteration status:ch1", "adenocarcinoma|kras|egfr"),
  build_pheno(meta_50081, "GSE50081", "histology:ch1", "adeno"),
  build_pheno(meta_37745, "GSE37745", "histology:ch1", "adeno"),
  build_pheno(meta_30219, "GSE30219", "histology:ch1", "adc|adeno")
)

clean_mat_cols <- sub("\\.CEL.*$", "", colnames(mapped_matrix), ignore.case = TRUE)
match_idx <- match(clean_mat_cols, target_pheno$SuppFile)
na_idx    <- is.na(match_idx)
match_idx[na_idx] <- match(clean_mat_cols[na_idx], target_pheno$GSM)

valid_cols <- which(!is.na(match_idx))
pheno <- target_pheno[match_idx[valid_cols], ]
pheno$Sample <- colnames(mapped_matrix)[valid_cols]

message("Mapped LUAD samples by cohort:")
print(table(pheno$Batch))

# 4. Subset expression matrix and remove zero-variance genes
filtered_matrix <- mapped_matrix[, pheno$Sample, drop = FALSE]
rm(mapped_matrix, meta_31210, meta_50081, meta_37745, meta_30219)
gc()

row_vars <- apply(filtered_matrix, 1, var, na.rm = TRUE)
clean_matrix <- filtered_matrix[is.finite(row_vars) & row_vars > 0, , drop = FALSE]

# 5. Execute ComBat empirical Bayes batch harmonization
mod <- model.matrix(~1, data = pheno)
harmonized_matrix <- ComBat(dat = clean_matrix, batch = pheno$Batch, mod = mod, par.prior = TRUE)

# 6. PCA quality-control check before vs. after ComBat
pca_raw  <- prcomp(t(clean_matrix), scale. = TRUE)
pca_norm <- prcomp(t(harmonized_matrix), scale. = TRUE)

df_raw  <- merge(data.frame(Sample = rownames(pca_raw$x),  PC1 = pca_raw$x[, 1],  PC2 = pca_raw$x[, 2]),  pheno, by = "Sample")
df_norm <- merge(data.frame(Sample = rownames(pca_norm$x), PC1 = pca_norm$x[, 1], PC2 = pca_norm$x[, 2]), pheno, by = "Sample")

p_before <- ggplot(df_raw, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 2, alpha = 0.7) + theme_classic(base_size = 12) +
  ggtitle("Before ComBat (Pre-Harmonization fRMA)") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))

p_after <- ggplot(df_norm, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 2, alpha = 0.7) + theme_classic(base_size = 12) +
  ggtitle("After ComBat (Harmonized)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.title = element_blank())

ggsave("Figure_2_Batch_Correction_PCA_All4.png", plot = p_before | p_after, width = 14, height = 5, dpi = 300)

save(harmonized_matrix, pheno, file = "Phase2_Final_Harmonized_Matrix.RData")
message("Saved Phase2_Final_Harmonized_Matrix.RData")