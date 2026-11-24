# ==============================================================================
# 01_preprocessing_frma.R
# Download raw Affymetrix GPL570 .CEL archives for four GEO LUAD cohorts,
# execute parallel frozen Robust Multi-array Analysis (fRMA), and save the
# unified probe-level expression matrix.
# ==============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(affy)
  library(frma)
  library(hgu133plus2frmavecs)
  library(parallel)
  library(doParallel)
  library(foreach)
})

options(timeout = 3600)

datasets <- c("GSE31210", "GSE50081", "GSE37745", "GSE30219")

# 1. Download and extract raw .CEL archives
for (geo_id in datasets) {
  tar_file <- file.path(geo_id, paste0(geo_id, "_RAW.tar"))
  cel_dir  <- file.path(geo_id, "CEL")

  if (!file.exists(tar_file)) {
    message("Downloading raw .CEL archive for ", geo_id, "...")
    getGEOSuppFiles(geo_id, makeDirectory = TRUE)
  }

  existing_cels <- list.files(cel_dir, pattern = "\\.CEL(\\.gz)?$", ignore.case = TRUE)
  if (!dir.exists(cel_dir) || length(existing_cels) == 0) {
    message("Extracting .CEL files for ", geo_id, "...")
    dir.create(cel_dir, showWarnings = FALSE, recursive = TRUE)
    untar(tar_file, exdir = cel_dir, tar = "internal")
  }
}

# 2. Locate all extracted .CEL files
all_cels <- unlist(lapply(datasets, function(d) {
  list.files(file.path(d, "CEL"), pattern = "\\.CEL(\\.gz)?$",
             full.names = TRUE, ignore.case = TRUE)
}))
stopifnot("No .CEL files found in cohort directories." = length(all_cels) > 0)

# 3. Parallel chunked fRMA normalization
n_cores <- min(10L, max(1L, detectCores() - 1L))
cl <- makeCluster(n_cores)
registerDoParallel(cl)

chunk_size <- 10L
cel_chunks <- split(all_cels, ceiling(seq_along(all_cels) / chunk_size))

message("Running fRMA normalization across ", length(all_cels), " arrays (", n_cores, " cores)...")
results <- foreach(chunk = cel_chunks,
                   .packages = c("affy", "frma", "hgu133plus2frmavecs")) %dopar% {
  raw_batch  <- ReadAffy(filenames = chunk)
  norm_batch <- frma(raw_batch)
  out <- list(
    raw  = log2(exprs(raw_batch)),
    norm = exprs(norm_batch)
  )
  rm(raw_batch, norm_batch)
  gc()
  out
}
stopCluster(cl)

combined_raw  <- do.call(cbind, lapply(results, `[[`, "raw"))
combined_frma <- do.call(cbind, lapply(results, `[[`, "norm"))
rm(results)
gc()

# 4. Generate QC boxplot before vs. after fRMA (subsampled for plotting efficiency)
set.seed(42)
idx_raw  <- sample(seq_len(nrow(combined_raw)),  min(10000L, nrow(combined_raw)))
idx_frma <- sample(seq_len(nrow(combined_frma)), min(10000L, nrow(combined_frma)))

cohort_colors <- ifelse(
  grepl("GSE31210", all_cels), "#E64B35FF",
  ifelse(grepl("GSE50081", all_cels), "#4DBBD5FF",
         ifelse(grepl("GSE37745", all_cels), "#00A087FF", "#3C5488FF"))
)

png("Figure_QC_fRMA_Combined.png", width = 3500, height = 1500, res = 300)
par(mfrow = c(1, 2), mar = c(4, 5, 4, 2), bg = "white")
boxplot(combined_raw[idx_raw, ],
        main = "All Cohorts - Before fRMA (Raw)",
        col = cohort_colors, outline = FALSE, xaxt = "n",
        cex.main = 1.5, cex.lab = 1.2, ylab = "Log2 Probe Intensity")
boxplot(combined_frma[idx_frma, ],
        main = "All Cohorts - After fRMA (Normalized)",
        col = cohort_colors, outline = FALSE, xaxt = "n",
        cex.main = 1.5, cex.lab = 1.2, ylab = "fRMA Expression Level")
legend("bottomright", legend = datasets,
       fill = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF"),
       cex = 1.2, bg = "white", box.col = "black")
dev.off()

# 5. Save normalized expression matrix
save(combined_frma, file = "fRMA_Normalized_Data.RData")
message("Saved fRMA_Normalized_Data.RData")