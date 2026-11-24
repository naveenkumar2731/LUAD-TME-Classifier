# ==============================================================================
# 03_immune_deconvolution_clustering.R
# Part A: Export LM22-intersected non-negative expression chunks for CIBERSORTx.
# Part B: Merge CIBERSORTx LM22 outputs, Z-score normalize cell fractions,
#         perform Ward.D2 hierarchical clustering (k = 2), and assign
#         Immune-Hot vs. Immune-Cold phenotypes by CD8+ T-cell enrichment.
# ==============================================================================

suppressPackageStartupMessages({
  library(pheatmap)
  library(stats)
})

load("Phase2_Final_Harmonized_Matrix.RData")

# ------------------------------------------------------------------------------
# PART A: Export non-negative LM22 subsets for CIBERSORTx deconvolution
# ------------------------------------------------------------------------------
ciber_files <- c("CIBERSORTx_Job4_Results.csv",
                 "CIBERSORTx_Job5_Results.csv",
                 "CIBERSORTx_Job6_Results.csv")

if (!all(file.exists(ciber_files))) {
  lm22 <- read.delim("LM22.txt", row.names = 1, stringsAsFactors = FALSE)
  mat_mini <- harmonized_matrix[intersect(rownames(lm22), rownames(harmonized_matrix)), , drop = FALSE]
  mat_mini[mat_mini < 0] <- 0

  write.table(mat_mini[, 1:100],                  "Phase3_Chunk1.txt", sep = "\t", quote = FALSE, col.names = NA)
  write.table(mat_mini[, 101:200],                "Phase3_Chunk2.txt", sep = "\t", quote = FALSE, col.names = NA)
  write.table(mat_mini[, 201:ncol(mat_mini)],     "Phase3_Chunk3.txt", sep = "\t", quote = FALSE, col.names = NA)

  stop("Exported Phase3_Chunk1-3.txt. Run CIBERSORTx (LM22) and place result CSVs in project root, then re-run.")
}

# ------------------------------------------------------------------------------
# PART B: Merge CIBERSORTx outputs and stratify Immune-Hot vs. Immune-Cold
# ------------------------------------------------------------------------------
results_combined <- do.call(rbind, lapply(ciber_files, function(f) {
  read.csv(f, row.names = 1, check.names = FALSE)
}))

# Retain the 22 immune cell fractions (excluding P-value, Correlation, RMSE)
immune_fractions <- results_combined[, 1:22]

common_samples   <- intersect(rownames(immune_fractions), pheno$Sample)
immune_fractions <- immune_fractions[common_samples, , drop = FALSE]
pheno_immune     <- pheno[match(common_samples, pheno$Sample), ]
harmonized_matrix <- harmonized_matrix[, common_samples, drop = FALSE]

# Z-score each lineage across patients prior to Euclidean Ward.D2 clustering
scaled_fractions <- scale(immune_fractions)
dist_mat <- dist(scaled_fractions, method = "euclidean")
hc       <- hclust(dist_mat, method = "ward.D2")
clusters <- cutree(hc, k = 2)

# Designate cluster with higher mean standardized CD8+ T-cell infiltration as Immune-Hot
cd8_means      <- aggregate(scaled_fractions[, "T cells CD8"], list(Cluster = clusters), mean)
hot_cluster_id <- cd8_means$Cluster[which.max(cd8_means$x)]

pheno_immune$Immune_Phenotype <- factor(
  ifelse(clusters == hot_cluster_id, "Immune-Hot", "Immune-Cold"),
  levels = c("Immune-Cold", "Immune-Hot")
)

message("Immune phenotype distribution:")
print(table(pheno_immune$Immune_Phenotype))

# Generate Z-scored LM22 immune landscape heatmap
annotation_col <- data.frame(Phenotype = pheno_immune$Immune_Phenotype)
rownames(annotation_col) <- pheno_immune$Sample
pheno_colors <- list(Phenotype = c("Immune-Hot" = "#E64B35FF", "Immune-Cold" = "#3C5488FF"))

png("Figure_3_Immune_Landscape_Heatmap_Corrected.png", width = 12, height = 8, units = "in", res = 300)
pheatmap(t(immune_fractions),
         scale = "row",
         annotation_col = annotation_col,
         annotation_colors = pheno_colors,
         show_colnames = FALSE,
         clustering_method = "ward.D2",
         color = colorRampPalette(c("#3C5488FF", "white", "#E64B35FF"))(100),
         main = "LM22 Immune Infiltration Profiles (Z-Scored)")
dev.off()

save(harmonized_matrix, pheno_immune, immune_fractions, file = "Phase3_Immune_Phenotypes.RData")
message("Saved Phase3_Immune_Phenotypes.RData")