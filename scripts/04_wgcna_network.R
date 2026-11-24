# ==============================================================================
# 04_wgcna_network.R
# Construct a signed Weighted Gene Co-expression Network (WGCNA) on the top
# 5,000 most variable genes, identify phenotype-associated modules, and export
# candidate genes from the module most strongly correlated with immune state.
# ==============================================================================

suppressPackageStartupMessages({
  library(WGCNA)
})

allowWGCNAThreads()

# 1. Load Phase 3 expression matrix and immune phenotypes
load("Phase3_Immune_Phenotypes.RData")

datExpr <- t(harmonized_matrix)

# 2. Retain top 5,000 most variable genes
gene_variances    <- apply(datExpr, 2, var, na.rm = TRUE)
top_variant_genes <- order(gene_variances, decreasing = TRUE)[1:5000]
datExpr_filtered  <- datExpr[, top_variant_genes]

# Encode binary trait (Immune-Cold = 1, Immune-Hot = 0)
datTraits <- data.frame(
  Immune_Cold = ifelse(pheno_immune$Immune_Phenotype == "Immune-Cold", 1, 0),
  row.names   = pheno_immune$Sample
)

# 3. Evaluate soft-thresholding powers for scale-free topology (signed network)
powers <- c(1:10, seq(from = 12, to = 20, by = 2))
sft    <- pickSoftThreshold(datExpr_filtered, powerVector = powers,
                            networkType = "signed", verbose = 0)

save(datExpr_filtered, datTraits, sft, file = "Phase4_WGCNA_Input.RData")

# 4. Construct signed adjacency and Topological Overlap Matrix (TOM) at beta = 9
softPower <- 9
adjacency <- adjacency(datExpr_filtered, power = softPower, type = "signed")
TOM       <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM   <- 1 - TOM
rm(adjacency)
gc()

# 5. Hierarchical clustering, dynamic tree cut, and module merging (r >= 0.75)
geneTree      <- hclust(as.dist(dissTOM), method = "average")
dynamicMods   <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                               deepSplit = 2, pamRespectsDendro = FALSE,
                               minClusterSize = 30)
dynamicColors <- labels2colors(dynamicMods)

merge        <- mergeCloseModules(datExpr_filtered, dynamicColors,
                                  cutHeight = 0.25, verbose = 0)
mergedColors <- merge$colors
mergedMEs    <- merge$newMEs

# 6. Exclude unassigned grey module and compute module-trait correlations
nSamples                <- nrow(datExpr_filtered)
mergedMEs_clean         <- mergedMEs[, names(mergedMEs) != "MEgrey", drop = FALSE]
moduleTraitCor_clean    <- cor(mergedMEs_clean, datTraits$Immune_Cold, use = "p")
moduleTraitPvalue_clean <- corPvalueStudent(moduleTraitCor_clean, nSamples)

# 7. Extract candidate genes from module with strongest absolute trait correlation
top_module_idx  <- which.max(abs(moduleTraitCor_clean[, 1]))
top_module_name <- substring(rownames(moduleTraitCor_clean)[top_module_idx], 3)
top_cor_val     <- moduleTraitCor_clean[top_module_idx, 1]
top_pval        <- moduleTraitPvalue_clean[top_module_idx, 1]

target_genes <- colnames(datExpr_filtered)[mergedColors == top_module_name]
message(sprintf("Top module: '%s' (r = %.3f, p = %.2e, n = %d candidate genes)",
                top_module_name, top_cor_val, top_pval, length(target_genes)))

save(mergedMEs_clean, mergedColors, target_genes, top_module_name,
     moduleTraitCor_clean, moduleTraitPvalue_clean,
     file = "Phase4_WGCNA_Modules_Corrected.RData")
write.csv(data.frame(Gene = target_genes),
          file = "Phase4_Target_Module_Genes.csv", row.names = FALSE)
message("Saved Phase4_WGCNA_Modules_Corrected.RData")