# Five-Gene Transcriptomic Classifier for LUAD TME Phenotyping

Reproducible R pipeline for discovering and validating the locked five-gene humoral immune signature (**IGHM, CD38, IGHG3, P2RX5, FAM30A**) distinguishing Immune-Hot and Immune-Cold lung adenocarcinoma (LUAD) tumors.

> **Note:** Raw datasets, intermediate `.RData` files, generated tables, and figures are excluded from this repository. All outputs can be reproduced directly from public repositories by executing the scripts sequentially from the project root directory.

## Pipeline Scripts (`scripts/`)
1. `01_preprocessing_frma.R` — Downloads raw Affymetrix GPL570 `.CEL` archives (`GSE30219`, `GSE31210`, `GSE37745`, `GSE50081`) and performs parallel frozen Robust Multi-array Analysis (fRMA).
2. `02_annotation_and_combat.R` — Maps probes to HGNC symbols, filters histologically confirmed LUAD samples, and performs ComBat empirical Bayes batch harmonization.
3. `03_immune_deconvolution_clustering.R` — Exports non-negative LM22-intersected chunks for CIBERSORTx deconvolution, merges outputs, and performs Z-scored Ward.D2 hierarchical clustering (`k = 2`) to define Immune-Hot vs. Immune-Cold phenotypes.
4. `04_wgcna_network.R` — Constructs a signed WGCNA network (β = 9) on the top 5,000 most variable genes and extracts the phenotype-associated `MEmagenta` module.
5. `05_triml_feature_selection.R` — Runs LASSO, Random Forest, and Linear SVM-RFE on candidate module genes to isolate the locked five-gene consensus signature.
6. `06_nested_cv_and_loco.R` — Executes strictly nested 10-fold CV (GLM, RF, SVM-Radial, early-stopped ANN), end-to-end 5-fold nested discovery CV, and pre-ComBat leave-one-cohort-out (LOCO) validation.
7. `07_tcga_external_validation.R` — Downloads GDC TCGA-LUAD STAR-TPM (`n = 516`) and clinical survival data (`n = 495`) from UCSC Xena, runs local LM22 ν-SVR deconvolution, evaluates locked GEO classifiers zero-shot, tests orthogonal TME biomarkers (Ayers 18-gene GEP, TLS, CYT, checkpoints), and performs stage-adjusted Cox survival analysis.
8. `08_threshold_sensitivity.R` — Evaluates phenotype robustness across alternative clustering algorithms, single-lineage CD8+ vs. multi-lineage effector quantiles, 3-tier stratification (`k = 3`), and probability cutoffs.
9. `09_statistical_report.R` — Generates statistical summary tables (cohort χ² independence, 22-cell LM22 Wilcoxon/FDR tests, 5-gene multivariable GLM/VIF/Brier calibration, pairwise DeLong ROC tests, and package versions).
10. `10_generate_publication_figures.R` — Regenerates the high-resolution (600 DPI) manuscript figures.

## Public Data Sources
- **GEO Discovery Cohorts (Affymetrix GPL570):** `GSE30219`, `GSE31210`, `GSE37745`, `GSE50081`
- **External Validation Cohort (UCSC Xena):** GDC TCGA-LUAD STAR-TPM (`TCGA-LUAD.star_tpm.tsv.gz`, GENCODE v36) and clinical survival matrix (`LUAD_clinicalMatrix`)
- **Immune Signature Matrix:** Obtain `LM22.txt` from the [CIBERSORTx portal](https://cibersortx.stanford.edu/) (free for academic use) and place it in the project root directory prior to running scripts `03` and `07`.

## Software Environment
Verified in **R 4.6.1** (`GEOquery 2.80.0`, `affy 1.91.0`, `frma 1.65.0`, `hgu133plus2.db 3.13.0`, `sva 3.61.0`, `WGCNA 1.74`, `glmnet 5.0`, `randomForest 4.7.1.2`, `e1071 1.7.17`, `keras 2.16.1`, `pROC 1.19.1`, `survival 3.8.12`, `ggplot2 4.0.3`).
