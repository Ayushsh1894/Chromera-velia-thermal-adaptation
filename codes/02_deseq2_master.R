# =============================================================================
# Master DESeq2 Analysis  —  Chromera velia thermal-adaptation RNA-seq
# Sharma et al.  |  Two strains (OTIK, SH) x two temperatures (T20, T26)
# =============================================================================
#
# PURPOSE
#   Builds the core DESeq2 objects and differential-expression contrasts that
#   every downstream script in this repository depends on, and saves them to
#   "DEA_results.RData".
#
#   Objects saved:
#     dds         - main DESeqDataSet, design = ~ strain + temperature
#     dds_sh      - SH-only subset,   design = ~ temperature (thermal contrast)
#     dds_otik    - OTIK-only subset, design = ~ temperature (thermal contrast)
#     vsd         - variance-stabilised transform of dds (for PCA / heatmaps)
#     res_strain  - overall strain effect (OTIK vs SH), controlling for temp
#     res_sh      - SH thermal response   (T26 vs T20)
#     res_otik    - OTIK thermal response (T26 vs T20)
#     sig_strain / sig_sh / sig_otik - significant subsets (padj<0.05, |LFC|>1)
#     metadata    - sample table
#
# INPUTS (in the working directory)
#     Count_matrix.xlsx            - gene x sample raw counts (30492 x 12);
#                                    first column "Geneid" (featureCounts output)
#     experimental_design_temp.csv - sample metadata: sample, strain, temperature
#
# CONVENTIONS
#   Reference levels: strain = SH, temperature = T20.
#     res_strain          : positive log2FC = higher in OTIK
#     res_sh / res_otik   : positive log2FC = up at 26 C (heat-induced)
#   Significance threshold throughout: padj < 0.05 AND |log2FoldChange| > 1.
#   Pre-filter: keep genes with >= 10 counts in >= 3 samples (smallest group n).
#
# NORMALISATION NOTE
#   Differential expression uses DESeq2 median-of-ratios normalisation.
#   Presence/absence (core transcriptome, script 01) instead uses edgeR CPM;
#   the two normalisations answer different questions and are not interchangeable.
#
# WITHIN-STRAIN CONTRASTS
#   Thermal responses are computed on per-strain subsets (each with ~ temperature)
#   rather than from an interaction term on the full model. This avoids
#   model-fitting issues arising from the very large strain effect.
# =============================================================================


# ---- 1. Packages -------------------------------------------------------------

required_cran <- c("readxl", "ggplot2", "dplyr", "tibble", "pheatmap",
                   "RColorBrewer")
required_bioc <- c("DESeq2", "EnhancedVolcano")

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

library(DESeq2)
library(readxl)
library(ggplot2)
library(dplyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)


# ---- 2. Paths ----------------------------------------------------------------
# EDIT THIS to the folder containing the two input files. Kept relative so the
# repository is portable; set to "." if inputs are in the working directory.

project_dir <- "."
setwd(project_dir)

results_dir <- file.path(project_dir, "Results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

cat("Working directory:", getwd(), "\n")
cat("Results ->", results_dir, "\n\n")


# ---- 3. Load Count Matrix ----------------------------------------------------

counts <- read_excel("Count_matrix.xlsx") |> as.data.frame()
rownames(counts) <- counts$Geneid
counts$Geneid <- NULL
counts <- as.matrix(counts)
mode(counts) <- "integer"          # DESeq2 requires integer counts

cat("Count matrix:", nrow(counts), "genes x", ncol(counts), "samples\n")


# ---- 4. Load Metadata --------------------------------------------------------

metadata <- read.csv("experimental_design_temp.csv", stringsAsFactors = FALSE)
rownames(metadata) <- metadata$sample

# Match count column order to metadata row order
counts <- counts[, rownames(metadata)]
stopifnot(all(colnames(counts) == rownames(metadata)))

# Reference levels: SH and T20
metadata$strain      <- factor(metadata$strain,      levels = c("SH", "OTIK"))
metadata$temperature <- factor(metadata$temperature, levels = c("T20", "T26"))

cat("Strain levels:     ", paste(levels(metadata$strain), collapse = ", "),
    "(first = reference)\n")
cat("Temperature levels:", paste(levels(metadata$temperature), collapse = ", "),
    "(first = reference)\n\n")


# ---- 5. Full DESeq2 Model:  ~ strain + temperature --------------------------

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = metadata,
                              design    = ~ strain + temperature)

# Pre-filter: keep genes with >= 10 counts in at least 3 samples (smallest group)
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
cat("After low-count filtering:", nrow(dds), "genes retained\n")

dds <- DESeq(dds)


# ---- 6. PCA on vst-transformed counts ---------------------------------------

vsd <- vst(dds, blind = TRUE)

pca_data    <- plotPCA(vsd, intgroup = c("strain", "temperature"),
                       returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = strain, shape = temperature)) +
  geom_point(size = 5, alpha = 0.85) +
  scale_color_manual(values = c("SH" = "#1f77b4", "OTIK" = "#d62728")) +
  xlab(paste0("PC1 (", percent_var[1], "%)")) +
  ylab(paste0("PC2 (", percent_var[2], "%)")) +
  ggtitle("PCA: Strain \u00d7 Temperature") +
  theme_bw(base_size = 13) +
  theme(legend.position = "right")

ggsave(file.path(results_dir, "01_PCA_strain_temperature.png"),
       p_pca, width = 7, height = 5, dpi = 300)


# ---- 7. Contrast 1: Strain Main Effect (OTIK vs SH) -------------------------
# Controls for temperature via the full model design.

res_strain <- results(dds, contrast = c("strain", "OTIK", "SH"), alpha = 0.05)
res_strain <- res_strain[order(res_strain$padj), ]

cat("\n--- Strain effect (OTIK vs SH) ---\n")
summary(res_strain)

write.csv(as.data.frame(res_strain) |> rownames_to_column("Geneid"),
          file.path(results_dir, "02_DEGs_Strain_OTIKvsSH_full.csv"),
          row.names = FALSE)

sig_strain <- subset(res_strain, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(sig_strain) |> rownames_to_column("Geneid"),
          file.path(results_dir, "02_DEGs_Strain_OTIKvsSH_significant.csv"),
          row.names = FALSE)


# ---- 8. Within-Strain Temperature Contrasts ---------------------------------
# Subset each strain and run a separate DESeq2 model with ~ temperature.
# Positive log2FC = up at 26 C.

run_within_strain <- function(strain_name) {
  samples <- rownames(metadata)[metadata$strain == strain_name]
  dds_sub <- DESeqDataSetFromMatrix(countData = counts[, samples],
                                    colData   = metadata[samples, ],
                                    design    = ~ temperature)
  dds_sub <- dds_sub[rowSums(counts(dds_sub) >= 10) >= 3, ]
  dds_sub <- DESeq(dds_sub)
  res_sub <- results(dds_sub, contrast = c("temperature", "T26", "T20"),
                     alpha = 0.05)
  res_sub <- res_sub[order(res_sub$padj), ]
  list(dds = dds_sub, res = res_sub)
}

# SH thermal response
sh_out <- run_within_strain("SH")
dds_sh <- sh_out$dds
res_sh <- sh_out$res
cat("\n--- SH: T26 vs T20 ---\n"); summary(res_sh)

write.csv(as.data.frame(res_sh) |> rownames_to_column("Geneid"),
          file.path(results_dir, "03_DEGs_SH_T26vsT20_full.csv"),
          row.names = FALSE)
sig_sh <- subset(res_sh, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(sig_sh) |> rownames_to_column("Geneid"),
          file.path(results_dir, "03_DEGs_SH_T26vsT20_significant.csv"),
          row.names = FALSE)

# OTIK thermal response
otik_out <- run_within_strain("OTIK")
dds_otik <- otik_out$dds
res_otik <- otik_out$res
cat("\n--- OTIK: T26 vs T20 ---\n"); summary(res_otik)

write.csv(as.data.frame(res_otik) |> rownames_to_column("Geneid"),
          file.path(results_dir, "04_DEGs_OTIK_T26vsT20_full.csv"),
          row.names = FALSE)
sig_otik <- subset(res_otik, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(sig_otik) |> rownames_to_column("Geneid"),
          file.path(results_dir, "04_DEGs_OTIK_T26vsT20_significant.csv"),
          row.names = FALSE)


# ---- 9. Top-50 Variable-Gene Heatmap (vst, row-centred) ---------------------

top50 <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat   <- assay(vsd)[top50, ]
mat   <- mat - rowMeans(mat)
ann_col <- as.data.frame(colData(vsd)[, c("strain", "temperature")])

pheatmap(mat,
         annotation_col = ann_col,
         show_rownames  = FALSE,
         color          = colorRampPalette(rev(brewer.pal(11, "RdBu")))(255),
         main           = "Top 50 most variable genes (vst, row-centred)",
         filename       = file.path(results_dir, "10_Heatmap_top50_variable_genes.png"),
         width = 8, height = 8)


# ---- 10. Summary Table of DEG Counts ----------------------------------------

summary_stats <- data.frame(
  Contrast = c("Strain (OTIK vs SH, temp-controlled)",
               "SH thermal (T26 vs T20)",
               "OTIK thermal (T26 vs T20)"),
  Total_tested = c(sum(!is.na(res_strain$padj)),
                   sum(!is.na(res_sh$padj)),
                   sum(!is.na(res_otik$padj))),
  Sig_padj_LFC = c(nrow(sig_strain), nrow(sig_sh), nrow(sig_otik))
)
write.csv(summary_stats,
          file.path(results_dir, "11_DEG_summary.csv"), row.names = FALSE)
cat("\n=== Summary ===\n"); print(summary_stats, row.names = FALSE); cat("\n")


# ---- 11. Save All Objects (loaded by every downstream script) ---------------

save(dds, dds_sh, dds_otik, vsd,
     res_strain, res_sh, res_otik,
     sig_strain, sig_sh, sig_otik,
     metadata,
     file = file.path(results_dir, "DEA_results.RData"))

cat("Saved:", file.path(results_dir, "DEA_results.RData"), "\n")


# ---- 12. Session Info -------------------------------------------------------

writeLines(capture.output(sessionInfo()),
           file.path(results_dir, "sessionInfo_02_deseq2_master.txt"))

cat("\n=== Master DESeq2 analysis complete ===\n")
