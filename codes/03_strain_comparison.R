# =============================================================================
# Strain DEA at Fixed Temperatures (OTIK vs SH)  |  Sharma et al.
# =============================================================================
# Two strain contrasts, each holding temperature constant (OTIK vs SH at 26 C
# and at 20 C; reference = SH -> positive log2FC = higher in OTIK). Then compares
# the two contrasts (Venn; T26-vs-T20 log2FC scatter).
#
# INPUT   Results/DEA_results.RData  (from 02_deseq2_master.R)
# OUTPUT  Results/Strain_at_T26/, Results/Strain_at_T20/, files 15-19 in Results/,
#         and Strain_at_temperature_results.RData
# =============================================================================

required_cran <- c("ggplot2", "dplyr", "tibble", "ggvenn", "ggrepel")
required_bioc <- c("DESeq2", "EnhancedVolcano")
for (pkg in required_cran) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in required_bioc) if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE, ask = FALSE)
library(DESeq2); library(ggplot2); library(dplyr); library(tibble)
library(ggvenn); library(ggrepel); library(EnhancedVolcano)

project_dir <- "."
results_dir <- file.path(project_dir, "Results")
dir_T26 <- file.path(results_dir, "Strain_at_T26"); dir_T20 <- file.path(results_dir, "Strain_at_T20")
if (!dir.exists(dir_T26)) dir.create(dir_T26, recursive = TRUE)
if (!dir.exists(dir_T20)) dir.create(dir_T20, recursive = TRUE)

load(file.path(results_dir, "DEA_results.RData"))
counts_mat <- counts(dds); meta <- as.data.frame(colData(dds))

run_strain_at_temp <- function(temp_value, out_dir, label) {
  keep_samples <- rownames(meta)[meta$temperature == temp_value]
  dds_sub <- DESeqDataSetFromMatrix(countData = counts_mat[, keep_samples],
                                    colData = meta[keep_samples, ], design = ~ strain)
  dds_sub <- dds_sub[rowSums(counts(dds_sub) >= 10) >= 3, ]
  dds_sub <- DESeq(dds_sub)
  res <- results(dds_sub, contrast = c("strain", "OTIK", "SH"), alpha = 0.05)
  res <- res[order(res$padj), ]
  write.csv(as.data.frame(res) |> rownames_to_column("Geneid"),
            file.path(out_dir, paste0("DEGs_", label, "_OTIKvsSH_full.csv")), row.names = FALSE)
  sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
  write.csv(as.data.frame(sig) |> rownames_to_column("Geneid"),
            file.path(out_dir, paste0("DEGs_", label, "_OTIKvsSH_significant.csv")), row.names = FALSE)
  cat(label, "significant DEGs:", nrow(sig), "\n")
  p_volcano <- EnhancedVolcano(res, lab = rownames(res), x = "log2FoldChange", y = "padj",
                               pCutoff = 0.05, FCcutoff = 1,
                               title = paste0("OTIK vs SH at ", temp_value),
                               subtitle = "Positive log2FC = higher in OTIK",
                               pointSize = 1.8, labSize = 3, colAlpha = 0.7, legendPosition = "right")
  ggsave(file.path(out_dir, paste0("Volcano_", label, "_OTIKvsSH.png")), p_volcano, width = 9, height = 7, dpi = 300)
  list(dds = dds_sub, res = res, sig = sig)
}

out_T26 <- run_strain_at_temp("T26", dir_T26, "T26")
out_T20 <- run_strain_at_temp("T20", dir_T20, "T20")
res_T26 <- out_T26$res; res_T20 <- out_T20$res
sig_T26 <- out_T26$sig; sig_T20 <- out_T20$sig

venn_list <- list(`At 26 C` = rownames(sig_T26), `At 20 C` = rownames(sig_T20))
p_venn <- ggvenn(venn_list, fill_color = c("#d62728", "#1f77b4"),
                 stroke_size = 0.5, set_name_size = 6, text_size = 5) +
  ggtitle("Strain DEGs (OTIK vs SH) at each temperature")
ggsave(file.path(results_dir, "15_Venn_Strain_at_T26_vs_T20.png"), p_venn, width = 6, height = 5, dpi = 300)

shared_both <- intersect(rownames(sig_T26), rownames(sig_T20))
T26_only <- setdiff(rownames(sig_T26), rownames(sig_T20))
T20_only <- setdiff(rownames(sig_T20), rownames(sig_T26))
write.csv(data.frame(Geneid = shared_both), file.path(results_dir, "16_StrainDEGs_constitutive_bothTemps.csv"), row.names = FALSE)
write.csv(data.frame(Geneid = T26_only), file.path(results_dir, "16_StrainDEGs_T26only_emergesUnderHeat.csv"), row.names = FALSE)
write.csv(data.frame(Geneid = T20_only), file.path(results_dir, "16_StrainDEGs_T20only_emergesUnderCool.csv"), row.names = FALSE)

T26_df <- as.data.frame(res_T26) |> rownames_to_column("Geneid") |> select(Geneid, log2FC_T26 = log2FoldChange, padj_T26 = padj)
T20_df <- as.data.frame(res_T20) |> rownames_to_column("Geneid") |> select(Geneid, log2FC_T20 = log2FoldChange, padj_T20 = padj)
combined <- inner_join(T26_df, T20_df, by = "Geneid") |>
  filter(!is.na(log2FC_T26) & !is.na(log2FC_T20)) |>
  mutate(sig_T26 = !is.na(padj_T26) & padj_T26 < 0.05 & abs(log2FC_T26) > 1,
         sig_T20 = !is.na(padj_T20) & padj_T20 < 0.05 & abs(log2FC_T20) > 1,
         category = case_when(
           sig_T26 &  sig_T20 ~ "Sig in both (constitutive)",
           sig_T26 & !sig_T20 ~ "Sig at T26 only (heat-specific)",
           !sig_T26 & sig_T20 ~ "Sig at T20 only (cool-specific)",
           TRUE               ~ "Not significant"),
         strain_diff_change = abs(log2FC_T26 - log2FC_T20))
write.csv(combined, file.path(results_dir, "17_Strain_DEGs_combined_T26_T20.csv"), row.names = FALSE)

top_changers <- combined |> filter(category != "Not significant") |> arrange(desc(strain_diff_change)) |> head(20)
cat_colors <- c("Sig in both (constitutive)" = "#7f7f7f", "Sig at T26 only (heat-specific)" = "#d62728",
                "Sig at T20 only (cool-specific)" = "#1f77b4", "Not significant" = "#e0e0e0")
p_scatter <- ggplot(combined, aes(x = log2FC_T20, y = log2FC_T26, color = category)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey30") +
  geom_point(alpha = 0.6, size = 1.4) +
  geom_text_repel(data = top_changers, aes(label = Geneid), size = 3, color = "black",
                  max.overlaps = 30, box.padding = 0.4, segment.color = "grey50") +
  scale_color_manual(values = cat_colors) +
  labs(title = "Strain difference (OTIK vs SH) across temperatures",
       x = "log2FC (OTIK vs SH) at 20 C", y = "log2FC (OTIK vs SH) at 26 C", color = "Category") +
  theme_bw(base_size = 12) + theme(legend.position = "right")
ggsave(file.path(results_dir, "18_Strain_DEGs_T26vsT20_scatter.png"), p_scatter, width = 11, height = 7, dpi = 300)

summary_df <- data.frame(
  Contrast = c("OTIK vs SH at 26 C", "OTIK vs SH at 20 C"),
  Total_tested = c(nrow(res_T26), nrow(res_T20)),
  Sig_padj_FC1 = c(nrow(sig_T26), nrow(sig_T20)),
  Up_in_OTIK = c(sum(sig_T26$log2FoldChange > 0, na.rm = TRUE), sum(sig_T20$log2FoldChange > 0, na.rm = TRUE)),
  Up_in_SH = c(sum(sig_T26$log2FoldChange < 0, na.rm = TRUE), sum(sig_T20$log2FoldChange < 0, na.rm = TRUE)))
write.csv(summary_df, file.path(results_dir, "19_Strain_at_temperature_summary.csv"), row.names = FALSE)
print(summary_df)
cat("Constitutive:", length(shared_both), "| T26-only:", length(T26_only), "| T20-only:", length(T20_only), "\n")

save(out_T26, out_T20, res_T26, res_T20, sig_T26, sig_T20, combined,
     file = file.path(results_dir, "Strain_at_temperature_results.RData"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "sessionInfo_03_strain_comparison.txt"))
cat("\n=== Analysis complete ===\n")
