# =============================================================================
# Concordance Analysis of Shared Thermal-Response DEGs  |  Sharma et al.
# =============================================================================
# For genes significant in BOTH strains' thermal responses (T26 vs T20),
# classify by sign of log2FC in the two strains:
#   Concordant_UP / Concordant_DOWN  - same direction (conserved response)
#   Discordant_*                     - opposite direction (candidate divergent
#                                      regulatory adaptation)
#
# INPUT   Results/DEA_results.RData  (res_sh, res_otik)
# OUTPUT  Files 12-14 in Results/ ; Concordance_results.RData
#
# Reference temperature = T20, so positive log2FC = up at 26 C. A Concordant_DOWN
# gene is lower at 26 C (higher at 20 C) in both strains.
# =============================================================================

required_cran <- c("ggplot2", "dplyr", "tibble", "ggrepel")
for (pkg in required_cran) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
library(ggplot2); library(dplyr); library(tibble); library(ggrepel)

project_dir <- "."
results_dir <- file.path(project_dir, "Results")
load(file.path(results_dir, "DEA_results.RData"))  # res_sh, res_otik

sh_df   <- as.data.frame(res_sh)   |> rownames_to_column("Geneid")
otik_df <- as.data.frame(res_otik) |> rownames_to_column("Geneid")
sh_sig   <- sh_df   |> filter(padj < 0.05 & abs(log2FoldChange) > 1)
otik_sig <- otik_df |> filter(padj < 0.05 & abs(log2FoldChange) > 1)
shared_ids <- intersect(sh_sig$Geneid, otik_sig$Geneid)
cat("Shared thermal DEGs:", length(shared_ids), "\n")

shared <- sh_sig |> filter(Geneid %in% shared_ids) |>
  select(Geneid, log2FC_SH = log2FoldChange, padj_SH = padj) |>
  inner_join(otik_sig |> filter(Geneid %in% shared_ids) |>
               select(Geneid, log2FC_OTIK = log2FoldChange, padj_OTIK = padj), by = "Geneid") |>
  mutate(
    category = case_when(
      log2FC_SH > 0 & log2FC_OTIK > 0 ~ "Concordant_UP",
      log2FC_SH < 0 & log2FC_OTIK < 0 ~ "Concordant_DOWN",
      log2FC_SH > 0 & log2FC_OTIK < 0 ~ "Discordant_SHup_OTIKdown",
      log2FC_SH < 0 & log2FC_OTIK > 0 ~ "Discordant_SHdown_OTIKup"),
    interpretation = case_when(
      category == "Concordant_UP"            ~ "Both strains up at 26 C",
      category == "Concordant_DOWN"          ~ "Both strains down at 26 C (up at 20 C)",
      category == "Discordant_SHup_OTIKdown" ~ "SH up at 26 C, OTIK down at 26 C",
      category == "Discordant_SHdown_OTIKup" ~ "SH down at 26 C, OTIK up at 26 C"),
    magnitude = sqrt(log2FC_SH^2 + log2FC_OTIK^2),
    strain_difference = abs(log2FC_SH - log2FC_OTIK))

summary_table <- shared |> count(category, interpretation, name = "n_genes") |>
  mutate(percent = round(100 * n_genes / sum(n_genes), 1)) |> arrange(desc(n_genes))
print(summary_table)
write.csv(summary_table, file.path(results_dir, "12_Concordance_summary.csv"), row.names = FALSE)
write.csv(shared |> arrange(category, desc(magnitude)),
          file.path(results_dir, "12_Shared_thermal_DEGs_categorized.csv"), row.names = FALSE)
for (cat_name in unique(shared$category)) {
  write.csv(shared |> filter(category == cat_name) |> arrange(desc(magnitude)),
            file.path(results_dir, paste0("13_Genes_", cat_name, ".csv")), row.names = FALSE)
}

top_discordant <- shared |> filter(grepl("Discordant", category)) |>
  arrange(desc(strain_difference)) |> head(20)
cat_colors <- c("Concordant_UP" = "#d62728", "Concordant_DOWN" = "#1f77b4",
                "Discordant_SHup_OTIKdown" = "#ff7f0e", "Discordant_SHdown_OTIKup" = "#9467bd")
p_concord <- ggplot(shared, aes(x = log2FC_SH, y = log2FC_OTIK, color = category)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey30") +
  geom_point(alpha = 0.6, size = 1.6) +
  geom_text_repel(data = top_discordant, aes(label = Geneid), size = 3, color = "black",
                  max.overlaps = 30, box.padding = 0.4, segment.color = "grey50") +
  scale_color_manual(values = cat_colors) +
  labs(title = "Concordance of shared thermal DEGs",
       subtitle = paste0(nrow(shared), " genes DE in both SH and OTIK (T26 vs T20)"),
       x = "log2FC in SH (T26 vs T20)", y = "log2FC in OTIK (T26 vs T20)", color = "Category") +
  theme_bw(base_size = 12) + theme(legend.position = "right")
ggsave(file.path(results_dir, "14_Concordance_scatter.png"), p_concord, width = 11, height = 7, dpi = 300)

save(shared, summary_table, top_discordant, file = file.path(results_dir, "Concordance_results.RData"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "sessionInfo_05_thermal_concordance.txt"))
cat("\n=== Concordance analysis complete ===\n")
