# =============================================================================
# Direction-Resolved Categorization of Strain DEGs  |  Sharma et al.
# =============================================================================
# Splits significant strain DEGs into 8 categories by (a) which temperatures they
# are significant at and (b) which strain is higher. Produces the strain-difference
# gene lists submitted to CryptoDB (Strain_OTI_heat_up, etc.).
#
# INPUT   Results/Strain_at_temperature_results.RData  (from script 03)
# OUTPUT  Files 20-24 in Results/ ; Directional_categories_results.RData
# =============================================================================

required_cran <- c("ggplot2", "dplyr", "tibble", "UpSetR")
for (pkg in required_cran) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
library(ggplot2); library(dplyr); library(tibble); library(UpSetR)

project_dir <- "."
results_dir <- file.path(project_dir, "Results")
load(file.path(results_dir, "Strain_at_temperature_results.RData"))  # -> combined

final <- combined |>
  mutate(final_category = case_when(
    sig_T20  & sig_T26  & log2FC_T20 > 0 & log2FC_T26 > 0 ~ "Constitutive_OTIK_high",
    sig_T20  & sig_T26  & log2FC_T20 < 0 & log2FC_T26 < 0 ~ "Constitutive_SH_high",
    sig_T20  & sig_T26  & log2FC_T20 > 0 & log2FC_T26 < 0 ~ "Direction_flip_OTIK20_SH26",
    sig_T20  & sig_T26  & log2FC_T20 < 0 & log2FC_T26 > 0 ~ "Direction_flip_SH20_OTIK26",
    !sig_T20 & sig_T26  & log2FC_T26 > 0                  ~ "Heat_specific_OTIK_up",
    !sig_T20 & sig_T26  & log2FC_T26 < 0                  ~ "Heat_specific_SH_up",
    sig_T20  & !sig_T26 & log2FC_T20 > 0                  ~ "Cool_specific_OTIK_up",
    sig_T20  & !sig_T26 & log2FC_T20 < 0                  ~ "Cool_specific_SH_up",
    TRUE                                                  ~ "Not_significant"))

category_counts <- final |> filter(final_category != "Not_significant") |>
  count(final_category, sort = TRUE) |> mutate(percent = round(100 * n / sum(n), 1))
print(category_counts)
write.csv(category_counts, file.path(results_dir, "20_Directional_categories_counts.csv"), row.names = FALSE)

for (cat_name in setdiff(unique(final$final_category), "Not_significant")) {
  cat_genes <- final |> filter(final_category == cat_name) |>
    select(Geneid, log2FC_T20, padj_T20, log2FC_T26, padj_T26) |>
    arrange(desc(pmax(abs(log2FC_T20), abs(log2FC_T26))))
  write.csv(cat_genes, file.path(results_dir, paste0("21_Genes_", cat_name, ".csv")), row.names = FALSE)
}

cat_colors <- c(
  "Constitutive_OTIK_high" = "#c0392b", "Constitutive_SH_high" = "#2874a6",
  "Heat_specific_OTIK_up" = "#e74c3c", "Heat_specific_SH_up" = "#5dade2",
  "Cool_specific_OTIK_up" = "#f1948a", "Cool_specific_SH_up" = "#aed6f1",
  "Direction_flip_OTIK20_SH26" = "#7f8c8d", "Direction_flip_SH20_OTIK26" = "#bdc3c7")

p_bar <- ggplot(category_counts, aes(x = reorder(final_category, n), y = n, fill = final_category)) +
  geom_col() + geom_text(aes(label = paste0(n, " (", percent, "%)")), hjust = -0.1, size = 4) +
  coord_flip() + scale_fill_manual(values = cat_colors) +
  labs(title = "Direction-resolved categorization of strain DEGs",
       subtitle = paste0("Total: ", sum(category_counts$n), " significant genes"),
       x = NULL, y = "Number of genes") +
  theme_bw(base_size = 12) + theme(legend.position = "none", plot.title.position = "plot") +
  expand_limits(y = max(category_counts$n) * 1.25)
ggsave(file.path(results_dir, "22_DirectionalCategories_barplot.png"), p_bar, width = 11, height = 6, dpi = 300)

gene_sets <- list(
  "OTIK higher at 20 C" = final$Geneid[final$sig_T20 & final$log2FC_T20 > 0],
  "SH higher at 20 C"   = final$Geneid[final$sig_T20 & final$log2FC_T20 < 0],
  "OTIK higher at 26 C" = final$Geneid[final$sig_T26 & final$log2FC_T26 > 0],
  "SH higher at 26 C"   = final$Geneid[final$sig_T26 & final$log2FC_T26 < 0])
png(file.path(results_dir, "23_UpSet_4_directional_sets.png"), width = 12, height = 7, units = "in", res = 300)
upset(fromList(gene_sets), order.by = "freq", nsets = 4, nintersects = NA,
      sets.bar.color = c("#c0392b", "#e74c3c", "#2874a6", "#5dade2"),
      main.bar.color = "#34495e", matrix.color = "#34495e", shade.color = "#ecf0f1",
      mainbar.y.label = "Genes in intersection", sets.x.label = "Genes per set",
      text.scale = c(1.5, 1.3, 1.3, 1.1, 1.4, 1.2), point.size = 3, line.size = 1)
dev.off()

sig_only <- final |> filter(final_category != "Not_significant")
q_ur <- sum(sig_only$log2FC_T20 > 0 & sig_only$log2FC_T26 > 0)
q_ul <- sum(sig_only$log2FC_T20 < 0 & sig_only$log2FC_T26 > 0)
q_ll <- sum(sig_only$log2FC_T20 < 0 & sig_only$log2FC_T26 < 0)
q_lr <- sum(sig_only$log2FC_T20 > 0 & sig_only$log2FC_T26 < 0)
max_x <- max(abs(sig_only$log2FC_T20)) * 1.05; max_y <- max(abs(sig_only$log2FC_T26)) * 1.05
p_quad <- ggplot(sig_only, aes(x = log2FC_T20, y = log2FC_T26, color = final_category)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey30", linewidth = 0.4) +
  geom_point(alpha = 0.55, size = 1.1) + scale_color_manual(values = cat_colors, name = "Category") +
  annotate("label", x = max_x*0.7, y = max_y*0.95, label = paste0("OTIK higher both\nn = ", q_ur), fill = "white", size = 3.6, alpha = 0.85) +
  annotate("label", x = -max_x*0.7, y = max_y*0.95, label = paste0("Flip: SH@20, OTIK@26\nn = ", q_ul), fill = "white", size = 3.6, alpha = 0.85) +
  annotate("label", x = -max_x*0.7, y = -max_y*0.95, label = paste0("SH higher both\nn = ", q_ll), fill = "white", size = 3.6, alpha = 0.85) +
  annotate("label", x = max_x*0.7, y = -max_y*0.95, label = paste0("Flip: OTIK@20, SH@26\nn = ", q_lr), fill = "white", size = 3.6, alpha = 0.85) +
  labs(title = "Strain DEGs by directional category",
       x = "log2FC (OTIK vs SH) at 20 C", y = "log2FC (OTIK vs SH) at 26 C") +
  xlim(-max_x, max_x) + ylim(-max_y, max_y) +
  theme_bw(base_size = 12) + theme(plot.title.position = "plot", legend.position = "right") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
ggsave(file.path(results_dir, "24_Annotated_quadrant_scatter.png"), p_quad, width = 12, height = 8, dpi = 300)

save(final, category_counts, gene_sets, file = file.path(results_dir, "Directional_categories_results.RData"))
writeLines(capture.output(sessionInfo()), file.path(results_dir, "sessionInfo_04_directional_categories.txt"))
cat("\n=== Done ===\n")
