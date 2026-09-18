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

T26_df <- as.data.frame(res_T26) |> rownames_to_column("Geneid") |>
  select(Geneid, log2FC_T26 = log2FoldChange, padj_T26 = padj)
T20_df <- as.data.frame(res_T20) |> rownames_to_column("Geneid") |>
  select(Geneid, log2FC_T20 = log2FoldChange, padj_T20 = padj)

combined <- inner_join(T26_df, T20_df, by = "Geneid") |>
  filter(!is.na(log2FC_T26) & !is.na(log2FC_T20)) |>
  mutate(
    # Directional sig flags - THE FIX: no abs(), direction checked explicitly
    sig_T26_OTIup = !is.na(padj_T26) & padj_T26 < 0.05 & log2FC_T26 >  1,
    sig_T26_SHup  = !is.na(padj_T26) & padj_T26 < 0.05 & log2FC_T26 < -1,
    sig_T20_OTIup = !is.na(padj_T20) & padj_T20 < 0.05 & log2FC_T20 >  1,
    sig_T20_SHup  = !is.na(padj_T20) & padj_T20 < 0.05 & log2FC_T20 < -1,
    sig_T26_any = sig_T26_OTIup | sig_T26_SHup,
    sig_T20_any = sig_T20_OTIup | sig_T20_SHup,
    category = case_when(
      sig_T26_any & sig_T20_any & sig_T26_OTIup & sig_T20_OTIup ~ "Sig in both, OTI higher (constitutive)",
      sig_T26_any & sig_T20_any & sig_T26_SHup  & sig_T20_SHup  ~ "Sig in both, SH higher (constitutive)",
      sig_T26_any & sig_T20_any                                  ~ "Sig in both, direction reversal",
      sig_T26_any & !sig_T20_any & sig_T26_OTIup                 ~ "Heat-specific, OTI higher",
      sig_T26_any & !sig_T20_any & sig_T26_SHup                  ~ "Heat-specific, SH higher",
      sig_T20_any & !sig_T26_any & sig_T20_SHup                  ~ "Cool-specific, SH higher",
      sig_T20_any & !sig_T26_any & sig_T20_OTIup                 ~ "Cool-specific, OTI higher",
      TRUE ~ "Not significant"),
    # collapse for plotting: direction split within "sig in both" isn't the point of
    # this panel (position on the diagonal already shows it), so merge for color,
    # but keep the reversal flag separate since that IS a distinct, rare story
    category_plot = case_when(
      category %in% c("Sig in both, OTI higher (constitutive)",
                      "Sig in both, SH higher (constitutive)") ~ "Sig in both (constitutive)",
      TRUE ~ category),
    strain_diff_change = abs(log2FC_T26 - log2FC_T20))

write.csv(combined, file.path(results_dir, "17_Strain_DEGs_combined_T26_T20.csv"), row.names = FALSE)

cat("\n=== Corrected category counts ===\n")
print(table(combined$category_plot))

# Correlation stats for the figure and text (Tomas's point 3 - previously missing)
spearman_test <- cor.test(combined$log2FC_T20, combined$log2FC_T26, method = "spearman")
pearson_test  <- cor.test(combined$log2FC_T20, combined$log2FC_T26, method = "pearson")
r2_val <- pearson_test$estimate^2
cat(sprintf("\nSpearman rho = %.3f (p = %.2e)\n", spearman_test$estimate, spearman_test$p.value))
cat(sprintf("Pearson r = %.3f, R2 = %.3f (p = %.2e)\n", pearson_test$estimate, r2_val, pearson_test$p.value))
write.csv(data.frame(statistic = c("Spearman_rho","Pearson_r","R2"),
                     value = c(spearman_test$estimate, pearson_test$estimate, r2_val)),
          file.path(results_dir, "17b_Strain_DEGs_correlation_stats.csv"), row.names = FALSE)

top_changers <- combined |>
  filter(category_plot == "Heat-specific, OTI higher") |>
  arrange(desc(strain_diff_change)) |> head(20)

cat_colors <- c(
  "Not significant"                        = "#ececec",
  "Sig in both (constitutive)"             = "#6e6e6e",
  "Sig in both, direction reversal"        = "#8e44ad",
  "Heat-specific, OTI higher"              = "#d62728",
  "Heat-specific, SH higher"               = "#e8964f",
  "Cool-specific, SH higher"               = "#1f77b4",
  "Cool-specific, OTI higher"              = "#4fc3c9")

stats_label <- sprintf("Spearman rho = %.3f\nPearson r = %.3f\nR2 = %.3f\nn = %s genes",
                       spearman_test$estimate, pearson_test$estimate, r2_val,
                       format(nrow(combined), big.mark = ","))

p_scatter <- ggplot(combined, aes(x = log2FC_T20, y = log2FC_T26, color = category_plot)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey30") +
  geom_point(alpha = 0.6, size = 1.4) +
  geom_text_repel(data = top_changers, aes(label = Geneid), size = 3, color = "black",
                  max.overlaps = 30, box.padding = 0.4, segment.color = "grey50") +
  annotate("label", x = -Inf, y = Inf, label = stats_label, hjust = -0.05, vjust = 1.05,
           size = 3.2, label.size = 0.3, fill = "white", alpha = 0.85) +
  scale_color_manual(values = cat_colors) +
  labs(title = "Strain difference (OTI vs SH) across temperatures",
       subtitle = "Each gene's log2FC at 20°C vs 26°C. Positive = higher in OTI. Categories are sign-matched.",
       x = "log2FC (OTI vs SH) at 20°C", y = "log2FC (OTI vs SH) at 26°C", color = "Category") +
  theme_bw(base_size = 12) + theme(legend.position = "right")

ggsave(file.path(results_dir, "18_Strain_DEGs_T26vsT20_scatter.png"), p_scatter, width = 12, height = 7.5, dpi = 300)

cat("\n=== Done. Replaces old 03_strain_comparison.R section. ===\n")