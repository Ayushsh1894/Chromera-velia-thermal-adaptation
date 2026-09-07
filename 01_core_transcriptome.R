# =============================================================================
# Core Transcriptome Analysis (Presence / Absence)
# Chromera velia: SH and OTIK at 20 C and 26 C  |  Sharma et al.
# =============================================================================
# GOAL: Identify which genes are detectably expressed in each condition,
# independent of any reference strain. Presence/absence, NOT differential expr.
#
# METHOD: Convert raw counts to CPM (edgeR). A gene is "expressed in a condition"
# if CPM >= 1 in ALL 3 replicates (primary). A 2-of-3 version is a sensitivity
# check.
#
# INPUTS (working directory)
#   Count_matrix.xlsx            - gene x sample raw counts; first col "Geneid"
#   experimental_design_temp.csv - sample, strain, temperature
#
# OUTPUTS -> Results/Core_transcriptome/  (gene lists, presence/absence matrix,
#   summary tables, bar chart, 4-way Venn, UpSet, gene-sharing histogram)
#
# NORMALISATION NOTE: edgeR CPM here (presence/absence only); DESeq2 for DE.
# =============================================================================

required_cran <- c("readxl", "edgeR", "ggplot2", "dplyr", "tibble", "UpSetR", "ggvenn")
for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "edgeR") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install("edgeR", update = FALSE, ask = FALSE)
    } else install.packages(pkg)
  }
}
library(readxl); library(edgeR); library(ggplot2); library(dplyr)
library(tibble); library(UpSetR); library(ggvenn)

project_dir <- "."
setwd(project_dir)
results_dir <- file.path(project_dir, "Results")
core_dir    <- file.path(results_dir, "Core_transcriptome")
if (!dir.exists(core_dir)) dir.create(core_dir, recursive = TRUE)

counts <- read_excel("Count_matrix.xlsx") |> as.data.frame()
rownames(counts) <- counts$Geneid; counts$Geneid <- NULL
counts <- as.matrix(counts); mode(counts) <- "integer"

metadata <- read.csv("experimental_design_temp.csv", stringsAsFactors = FALSE)
rownames(metadata) <- metadata$sample
counts <- counts[, rownames(metadata)]
stopifnot(all(colnames(counts) == rownames(metadata)))
metadata$condition <- paste(metadata$strain, metadata$temperature, sep = "_")
conditions <- unique(metadata$condition)
cat("Conditions:", paste(conditions, collapse = ", "), "\n")
cat("Total genes:", nrow(counts), "\n\n")

cpm_matrix <- cpm(counts)
CPM_THRESHOLD <- 1
get_expressed_genes <- function(condition_name, min_replicates) {
  samples <- rownames(metadata)[metadata$condition == condition_name]
  pass_count <- rowSums(cpm_matrix[, samples, drop = FALSE] >= CPM_THRESHOLD)
  names(pass_count)[pass_count >= min_replicates]
}
expressed_strict  <- lapply(conditions, get_expressed_genes, min_replicates = 3)
names(expressed_strict) <- conditions
expressed_relaxed <- lapply(conditions, get_expressed_genes, min_replicates = 2)
names(expressed_relaxed) <- conditions

core_genes <- Reduce(intersect, expressed_strict)
pan_genes  <- Reduce(union,     expressed_strict)
cat("Core transcriptome:", length(core_genes), "genes\n")
cat("Pan-transcriptome:", length(pan_genes), "genes\n\n")

condition_specific <- list()
for (cond in conditions) {
  others_union <- Reduce(union, expressed_strict[setdiff(conditions, cond)])
  condition_specific[[cond]] <- setdiff(expressed_strict[[cond]], others_union)
}

for (cond in conditions)
  write.csv(data.frame(Geneid = expressed_strict[[cond]]),
            file.path(core_dir, paste0("expressed_", cond, "_strict.csv")), row.names = FALSE)
write.csv(data.frame(Geneid = core_genes), file.path(core_dir, "core_transcriptome.csv"), row.names = FALSE)
write.csv(data.frame(Geneid = pan_genes),  file.path(core_dir, "pan_transcriptome.csv"), row.names = FALSE)
for (cond in conditions)
  write.csv(data.frame(Geneid = condition_specific[[cond]]),
            file.path(core_dir, paste0("specific_to_", cond, ".csv")), row.names = FALSE)

membership <- data.frame(Geneid = pan_genes)
for (cond in conditions) membership[[cond]] <- as.integer(membership$Geneid %in% expressed_strict[[cond]])
membership$n_conditions <- rowSums(membership[, conditions])
write.csv(membership, file.path(core_dir, "gene_presence_absence_matrix.csv"), row.names = FALSE)

summary_df <- data.frame(
  Condition = conditions,
  Expressed_strict  = sapply(conditions, function(c) length(expressed_strict[[c]])),
  Expressed_relaxed = sapply(conditions, function(c) length(expressed_relaxed[[c]])),
  Specific_strict   = sapply(conditions, function(c) length(condition_specific[[c]])))
summary_extra <- data.frame(
  Condition = c("CORE (all 4)", "PAN (any)"),
  Expressed_strict  = c(length(core_genes), length(pan_genes)),
  Expressed_relaxed = c(length(Reduce(intersect, expressed_relaxed)), length(Reduce(union, expressed_relaxed))),
  Specific_strict   = c(NA, NA))
write.csv(rbind(summary_df, summary_extra),
          file.path(core_dir, "core_transcriptome_summary.csv"), row.names = FALSE)

bar_data <- summary_df |> select(Condition, Expressed_strict) |>
  mutate(Condition = factor(Condition, levels = sort(conditions)))
p_bar <- ggplot(bar_data, aes(x = Condition, y = Expressed_strict, fill = Condition)) +
  geom_col(width = 0.65) + geom_text(aes(label = Expressed_strict), vjust = -0.5, size = 4.5) +
  geom_hline(yintercept = length(core_genes), linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("#9ecae1", "#1f77b4", "#f4a582", "#d62728")) +
  labs(title = "Expressed genes per condition",
       subtitle = paste0("CPM >= 1 in all 3 replicates. Pan-transcriptome = ", length(pan_genes), " genes."),
       x = NULL, y = "Number of expressed genes") +
  theme_bw(base_size = 12) + theme(legend.position = "none", plot.title.position = "plot") +
  expand_limits(y = max(bar_data$Expressed_strict) * 1.1)
ggsave(file.path(core_dir, "fig1_expressed_genes_barchart.png"), p_bar, width = 8, height = 6, dpi = 300)

p_venn <- ggvenn(expressed_strict, fill_color = c("#9ecae1", "#1f77b4", "#f4a582", "#d62728"),
                 stroke_size = 0.4, set_name_size = 4.5, text_size = 3.2, show_percentage = FALSE) +
  ggtitle("Core transcriptome - expressed genes across 4 conditions",
          subtitle = paste0("Center (all 4) = core = ", length(core_genes), " genes"))
ggsave(file.path(core_dir, "fig2_core_transcriptome_venn4.png"), p_venn, width = 9, height = 8, dpi = 300)

png(file.path(core_dir, "fig3_core_transcriptome_upset.png"), width = 12, height = 7, units = "in", res = 300)
upset(fromList(expressed_strict), order.by = "freq", nsets = 4, nintersects = NA,
      sets.bar.color = c("#9ecae1", "#1f77b4", "#f4a582", "#d62728"),
      main.bar.color = "#34495e", matrix.color = "#34495e",
      mainbar.y.label = "Genes in intersection", sets.x.label = "Expressed genes per condition",
      text.scale = c(1.5, 1.3, 1.3, 1.1, 1.4, 1.2), point.size = 3, line.size = 1)
dev.off()

dist_data <- membership |> count(n_conditions, name = "n_genes")
write.csv(dist_data, file.path(core_dir, "gene_sharing_distribution.csv"), row.names = FALSE)
p_dist <- ggplot(dist_data, aes(x = factor(n_conditions), y = n_genes, fill = factor(n_conditions))) +
  geom_col(width = 0.65) + geom_text(aes(label = n_genes), vjust = -0.5, size = 4.5) +
  scale_fill_brewer(palette = "YlGnBu") +
  labs(title = "Gene sharing across conditions",
       x = "Number of conditions a gene is expressed in", y = "Number of genes") +
  theme_bw(base_size = 12) + theme(legend.position = "none", plot.title.position = "plot") +
  expand_limits(y = max(dist_data$n_genes) * 1.1)
ggsave(file.path(core_dir, "fig4_gene_sharing_distribution.png"), p_dist, width = 8, height = 6, dpi = 300)

save(expressed_strict, expressed_relaxed, core_genes, pan_genes, condition_specific, membership,
     file = file.path(core_dir, "Core_transcriptome_results.RData"))
writeLines(capture.output(sessionInfo()), file.path(core_dir, "sessionInfo_01_core_transcriptome.txt"))
cat("=== Core transcriptome analysis complete ===\n")
