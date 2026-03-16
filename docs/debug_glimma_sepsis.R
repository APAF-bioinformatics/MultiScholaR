# debug_glimma_sepsis.R
library(MultiScholaR)
library(Glimma)
library(htmlwidgets)
devtools::load_all()

cp_file <- "tests/testdata/sepsis/proteomics/cp08_volcano_input.rds"
input_data <- readRDS(cp_file)

# Replicate internal logic of generateProtDAVolcanoPlotGlimma but minimally
da_results_list <- input_data$da_results_list
da_proteins_long <- da_results_list$da_proteins_long
selected_contrast <- "Sera_vs_RPMI"

plot_data <- da_proteins_long |> dplyr::filter(comparison == selected_contrast)

# 1. SIMPLEST ATTEMPT: x, y, anno (ID only)
message("Attempt 1: x, y, anno (ID only)...")
x <- as.numeric(plot_data$log2FC)
y <- -log10(as.numeric(plot_data$fdr_qvalue))
# Handle zeros
y[is.infinite(y)] <- max(y[!is.infinite(y)]) + 1

anno_simple <- data.frame(
  ID = as.character(plot_data$Protein.Ids),
  Gene = as.character(plot_data$gene_name),
  stringsAsFactors = FALSE
)

widget1 <- glimmaXY(x, y, anno = anno_simple)
saveWidget(widget1, "docs/debug_sepsis_1_anno.html", selfcontained = TRUE)

# 2. ATTEMPT 2: Add counts
message("Attempt 2: Add counts...")
theObject <- da_results_list$theObject
counts_mat <- as.matrix(theObject@protein_quant_table[,-1])
rownames(counts_mat) <- as.character(theObject@protein_quant_table[[1]])
# Sync to plot_data
counts_mat <- counts_mat[as.character(plot_data$Protein.Ids), ]

groups <- factor(theObject@design_matrix$group)

widget2 <- glimmaXY(x, y, anno = anno_simple, counts = counts_mat, groups = groups, transform.counts = "none")
saveWidget(widget2, "docs/debug_sepsis_2_counts.html", selfcontained = TRUE)

message("Done.")
