# KIT-093: the notebook's publication renderers (dev/distributable/plans/transcriptomics_workflow.rmd:313-395, 441-513,
# 565-583, 597-666, 674-689, 731-747, 755-775), moved into the package. Each takes kernel products and returns the figure
# with its data (`layers`: ggplot2::layer_data() of every layer; `data`: the plotted data frame) or the tables; none writes
# a file (KIT-095's publication operations do). The plotting code is the notebook's, with qualified calls; the volcano's
# y axis is labelled FDR (BH).

utils::globalVariables(c("K", "cc", "featureset", "Percentage", "Composite_Score", "Best_K", "PC1", "PC2", "Run",
    "rle_value", "significance", "n", "FDR", "logFC", "DisplaySymbol"))

#' layer_data() of every layer of a figure.
#' @noRd
trxPlotLayers <- function(p) lapply(seq_along(p$layers), function(i) ggplot2::layer_data(p, i))

#' A figure with its layers and plotted data.
#' @noRd
trxFigure <- function(p) list(plot = p, layers = trxPlotLayers(p), data = p$data)

#' The panels of a patchwork grid in reading order (a filler, patchwork's empty top level, is not a panel).
#' @noRd
trxGridPanels <- function(p) {
    inner <- if (inherits(p, "patchwork")) unlist(lapply(p$patches$plots, trxGridPanels), recursive = FALSE) else list()
    if (inherits(p, "plot_filler")) inner else c(inner, list(p))
}

#' A patchwork grid with every panel's layers and its title and caption.
#' @noRd
trxGridFigure <- function(grid) {
    list(plot = grid, layers = lapply(trxGridPanels(grid), trxPlotLayers), annotation = grid$patches$annotation[c("title", "caption")])
}

#' The notebook's `format_cancor_plot` (rmd:313-341).
#' @noRd
trxCancorPlot <- function(opt_res, panel_label = "A") {
    df <- opt_res$cancor_data
    best_k_val <- opt_res$best_k
    ggplot2::ggplot(df, ggplot2::aes(x = K, y = cc, color = featureset, group = featureset)) +
        ggplot2::geom_line(linewidth = 1.1) +
        ggplot2::geom_point(size = 3) +
        ggplot2::geom_vline(xintercept = best_k_val, linetype = "dashed", color = "#DC2626", linewidth = 0.8) +
        ggplot2::annotate("text", x = best_k_val + 0.15, y = 0.25,
            label = sprintf("Optimal k = %d\n(Separation = %.3f)", best_k_val, opt_res$separation_score),
            hjust = 0, size = 3.5, fontface = "bold", color = "#DC2626") +
        ggplot2::scale_color_manual(values = c("All" = "#2563EB", "Control" = "#10B981")) +
        ggplot2::scale_x_continuous(breaks = unique(df$K)) +
        ggplot2::coord_cartesian(ylim = c(0, 1.05)) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::labs(title = sprintf("Panel %s: Canonical Correlation Analysis", panel_label),
            subtitle = sprintf("Tested %d%% negative controls (n = %d) selected via ANOVA/GLM", opt_res$percentage, opt_res$n_controls),
            x = "Factor Dimensions (K)", y = expression(paste("Squared Canonical Correlation (", cancor^2, ")")), color = "Feature Set") +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.position = "right")
}

#' The notebook's `format_tuning_curve` (rmd:343-371).
#' @noRd
trxTuningCurve <- function(summary_df, opt_res, panel_label = "B") {
    best_pct_val <- opt_res$percentage
    score_y <- min(summary_df$Composite_Score) + 0.85 * diff(range(summary_df$Composite_Score))
    text_x <- if (best_pct_val >= 20) best_pct_val - 1.5 else best_pct_val + 1.5
    ggplot2::ggplot(summary_df, ggplot2::aes(x = Percentage, y = Composite_Score)) +
        ggplot2::geom_line(color = "#475569", linewidth = 1) +
        ggplot2::geom_point(ggplot2::aes(color = factor(Best_K)), size = 3) +
        ggplot2::geom_vline(xintercept = best_pct_val, linetype = "dashed", color = "#DC2626", linewidth = 0.8) +
        ggplot2::geom_point(data = summary_df[summary_df$Percentage == best_pct_val, ], color = "#DC2626", size = 4.5, shape = 18) +
        ggplot2::annotate("text", x = text_x, y = score_y,
            label = sprintf("Optimal: %d%% controls\nComposite Score = %.3f", best_pct_val, opt_res$composite_score),
            hjust = if (best_pct_val >= 20) 1 else 0, size = 3.2, fontface = "bold", color = "#DC2626") +
        ggplot2::scale_color_brewer(palette = "Set1", name = "Selected k") +
        ggplot2::scale_x_continuous(breaks = summary_df$Percentage) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::labs(title = sprintf("Panel %s: Control Percentage Optimization", panel_label),
            subtitle = "MultiScholaR Composite Score = Separation * (1 - k_penalty)",
            x = "Negative Controls Selected from Invariant Features (%)", y = "Composite Score") +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.position = "right")
}

#' The notebook's CCA panels, their grid and the workbook's two sheets (rmd:373-392).
#' @param best,grid_summary trxRuvGrid's `best` and `grid_summary`.
#' @param total_features The number of filtered genes.
#' @return A list: `cancor` and `tuning` (figures), `grid` (the cowplot grid and its layers) and `workbook`
#'   (`Optimal_Parameters`, `Grid_Optimization`).
#' @noRd
trxCcaPanels <- function(best, grid_summary, total_features) {
    requireOptionalPackage("cowplot", "the transcriptomics CCA panels")
    p_cca_a <- trxCancorPlot(best, "A")
    p_cca_b <- trxTuningCurve(grid_summary, best, "B")
    grid_cca <- cowplot::plot_grid(p_cca_a, p_cca_b, ncol = 2, labels = NULL)
    opt_summary_tbl <- tibble::tibble(
        Parameter = c("Optimal_K", "Optimal_Percentage", "Num_Negative_Controls", "Total_Features", "Separation_Score", "Composite_Score"),
        Value = c(best$best_k, best$percentage, best$n_controls, total_features, round(best$separation_score, 4), round(best$composite_score, 4)))
    list(cancor = trxFigure(p_cca_a), tuning = trxFigure(p_cca_b), grid = list(plot = grid_cca, layers = trxPlotLayers(grid_cca)),
        workbook = list(Optimal_Parameters = opt_summary_tbl, Grid_Optimization = grid_summary))
}

#' The notebook's `generatePcaPlot` on trxPcaStats' products (rmd:472-488).
#' @noRd
trxPcaPlot <- function(stats, title, group_col = "group") {
    ggplot2::ggplot(stats$pca_scores, ggplot2::aes(x = PC1, y = PC2, color = .data[[group_col]])) +
        ggplot2::geom_point(size = 3.5) +
        ggplot2::labs(title = title, x = stats$pca_axis_labels[["x"]], y = stats$pca_axis_labels[["y"]], color = "Group") +
        ggplot2::theme_bw()
}

#' The notebook's `generateRlePlot` on trxRleStats' products (rmd:441-470).
#' @noRd
trxRlePlot <- function(stats, title, group_col = "group") {
    ggplot2::ggplot(stats$rle_matrix, ggplot2::aes(x = Run, y = rle_value, fill = .data[[group_col]])) +
        ggplot2::geom_boxplot(outlier.shape = NA) +
        ggplot2::coord_cartesian(ylim = stats$rle_ylim) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        ggplot2::labs(title = title, x = "Sample", y = "RLE", fill = "Group") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
}

#' PCA and RLE panels of named stages, as figures named `pca_<stage>` then `rle_<stage>`.
#' @noRd
trxStagePanels <- function(mats, titles, design) {
    pca <- lapply(names(mats), function(s) trxFigure(trxPcaPlot(trxPcaStats(mats[[s]], design), paste("PCA:", titles[[s]]))))
    rle <- lapply(names(mats), function(s) trxFigure(trxRlePlot(trxRleStats(mats[[s]], design), paste("RLE:", titles[[s]]))))
    stats::setNames(c(pca, rle), c(paste0("pca_", names(mats)), paste0("rle_", names(mats))))
}

#' The notebook's normalisation-stage grid: PCA and RLE for raw, TMM and RUVs (rmd:490-506).
#' @param log_cpm_raw,log_cpm_tmm,log_cpm_ruvs The three stages' genes-by-runs matrices.
#' @param design The validated design.
#' @param best_k The RUVs k (in the titles).
#' @return A list: `panels` (six figures) and `grid` (trxGridFigure).
#' @noRd
trxNormStageGrid <- function(log_cpm_raw, log_cpm_tmm, log_cpm_ruvs, design, best_k) {
    panels <- trxStagePanels(list(raw = log_cpm_raw, tmm = log_cpm_tmm, ruvs = log_cpm_ruvs),
        list(raw = "Raw logCPM", tmm = "After TMM", ruvs = sprintf("After RUVs (k=%d)", best_k)), design)
    p <- lapply(panels, `[[`, "plot")
    grid <- (p$pca_raw | p$pca_tmm | p$pca_ruvs) / (p$rle_raw | p$rle_tmm | p$rle_ruvs) +
        patchwork::plot_annotation(title = "Comparison of Transcriptomics Normalization Stages",
            caption = "Top Row: PCA Plots. Bottom Row: RLE Plots.") +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "bottom")
    list(panels = panels, grid = trxGridFigure(grid))
}

#' The notebook's before/after RUVIII grid (rmd:731-740).
#' @param log_cpm_tmm The TMM logCPM (genes by runs).
#' @param newY trxRuviii's `newY` (runs by genes).
#' @param design The validated design.
#' @return A list: `panels` (pca_tmm, pca_post_ruviii, rle_tmm, rle_post_ruviii) and `grid`.
#' @noRd
trxRuviiiGrid <- function(log_cpm_tmm, newY, design) {
    panels <- trxStagePanels(list(tmm = log_cpm_tmm, post_ruviii = t(newY)),
        list(tmm = "After TMM", post_ruviii = "After RUVIII Normalization"), design)
    p <- lapply(panels, `[[`, "plot")
    grid <- (p$pca_tmm | p$pca_post_ruviii) / (p$rle_tmm | p$rle_post_ruviii) +
        patchwork::plot_annotation(title = "Pre vs Post RUVIII Normalization Diagnostics",
            caption = "Left: After TMM (Pre-RUVIII). Right: After RUVIII.") +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "bottom")
    list(panels = panels, grid = trxGridFigure(grid))
}

#' The notebook's significance classes (rmd:598-606, 631-639).
#' @noRd
trxSignificance <- function(de_df, fdr_thresh) {
    de_df |>
        dplyr::mutate(
            significance = dplyr::case_when(
                FDR < fdr_thresh & logFC > 0 ~ "Upregulated",
                FDR < fdr_thresh & logFC < 0 ~ "Downregulated",
                TRUE ~ "Not Significant"),
            significance = factor(significance, levels = c("Upregulated", "Downregulated", "Not Significant")))
}

TRX_SIGNIFICANCE_COLOURS <- c("Upregulated" = "#EF4444", "Downregulated" = "#3B82F6", "Not Significant" = "#9CA3AF")

#' The notebook's volcano (rmd:597-628): the top 20 significant genes by FDR labelled.
#' @param de_df trxAnnotateDe's `de_annotated`.
#' @param title The contrast name.
#' @return A figure (`plot`, `layers`, `data`).
#' @noRd
trxVolcano <- function(de_df, title, fdr_thresh = 0.05, top_n = 20) {
    df <- trxSignificance(de_df, fdr_thresh)
    labels_df <- df |>
        dplyr::filter(significance != "Not Significant") |>
        dplyr::arrange(FDR) |>
        utils::head(top_n)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = logFC, y = -log10(FDR), color = significance)) +
        ggplot2::geom_point(alpha = 0.6, size = 1.6) +
        ggplot2::scale_color_manual(values = TRX_SIGNIFICANCE_COLOURS) +
        ggplot2::geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed", color = "black") +
        ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey40") +
        ggrepel::geom_text_repel(data = labels_df, ggplot2::aes(label = DisplaySymbol), size = 3.2, max.overlaps = 25) +
        ggplot2::labs(title = paste("Volcano Plot:", title),
            subtitle = sprintf("Highlighted top %d significant genes by FDR (BH)", nrow(labels_df)),
            x = expression(Log[2] ~ Fold ~ Change), y = expression(-Log[10] ~ (FDR ~ (BH))), color = "Significance") +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(legend.position = "bottom")
    trxFigure(p)
}

#' The notebook's significance bar plot (rmd:630-653), zero-count classes kept (`.drop = FALSE`).
#' @inheritParams trxVolcano
#' @return A figure (`plot`, `layers`, `data`).
#' @noRd
trxSignificanceBar <- function(de_df, title, fdr_thresh = 0.05) {
    df <- trxSignificance(de_df, fdr_thresh) |> dplyr::count(significance, .drop = FALSE)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = significance, y = n, fill = significance)) +
        ggplot2::geom_bar(stat = "identity", width = 0.55) +
        ggplot2::geom_text(ggplot2::aes(label = n), vjust = -0.5, fontface = "bold", size = 3.5) +
        ggplot2::scale_fill_manual(values = TRX_SIGNIFICANCE_COLOURS) +
        ggplot2::labs(title = paste("Significance Summary:", title), x = "Class", y = "Gene Count") +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(legend.position = "none")
    trxFigure(p)
}

#' The notebook's DE tables (rmd:565-583): the full annotated table (the CSV's content) and the top genes (the TSV's).
#' @noRd
trxDeTables <- function(de_annotated, top_de) list(de_table = de_annotated, top_de = top_de)

#' The notebook's Glimma volcano inputs (rmd:674-686), the annotation aligned to `qlf`'s feature order (the order
#' glimmaVolcano binds it in). KIT-095 writes the HTML and records its outcome.
#' @param de_annotated trxAnnotateDe's `de_annotated`.
#' @param qlf trxDaContrast's `qlf`.
#' @param contrast_name The contrast name.
#' @return A list: `anno`, `display_columns`, `main`, `lfc`, `fdr`.
#' @noRd
trxGlimmaInputs <- function(de_annotated, qlf, contrast_name) {
    glimma_anno <- de_annotated |>
        dplyr::select(dplyr::all_of(c("GeneID", "GeneName", "DisplaySymbol"))) |>
        tibble::column_to_rownames("GeneID")
    list(anno = glimma_anno[rownames(qlf$table), , drop = FALSE], display_columns = c("GeneID", "GeneName", "DisplaySymbol"),
        main = paste(contrast_name, "(RUVs QL-FTest)"), lfc = 1, fdr = 0.05)
}

#' The notebook's summary text (rmd:756-775): the lines and the grid table appended under them.
#' @param best,grid_summary trxRuvGrid's `best` and `grid_summary`.
#' @param total_features The number of filtered genes.
#' @param design_formula,contrast_entry trxDaContrast's products.
#' @param de_annotated trxAnnotateDe's `de_annotated`.
#' @return A list: `lines` (the second is the timestamp) and `grid_table`.
#' @noRd
trxSummaryText <- function(best, total_features, design_formula, contrast_entry, de_annotated, grid_summary) {
    lines <- c(
        "=== MultiScholaR RUV Parameter Optimization Summary ===",
        paste("Timestamp:", Sys.time()),
        "",
        paste("Dataset: GSE152840 (Sepsis Transcriptomics Workshop)"),
        paste("Total Features Analyzed:", total_features),
        paste("Optimal In silico Negative Controls:", best$percentage, "% (n =", best$n_controls, ")"),
        paste("Optimal Factor Dimensions (k):", best$best_k),
        paste("CCA Separation Score:", round(best$separation_score, 4)),
        paste("MultiScholaR Composite Score:", round(best$composite_score, 4)),
        paste("Dynamic edgeR Model Formula:", design_formula),
        paste("Evaluated Contrast:", contrast_entry),
        paste("Differentially Expressed Genes (FDR < 0.05):", sum(de_annotated$FDR < 0.05)),
        "",
        "--- Full Grid Search Summary ---")
    list(lines = lines, grid_table = grid_summary)
}
