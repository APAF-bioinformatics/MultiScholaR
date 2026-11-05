# KIT-093: the statistics behind the notebook's PCA and RLE diagnostics (dev/distributable/plans/
# transcriptomics_workflow.rmd:441-498, 731-732), for the raw, TMM, RUVs and post-RUVIII stages. PCA is stats::prcomp with
# the notebook's defaults (count data are complete, so never pcaWithMissing).

#' PCA statistics of one stage (rmd:472-484).
#' @param mat A genes-by-runs expression matrix.
#' @param design The validated design (joined by `Run`).
#' @return A list: `pca_scores` (the scores of every PC with the design, the notebook's plotted data), `pca_variance`
#'   (every PC's variance fraction, `var` of the plotted scores over their sum, unrounded) and `pca_axis_labels` (the
#'   notebook's labels from `summary(prcomp)`).
#' @noRd
trxPcaStats <- function(mat, design) {
    pca_res <- stats::prcomp(t(mat))
    var_explained <- summary(pca_res)$importance[2, ] * 100
    pca_df <- as.data.frame(pca_res$x) |>
        tibble::rownames_to_column("Run") |>
        dplyr::left_join(design, by = "Run")
    scores <- as.matrix(pca_df[grepl("^PC[0-9]+$", names(pca_df))])
    v <- apply(scores, 2L, stats::var)
    list(pca_scores = pca_df, pca_variance = v / sum(v),
        pca_axis_labels = c(x = sprintf("PC1 (%.2f%%)", var_explained[1]), y = sprintf("PC2 (%.2f%%)", var_explained[2])))
}

#' RLE statistics of one stage (rmd:441-461).
#' @param mat A genes-by-runs expression matrix.
#' @param design The validated design (joined by `Run`).
#' @return A list: `rle_matrix` (the gene-median-centred matrix in long form with the design, the notebook's plotted
#'   data), `rle_whiskers` (per-run quartiles and 1.5 IQR whiskers) and `rle_ylim` (the plot's y limits: the whisker
#'   range plus 5% each side).
#' @noRd
trxRleStats <- function(mat, design) {
    gene_medians <- apply(mat, 1, stats::median)
    rle_matrix <- sweep(mat, 1, gene_medians, "-")
    rle_long <- as.data.frame(rle_matrix) |>
        tibble::rownames_to_column("gene") |>
        tidyr::pivot_longer(cols = !dplyr::all_of("gene"), names_to = "Run", values_to = "rle_value") |>
        dplyr::left_join(design, by = "Run")
    q1 <- tapply(rle_long$rle_value, rle_long$Run, stats::quantile, 0.25, na.rm = TRUE)
    q3 <- tapply(rle_long$rle_value, rle_long$Run, stats::quantile, 0.75, na.rm = TRUE)
    whiskers <- data.frame(Run = names(q1), q1 = unname(q1), q3 = unname(q3), iqr = unname(q3 - q1),
        upper_whisker = unname(q3 + 1.5 * (q3 - q1)), lower_whisker = unname(q1 - 1.5 * (q3 - q1)))
    ylim_vals <- c(min(whiskers$lower_whisker), max(whiskers$upper_whisker))
    buffer <- (ylim_vals[2] - ylim_vals[1]) * 0.05
    list(rle_matrix = rle_long, rle_whiskers = whiskers, rle_ylim = c(ylim_vals[1] - buffer, ylim_vals[2] + buffer))
}
