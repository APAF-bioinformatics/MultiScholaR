# KIT-093: the notebook's TMM pass on the filtered counts and its two logCPM matrices
# (dev/distributable/plans/transcriptomics_workflow.rmd:178-180, 490) as a kernel.

#' TMM factors and logCPM of the filtered counts.
#' @param filtered_counts The filtered counts (trxExpressionFilter's `filtered_counts`).
#' @param group The design's group factor.
#' @return A list: `norm_factors_filtered` (edgeR's samples table of the filtered DGEList after TMM: `group`,
#'   `lib.size`, `norm.factors`), `log_cpm_tmm` (`cpm(dge_filtered, log = TRUE)`) and `log_cpm_raw` (the notebook's
#'   "raw" diagnostic: logCPM of the filtered counts without TMM).
#' @noRd
trxTmmLogCpm <- function(filtered_counts, group) {
    validateTrxCounts(filtered_counts)
    dge_filtered <- trxFilteredDge(filtered_counts, group)
    list(norm_factors_filtered = dge_filtered$samples, log_cpm_tmm = edgeR::cpm(dge_filtered, log = TRUE),
        log_cpm_raw = edgeR::cpm(edgeR::DGEList(counts = filtered_counts), log = TRUE))
}

#' The notebook's filtered DGEList with recomputed TMM factors (rmd:178-179).
#' @noRd
trxFilteredDge <- function(filtered_counts, group) {
    requireOptionalPackage("edgeR", "the transcriptomics TMM normalisation")
    dge_filtered <- edgeR::DGEList(counts = filtered_counts, group = group)
    edgeR::calcNormFactors(dge_filtered, method = "TMM")
}
