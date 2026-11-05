# KIT-093: the notebook's edgeR QL differential expression with the RUVs covariates and its annotation
# (dev/distributable/plans/transcriptomics_workflow.rmd:144-148, 524-583) as kernels. The significance column is edgeR's
# Benjamini-Hochberg `FDR` from topTags, kept under that name; no q-value is computed.

#' The notebook's contrast entry: the first non-header, non-empty line of a contrast file (rmd:540-542).
#' @param path The contrast file.
#' @return The entry, `"name = expression"`. Refuses `trx_contrast_unparseable` when there is none.
#' @noRd
trxParseContrastFile <- function(path) {
    contrast_lines <- readLines(path, warn = FALSE)
    contrast_entry <- contrast_lines[!grepl("^contrast", contrast_lines) & nzchar(contrast_lines)][1]
    trxParseContrastEntry(contrast_entry)
    contrast_entry
}

#' Split a contrast entry as the notebook does (rmd:543-545); refuses `trx_contrast_unparseable`.
#' @noRd
trxParseContrastEntry <- function(contrast_entry) {
    refuse <- function() trxAbort(sprintf("the contrast entry '%s' is not 'name = expression'", contrast_entry), c("trx_contrast_unparseable", "trx_error"))
    if (length(contrast_entry) != 1L || is.na(contrast_entry)) refuse()
    contrast_parts <- strsplit(contrast_entry, "=")[[1]]
    if (length(contrast_parts) != 2L) refuse()
    parsed <- list(name = trimws(contrast_parts[1]), expr = trimws(contrast_parts[2]))
    if (!nzchar(parsed$name) || !nzchar(parsed$expr) || inherits(tryCatch(str2lang(parsed$expr), error = function(e) e), "error")) refuse()
    parsed
}

#' The notebook's contrast matrix (rmd:547-548); refuses `trx_contrast_not_estimable` for a non-column or a null contrast.
#' @noRd
trxContrastMatrix <- function(parsed, design_final) {
    refuse <- function(why) trxAbort(sprintf("the contrast '%s' is not estimable: %s", parsed$expr, why), c("trx_contrast_not_estimable", "trx_error"))
    contrasts_matrix <- tryCatch(limma::makeContrasts(contrasts = parsed$expr, levels = design_final),
        error = function(e) refuse(conditionMessage(e)))
    colnames(contrasts_matrix) <- parsed$name
    if (all(contrasts_matrix == 0)) refuse("every coefficient is zero")
    contrasts_matrix
}

#' The notebook's DE model with the RUVs covariates and one contrast (rmd:524-555).
#' @param ruvs_counts,pdata trxRuvs' `ruvs_counts` and `pdata`.
#' @param contrast_entry One `"name = expression"` entry (trxParseContrastFile).
#' @return A list: `design_final`, `design_formula`, `contrast_entry`, `contrasts_matrix`, `fit` (the robust DGEGLM),
#'   `qlf` (the DGELRT), `top_table` (`topTags(n = Inf, sort.by = "PValue")`, BH `FDR`) and `method`. Refuses
#'   `trx_da_rank_deficient`, `trx_contrast_unparseable` and `trx_contrast_not_estimable` before fitting.
#' @noRd
trxDaContrast <- function(ruvs_counts, pdata, contrast_entry) {
    validateTrxCounts(ruvs_counts)
    requireOptionalPackage("edgeR", "the transcriptomics differential expression")
    parsed <- trxParseContrastEntry(contrast_entry)
    dge_final <- edgeR::DGEList(counts = ruvs_counts, group = pdata$group)
    dge_final <- edgeR::calcNormFactors(dge_final, method = "TMM")
    w_terms <- colnames(pdata)[stringr::str_starts(colnames(pdata), "W_")]
    design_formula_str <- paste("~ 0 + group", paste(w_terms, collapse = " + "), sep = " + ")
    design_final <- stats::model.matrix(stats::as.formula(design_formula_str), data = pdata)
    colnames(design_final) <- make.names(colnames(design_final))
    if (qr(design_final)$rank < ncol(design_final)) {
        trxAbort("the final design is not of full column rank (W confounded with group)", c("trx_da_rank_deficient", "trx_error"))
    }
    contrasts_matrix <- trxContrastMatrix(parsed, design_final)
    dge_final <- edgeR::estimateDisp(dge_final, design_final, robust = TRUE)
    fit_final <- edgeR::glmQLFit(dge_final, design_final, robust = TRUE)
    lrt <- edgeR::glmQLFTest(fit_final, contrast = contrasts_matrix[, parsed$name])
    top_genes <- edgeR::topTags(lrt, n = Inf, sort.by = "PValue")$table
    list(design_final = design_final, design_formula = design_formula_str, contrast_entry = contrast_entry,
        contrasts_matrix = contrasts_matrix, fit = fit_final, qlf = lrt, top_table = top_genes, method = "edgeR_QL_robust")
}

#' The notebook's gene annotation of the DE table (rmd:144-148, 558-582).
#' @param top_table trxDaContrast's `top_table`.
#' @param count_table The count table (gene ID and gene name columns).
#' @param comparison The contrast name.
#' @param gene_id_column,annotation_id_column The count table's gene ID and gene name columns.
#' @return A list: `feature_meta` (GeneID, GeneName; GeneName falls back to GeneID when NA or empty), `de_annotated`
#'   (`top_table` joined by GeneID, one match per GeneID, with DisplaySymbol) and `top_de` (top 10 up and 10 down at
#'   FDR < 0.05, with `comparison`).
#' @noRd
trxAnnotateDe <- function(top_table, count_table, comparison, gene_id_column = "GeneID", annotation_id_column = "GeneName") {
    feature_meta <- count_table |>
        dplyr::select(GeneID = dplyr::all_of(gene_id_column), GeneName = dplyr::all_of(annotation_id_column)) |>
        dplyr::mutate(GeneName = ifelse(is.na(.data$GeneName) | .data$GeneName == "", .data$GeneID, .data$GeneName))
    de_annotated <- top_table |>
        tibble::rownames_to_column("GeneID") |>
        dplyr::left_join(feature_meta, by = "GeneID", relationship = "many-to-one") |>
        dplyr::mutate(DisplaySymbol = ifelse(!is.na(.data$GeneName) & .data$GeneName != "", .data$GeneName, .data$GeneID)) |>
        dplyr::relocate(dplyr::all_of(c("GeneID", "GeneName", "DisplaySymbol")))
    top <- function(direction, regulation) {
        de_annotated |>
            dplyr::filter(.data$FDR < 0.05, direction * .data$logFC > 0) |>
            dplyr::arrange(.data$FDR) |>
            utils::head(10) |>
            dplyr::mutate(regulation = regulation)
    }
    top_de <- dplyr::bind_rows(top(1, "Upregulated"), top(-1, "Downregulated")) |>
        dplyr::mutate(comparison = .env$comparison)
    list(feature_meta = feature_meta, de_annotated = de_annotated, top_de = top_de)
}
