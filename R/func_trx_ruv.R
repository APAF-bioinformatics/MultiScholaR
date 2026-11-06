# KIT-093: the notebook's RUV stages (dev/distributable/plans/transcriptomics_workflow.rmd:195-298 grid search, 406-426
# RUVs, 704-712 RUVIII) as side-effect-free kernels. The edgeR, RUVSeq and ruv calls are the notebook's, unchanged; there
# is no RNG on this path (ruv_cancorplot is svd plus cancor; RUVs, RUVIII with a fixed k and edgeR are deterministic).

TRX_CONTROL_SOURCE <- "inferred_least_significant"

#' The notebook's RUV parameter grid search.
#' @param filtered_counts,log_cpm_tmm The filtered counts and their TMM logCPM (trxTmmLogCpm).
#' @param design The validated design (`Run`, `group`, design order).
#' @param params Optional list: `percentage_range` (default `seq(5, 30, by = 5)`), `k_penalty_weight` (0.05).
#' @return A list: `prelim_order` (genes ordered least significant first by the non-robust preliminary QL fit),
#'   `max_acceptable_k`, `grid_members` (the successful members), `grid_failures` (percentage and condition class of each
#'   failed member; a declared addition), `grid_summary` (the notebook's 4-dp `df_summary`, display only), `best` (the
#'   first maximum of the unrounded composite score), `control_genes` (`best`'s) and `control_source`.
#' @noRd
trxRuvGrid <- function(filtered_counts, log_cpm_tmm, design, params = list()) {
    validateTrxCounts(filtered_counts)
    params <- utils::modifyList(list(percentage_range = seq(5, 30, by = 5), k_penalty_weight = 0.05), params)
    total_features <- nrow(log_cpm_tmm)
    if (total_features < 10L) {
        trxAbort(sprintf("%d genes are fewer than the grid's minimum of 10 controls", total_features), c("trx_grid_too_few_features", "trx_error"))
    }
    design_prelim <- stats::model.matrix(~ group, data = design)
    prelim_order <- trxPrelimOrder(trxFilteredDge(filtered_counts, design$group), design_prelim)
    max_acceptable_k <- if (ncol(filtered_counts) < 15) 3L else 4L
    setting <- list(log_cpm_tmm = log_cpm_tmm, order = prelim_order, X_bio = design_prelim[, -1, drop = FALSE],
        max_k = max_acceptable_k, weight = params$k_penalty_weight)
    results <- lapply(params$percentage_range, trxGridMember, setting = setting)
    failed <- vapply(results, function(r) is.null(r$percentage), logical(1))
    members <- results[!failed]
    if (!length(members)) trxAbort("every RUV grid member failed", c("trx_grid_no_member", "trx_error"))
    best <- members[[trxSelectGridWinner(members)]]
    list(prelim_order = prelim_order, max_acceptable_k = max_acceptable_k, grid_members = members,
        grid_failures = data.frame(percentage = params$percentage_range[failed],
            condition_class = vapply(results[failed], `[[`, "", "condition_class")),
        grid_summary = trxGridSummary(members), best = best, control_genes = best$control_genes, control_source = TRX_CONTROL_SOURCE)
}

#' The preliminary non-robust QL ranking (rmd:232-237): least significant first. The test is over every group coefficient of
#' the `~ group` design (KIT-093 revision 1.2: across all groups, as proteomics' ANOVA control selection); with two groups
#' that is the notebook's `coef = 2`.
#' @noRd
trxPrelimOrder <- function(dge_filtered, design_prelim) {
    dge_prelim <- edgeR::estimateDisp(dge_filtered, design_prelim)
    fit_prelim <- edgeR::glmQLFit(dge_prelim, design_prelim)
    lrt_prelim <- edgeR::glmQLFTest(fit_prelim, coef = 2:ncol(design_prelim))
    tt_prelim <- edgeR::topTags(lrt_prelim, n = nrow(dge_filtered), sort.by = "none")$table
    order(tt_prelim$PValue, decreasing = TRUE)
}

#' The CCA call of a grid member (the seam a test replaces to force a member failure).
#' @noRd
trxCancor <- function(Y, X, ctl) ruv::ruv_cancorplot(Y = Y, X = X, ctl = ctl)

#' One grid member (rmd:248-280); a failed member returns its condition class instead.
#' @noRd
trxGridMember <- function(pct, setting) {
    total_features <- nrow(setting$log_cpm_tmm)
    max_acceptable_k <- setting$max_k
    n_ctrl <- max(10, floor(total_features * (pct / 100)))
    ctl_idx <- rep(FALSE, total_features)
    ctl_idx[setting$order[seq_len(n_ctrl)]] <- TRUE
    cancor_p <- tryCatch(trxCancor(Y = t(setting$log_cpm_tmm), X = setting$X_bio, ctl = ctl_idx), error = function(e) e)
    if (inherits(cancor_p, "error")) return(list(condition_class = class(cancor_p)[[1L]]))
    if (is.null(cancor_p) || is.null(cancor_p$data)) return(list(condition_class = "trx_grid_member_no_data"))
    best_k <- trxFindBestK(cancor_p, max_k = max_acceptable_k)
    sep_score <- trxSeparationScore(cancor_p, "max_difference", max_k = max_acceptable_k)
    if (best_k <= max_acceptable_k) {
        k_penalty <- if (max_acceptable_k <= 1) 0 else setting$weight * (best_k - 1) / (max_acceptable_k - 1)
    } else {
        excess_k <- best_k - max_acceptable_k
        k_penalty <- setting$weight + (1 - setting$weight) * (1 - exp(-excess_k))
    }
    k_penalty <- max(0, min(1, k_penalty))
    list(percentage = pct, n_controls = n_ctrl, best_k = best_k, separation_score = sep_score, k_penalty = k_penalty,
        composite_score = sep_score * (1 - k_penalty), cancor_data = cancor_p$data, ctl_idx = ctl_idx,
        control_genes = rownames(setting$log_cpm_tmm)[ctl_idx])
}

#' The notebook's `find_best_k` (rmd:195-210).
#' @noRd
trxFindBestK <- function(cancor_p, max_k = 3L) {
    diffs <- trxCancorDiffs(cancor_p, max_k)
    if (is.null(diffs)) return(1L)
    valid_diffs <- diffs$diffs[is.finite(diffs$diffs)]
    if (length(valid_diffs) == 0) return(1L)
    max_diff <- max(valid_diffs, na.rm = TRUE)
    best_idx <- which(diffs$diffs == max_diff)
    as.integer(diffs$common_k[best_idx[1]])
}

#' The notebook's `calc_separation_score` (rmd:212-229).
#' @noRd
trxSeparationScore <- function(cancor_p, metric = "max_difference", max_k = 3L) {
    diffs <- trxCancorDiffs(cancor_p, max_k)
    if (is.null(diffs)) return(0)
    valid_diffs <- diffs$diffs[is.finite(diffs$diffs)]
    if (length(valid_diffs) == 0) return(0)
    switch(metric,
        "max_difference" = max(valid_diffs, na.rm = TRUE),
        "mean_difference" = mean(valid_diffs, na.rm = TRUE),
        max(valid_diffs, na.rm = TRUE))
}

#' The All-minus-Control canonical correlation differences shared by both notebook helpers; NULL without a common K.
#' @noRd
trxCancorDiffs <- function(cancor_p, max_k) {
    df <- cancor_p$data
    if (!is.null(max_k)) df <- df[df$K <= max_k, , drop = FALSE]
    all_data <- df[df$featureset == "All", ]
    ctrl_data <- df[df$featureset == "Control", ]
    common_k <- intersect(all_data$K, ctrl_data$K)
    if (length(common_k) == 0) return(NULL)
    list(common_k = common_k, diffs = all_data$cc[match(common_k, all_data$K)] - ctrl_data$cc[match(common_k, ctrl_data$K)])
}

#' The winner: the first maximum of the unrounded composite score (rmd:295).
#' @noRd
trxSelectGridWinner <- function(members) which.max(purrr::map_dbl(members, "composite_score"))

#' The notebook's 4-dp rounded grid summary (rmd:283-292), for display.
#' @noRd
trxGridSummary <- function(members) {
    purrr::map_dfr(members, function(r) {
        tibble::tibble(Percentage = r$percentage, Num_Controls = r$n_controls, Best_K = r$best_k,
            Separation_Score = round(r$separation_score, 4), K_Penalty = round(r$k_penalty, 4),
            Composite_Score = round(r$composite_score, 4))
    })
}

#' The notebook's RUVs normalisation (rmd:406-426).
#' @param filtered_counts The filtered counts (whole numbers; the notebook's `round()` is kept and is a no-op).
#' @param log_cpm_tmm Their TMM logCPM.
#' @param design The validated design.
#' @param control_genes,k The grid's control genes and best k.
#' @return A list: `W` (`W_1..W_k` from `pData`), `ruvs_counts` (`counts(set_ruvs)`), `log_cpm_ruvs` (visualisation only:
#'   `removeBatchEffect` with W, never the edgeR input), `pdata` (`pData(set_ruvs)`, the DA kernel's input) and
#'   `control_source`.
#' @noRd
trxRuvs <- function(filtered_counts, log_cpm_tmm, design, control_genes, k) {
    validateTrxCounts(filtered_counts)
    for (package in c("EDASeq", "RUVSeq", "Biobase", "BiocGenerics")) requireOptionalPackage(package, "the transcriptomics RUVs normalisation")
    set <- EDASeq::newSeqExpressionSet(counts = as.matrix(round(filtered_counts)),
        phenoData = data.frame(design, row.names = design$Run))
    differences <- RUVSeq::makeGroups(design$group)
    # RUVSeq's SeqExpressionSet method calls counts() unqualified and finds it only with BiocGenerics on the search path
    # (the notebook attaches it through RUVSeq's Depends), so BiocGenerics is attached for this call and detached after.
    set_ruvs <- withr::with_package("BiocGenerics", RUVSeq::RUVs(set, cIdx = control_genes, k = k, scIdx = differences))
    ruvs_pdata <- Biobase::pData(set_ruvs)
    w_matrix <- ruvs_pdata[, paste0("W_", seq_len(k)), drop = FALSE]
    design_groups_only <- stats::model.matrix(~ group, data = ruvs_pdata)
    log_cpm_corrected <- limma::removeBatchEffect(log_cpm_tmm, covariates = w_matrix, design = design_groups_only)
    list(W = w_matrix, ruvs_counts = BiocGenerics::counts(set_ruvs), log_cpm_ruvs = log_cpm_corrected, pdata = ruvs_pdata,
        control_source = TRX_CONTROL_SOURCE)
}

#' RUVIII on the TMM logCPM (rmd:704-712), `M` from the group (the notebook's default) or from technical replicates.
#' @param log_cpm_tmm The TMM logCPM (genes by runs).
#' @param design The validated design.
#' @param control_genes,k The grid's control genes and best k.
#' @param m_source `"group"` or `"technical_replicate"`.
#' @param technical_replicate_id The design column holding the biological-sample ID technical replicates share.
#' @return A list: `newY` (runs by genes), `M`, `m_source`, `m_label`, `control_genes`, `ctl_mask`, `control_source`, `k`.
#' @noRd
trxRuviii <- function(log_cpm_tmm, design, control_genes, k, m_source = c("group", "technical_replicate"),
                      technical_replicate_id = NA_character_) {
    m_source <- match.arg(m_source)
    M <- trxRuviiiM(design, m_source, technical_replicate_id)
    ctl_mask <- rep(FALSE, nrow(log_cpm_tmm))
    ctl_mask[rownames(log_cpm_tmm) %in% control_genes] <- TRUE
    res_ruviii <- ruv::RUVIII(Y = t(log_cpm_tmm), M = M, ctl = ctl_mask, k = k, return.info = TRUE)
    list(newY = res_ruviii$newY, M = M, m_source = m_source,
        m_label = if (m_source == "group") "M from treatment group (notebook default)" else "M from technical replicates",
        control_genes = control_genes, ctl_mask = ctl_mask, control_source = TRX_CONTROL_SOURCE, k = k)
}

#' The RUVIII replicate matrix. Technical-replicate sets must exist, be complete, stay within one group and include a set
#' of at least two runs (so the biological `replicates` number, shared across groups, is refused).
#' @noRd
trxRuviiiM <- function(design, m_source, technical_replicate_id) {
    if (m_source == "group") return(ruv::replicate.matrix(as.character(design$group)))
    refuse <- function(message) trxAbort(message, c("trx_ruviii_replicates_invalid", "trx_error"))
    column <- technical_replicate_id
    if (length(column) != 1L || is.na(column) || !column %in% names(design)) refuse("the technical-replicate column is absent")
    sets <- as.character(design[[column]])
    if (anyNA(sets)) refuse(sprintf("the technical-replicate column '%s' has NA", column))
    spans <- tapply(as.character(design$group), sets, function(g) length(unique(g)))
    if (any(spans > 1L)) refuse(sprintf("replicate set(s) span two groups: %s", paste(names(spans)[spans > 1L], collapse = ", ")))
    if (max(table(sets)) < 2L) refuse("no technical-replicate set has two or more runs")
    ruv::replicate.matrix(sets)
}
