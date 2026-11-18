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
    order(trxPrelimPValues(dge_filtered, design_prelim), decreasing = TRUE)
}

#' The p-values of the preliminary non-robust QL test, in gene order: what `trxPrelimOrder` ranks, kept for the
#' first-pass control rule (KIT-129 W4), which drops the most significant quarter.
#' @noRd
trxPrelimPValues <- function(dge_filtered, design_prelim) {
    dge_prelim <- edgeR::estimateDisp(dge_filtered, design_prelim)
    fit_prelim <- edgeR::glmQLFit(dge_prelim, design_prelim)
    lrt_prelim <- edgeR::glmQLFTest(fit_prelim, coef = 2:ncol(design_prelim))
    tt_prelim <- edgeR::topTags(lrt_prelim, n = nrow(dge_filtered), sort.by = "none")$table
    tt_prelim$PValue
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

# ---------------------------------------------------------------------------------------------------------------
# KIT-129 W4: the RUV designs beside today's grid, as KIT-124's simulator S3 ran them
# (docs/qa/workflow-execution/kit/kit-124/sim/s3_trx.R). The factors are estimated on upper-quartile normalised
# counts (the RUVSeq vignette's usage) and enter the edgeR model; the test reads the filtered counts. k is an
# argument: the rule that gives it is shared by every omic (KIT-129 W6) and is not computed here. Nothing below is
# called by today's path.
# ---------------------------------------------------------------------------------------------------------------

#' The declared RUV designs (KIT-129 section 1c): the recommended one (arm 7), a supplied control list (arm 8), the
#' first-pass nulls (arm 5) and no RUV (arm 1). Today's grid is not one of them: it is the default path.
#' @noRd
TRX_RUV_DESIGNS <- c("groups_all_features", "supplied_controls", "first_pass_nulls", "none")

#' The share of the first pass's most significant genes the first-pass control rule leaves out (arm 5).
#' @noRd
TRX_FIRST_PASS_DROPPED_SHARE <- 0.25

#' The fewest supplied control ids that must be genes of the filtered counts (KIT-129 audit F-06). The ticket states no
#' share; the rule is the conservative one: every unmatched id is warned of, and a list that leaves fewer than two
#' controls is refused (one control cannot carry a factor: k must be below the number of controls).
#' @noRd
TRX_MIN_SUPPLIED_CONTROLS <- 2L

#' A supplied control list in its canonical form (KIT-129 audit F-06, F-07): character, each id once, sorted in the C
#' locale (`method = "radix"`), so the list is the same on every machine and in whatever order it was given. An NA is
#' kept (last), for the controls' kernel to refuse.
#' @noRd
trxCanonicalControlIds <- function(control_ids) sort(unique(as.character(control_ids)), method = "radix", na.last = TRUE)

#' The record's `ruv_control_list_sha256` (KIT-129 section 1e, "hash of the sorted ids"): the sha256 of the canonical
#' list's ids as UTF-8 text, joined by a line feed, without a final one. It covers the list as supplied, matched or not.
#' @noRd
trxControlListSha256 <- function(control_ids) {
    digest::digest(paste(enc2utf8(trxCanonicalControlIds(control_ids)), collapse = "\n"), algo = "sha256", serialize = FALSE)
}

#' Refuse RUV inputs that are not one set of runs and genes (KIT-129 audit F-05): the design's `Run` must be the count
#' columns, and the logCPM's genes and runs the counts', identical in order and content. RUVSeq takes the replicate
#' sets by position, so a design in another order gives factors of the wrong runs without an error.
#' @param log_cpm_tmm The TMM logCPM; NULL where the caller has none.
#' @noRd
trxRuvCheckAlignment <- function(filtered_counts, design, log_cpm_tmm = NULL) {
    refuse <- function(message) trxAbort(message, c("trx_ruv_inputs_misaligned", "trx_error"))
    if (!identical(colnames(filtered_counts), as.character(design$Run))) {
        refuse("the design's runs are not the count columns in the same order")
    }
    if (is.null(log_cpm_tmm)) return(invisible(NULL))
    if (!identical(colnames(log_cpm_tmm), colnames(filtered_counts))) refuse("the logCPM's runs are not the count columns in the same order")
    if (!identical(rownames(log_cpm_tmm), rownames(filtered_counts))) refuse("the logCPM's genes are not the counts' genes in the same order")
    invisible(NULL)
}

#' The controls and the RUVSeq method of a declared design.
#' @param filtered_counts The filtered counts (genes by runs).
#' @param design The validated design (`Run`, `group`).
#' @param ruv_design One of `TRX_RUV_DESIGNS`.
#' @param control_ids The supplied control gene ids, for `supplied_controls`.
#' @return A list: `control_genes`, `control_rule` (`all_features`, `supplied_list`, `first_pass_drop_top_25` or NA),
#'   `method` (`RUVs`, `RUVg` or NA), `unmatched` (supplied ids that are not in the matrix) and `list_sha256`
#'   (`trxControlListSha256` of the supplied list; NA for the other designs). The controls of every rule are in the
#'   matrix's row order, so the order a list was given in, or the first pass's ranking, changes nothing. A supplied
#'   list is taken once per id (`trxCanonicalControlIds`). Refuses `trx_ruv_inputs_misaligned` (the design's runs are
#'   not the count columns) and `trx_ruv_controls_invalid` (a supplied id that is NA or empty; fewer than
#'   `TRX_MIN_SUPPLIED_CONTROLS` supplied ids in the matrix). Warns `trx_ruv_controls_unmatched` whenever a supplied
#'   id is not in the matrix, with their number and the first of them.
#' @noRd
trxRuvDesignControls <- function(filtered_counts, design, ruv_design, control_ids = character()) {
    trxRuvCheckAlignment(filtered_counts, design)
    genes <- rownames(filtered_counts)
    controls <- function(control_genes, control_rule, method, unmatched = character(), list_sha256 = NA_character_) {
        list(control_genes = control_genes, control_rule = control_rule, method = method, unmatched = unmatched, list_sha256 = list_sha256)
    }
    switch(ruv_design,
        groups_all_features = controls(genes, "all_features", "RUVs"),
        supplied_controls = {
            refuse <- function(message) trxAbort(message, c("trx_ruv_controls_invalid", "trx_error"))
            ids <- as.character(control_ids)
            blank <- is.na(ids) | !nzchar(trimws(ids))
            if (any(blank)) refuse(sprintf("%d of the %d supplied control ids are NA or empty", sum(blank), length(ids)))
            supplied <- trxCanonicalControlIds(ids)
            matched <- genes[genes %in% supplied]
            unmatched <- setdiff(supplied, genes)
            if (length(matched) < TRX_MIN_SUPPLIED_CONTROLS) {
                refuse(sprintf("%d of the %d supplied control ids are genes of the filtered counts; at least %d are needed", length(matched),
                    length(supplied), TRX_MIN_SUPPLIED_CONTROLS))
            }
            if (length(unmatched)) {
                rlang::warn(sprintf("%d of the %d supplied control ids are not genes of the filtered counts and are not used: %s%s",
                    length(unmatched), length(supplied), paste(utils::head(unmatched, 10L), collapse = ", "),
                    if (length(unmatched) > 10L) ", ..." else ""), class = c("trx_ruv_controls_unmatched", "trx_warning"))
            }
            controls(matched, "supplied_list", "RUVg", unmatched, trxControlListSha256(supplied))
        },
        first_pass_nulls = {
            requireOptionalPackage("edgeR", "the first pass of the RUV controls")
            first_p <- trxPrelimPValues(trxFilteredDge(filtered_counts, design$group), stats::model.matrix(~ group, data = design))
            dropped <- seq_len(ceiling(TRX_FIRST_PASS_DROPPED_SHARE * length(genes)))
            controls(genes[sort(order(first_p)[-dropped])], "first_pass_drop_top_25", "RUVg")
        },
        none = controls(character(), NA_character_, NA_character_),
        trxAbort(sprintf("RUV design '%s' has no control rule", ruv_design), c("trx_ruv_design_unsupported", "trx_error")))
}

#' Upper-quartile normalised counts, the scale RUVSeq's factors are estimated on (the vignette's
#' `betweenLaneNormalization(which = "upper")`, as every RUV arm of the simulator but today's).
#' @param filtered_counts The filtered counts (genes by runs).
#' @return The normalised counts, rounded. Refuses `trx_ruv_upper_quartile_invalid` when a run's upper quartile is 0.
#' @noRd
trxUpperQuartileCounts <- function(filtered_counts) {
    requireOptionalPackage("EDASeq", "the upper-quartile normalisation of the RUV step")
    uq <- EDASeq::betweenLaneNormalization(filtered_counts, which = "upper", round = TRUE)
    if (!all(is.finite(uq))) {
        trxAbort("a run's upper-quartile count is 0, so its counts cannot be upper-quartile normalised",
            c("trx_ruv_upper_quartile_invalid", "trx_error"))
    }
    uq
}

#' The RUV factors of a declared design (KIT-129 W4): `RUVSeq::RUVs` with the groups as replicate sets, or
#' `RUVSeq::RUVg`, on upper-quartile counts with the given controls and k.
#' @param filtered_counts The filtered counts (genes by runs); the test's input, returned unchanged.
#' @param log_cpm_tmm Their TMM logCPM.
#' @param design The validated design (`Run`, `group`). Columns named `W_...` are an earlier RUV step's factors (the
#'   package's rule, `trxDaFormula`) and are dropped, so the returned `pdata` carries this call's factors alone.
#' @param control_genes The control genes (`trxRuvDesignControls`).
#' @param k The number of factors, a whole number; 0 gives no factor.
#' @param method `"RUVs"` or `"RUVg"`; not read when k is 0.
#' @return `trxRuvs`' list without `control_source`: `W` (`W_1..W_k`, a data frame by run, no column when k is 0),
#'   `ruvs_counts` (the filtered counts), `log_cpm_ruvs` (visualisation and export only: `removeBatchEffect` with W;
#'   the TMM logCPM itself when k is 0) and `pdata` (the design with W, the DA kernel's input). Refuses
#'   `trx_ruv_inputs_misaligned` (`trxRuvCheckAlignment`), `trx_ruv_k_invalid` (k is not one whole number of 0 or
#'   more), `trx_ruv_controls_invalid` (a control that is not a gene of the matrix) and `trx_ruv_k_unavailable` (k
#'   leaves no residual degree of freedom beside the groups, is not below the number of controls, the runs of every
#'   group are identical under RUVs, or RUVSeq returned fewer factors: k is never reduced). The first two are a k
#'   above the limit and also carry `ruv_k_above_limit` and `ruv_k_refused` (KIT-129 item 2).
#' @noRd
trxRuvFactors <- function(filtered_counts, log_cpm_tmm, design, control_genes, k, method = c("RUVs", "RUVg")) {
    validateTrxCounts(filtered_counts)
    trxRuvCheckAlignment(filtered_counts, design, log_cpm_tmm)
    if (!checkmate::test_count(k)) {
        trxAbort(sprintf("k must be one whole number of 0 or more, not %s", substr(deparse1(k), 1L, 40L)), c("trx_ruv_k_invalid", "trx_error"))
    }
    k <- as.integer(round(k))
    design <- design[, !startsWith(names(design), "W_"), drop = FALSE]
    pdata <- data.frame(design, row.names = design$Run)
    if (k == 0L) {
        return(list(W = pdata[, character(), drop = FALSE], ruvs_counts = filtered_counts, log_cpm_ruvs = log_cpm_tmm, pdata = pdata))
    }
    method <- match.arg(method)
    requireOptionalPackage("RUVSeq", "the transcriptomics RUV factors")
    refuse_k <- function(why, above_limit = FALSE) trxAbort(sprintf("k = %d cannot be estimated: %s", k, why),
        c("trx_ruv_k_unavailable", "trx_error", if (above_limit) c("ruv_k_above_limit", "ruv_k_refused")))
    if (!all(control_genes %in% rownames(filtered_counts))) {
        trxAbort(sprintf("%d control gene(s) are not in the filtered counts", sum(!control_genes %in% rownames(filtered_counts))),
            c("trx_ruv_controls_invalid", "trx_error"))
    }
    residual_limit <- ncol(filtered_counts) - length(unique(as.character(design$group))) - 1L
    if (k > residual_limit) {
        refuse_k(sprintf("%d runs in %d groups leave a residual degree of freedom for at most %d factor(s)", ncol(filtered_counts),
            length(unique(as.character(design$group))), max(residual_limit, 0L)), above_limit = TRUE)
    }
    if (k >= length(control_genes)) {
        refuse_k(sprintf("there are %d control gene(s) and k must be below their number", length(control_genes)), above_limit = TRUE)
    }
    uq <- trxUpperQuartileCounts(filtered_counts)
    # RUVs estimates the factors from the differences of the runs within a group; when there is none it fails on an
    # empty index (KIT-129 audit F-08).
    runs_by_group <- split(seq_len(ncol(uq)), as.character(design$group))
    if (method == "RUVs" && all(vapply(runs_by_group, \(runs) all(uq[, runs] == uq[, runs[[1L]]]), logical(1)))) {
        refuse_k("the runs of every group are identical, so their differences hold no unwanted variation")
    }
    fit <- if (method == "RUVs") {
        RUVSeq::RUVs(uq, control_genes, k, RUVSeq::makeGroups(design$group))
    } else {
        RUVSeq::RUVg(uq, control_genes, k)
    }
    # RUVSeq drops the factors whose singular value is below its tolerance without saying so.
    if (ncol(fit$W) != k) refuse_k(sprintf("RUVSeq::%s returned %d factor(s)", method, ncol(fit$W)))
    w_matrix <- data.frame(fit$W, row.names = design$Run)
    names(w_matrix) <- paste0("W_", seq_len(k))
    pdata <- cbind(pdata, w_matrix)
    log_cpm_corrected <- limma::removeBatchEffect(log_cpm_tmm, covariates = w_matrix, design = stats::model.matrix(~ group, data = pdata))
    list(W = w_matrix, ruvs_counts = filtered_counts, log_cpm_ruvs = log_cpm_corrected, pdata = pdata)
}
