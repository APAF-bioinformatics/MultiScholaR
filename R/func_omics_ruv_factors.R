# KIT-127: lipidomics and metabolomics RUV with the factors in the model (KIT-124's gate, 2026-11-18). One kernel for
# both omics. Per assay it computes the RUV-III factors W from the features observed in every sample, with the groups
# as replicate sets and either every such feature or a supplied list (internal standards) as controls: arms 7 and 8
# of KIT-124's simulator S4. The differential test then reads the uncorrected values with W in the limma model
# (`runLipidsDA` and `runMetabolitesDA`, `ruv_factors`). The corrected matrix is still produced, for figures, export
# and integration, with each feature's level kept (KIT-123 R2's `restoreRuvFeatureLevels`). Today's design (the
# sweep's controls, `ruvIII_C_Varying`, the test on the corrected matrix) is not routed through this file and is
# unchanged.

#' The designs of the factors-in-the-model path, by the arm of KIT-124's simulator S4 each one is.
#' @noRd
omicsRuvFactorsDesigns <- c(groups_all_features = "arm 7", supplied_controls = "arm 8")

#' Stop for a wrong argument of the factors-in-the-model kernel.
#' @noRd
abortOmicsRuvFactors <- function(message) {
    rlang::abort(message, class = "omics_ruv_factors_invalid", call = NULL)
}

#' The names of the assay list and feature id slots of a lipidomics or metabolomics object.
#' @noRd
omicsRuvAssaySlots <- function(theObject) {
    if (inherits(theObject, "LipidomicsAssayData")) return(list(data = "lipid_data", id = "lipid_id_column"))
    if (inherits(theObject, "MetaboliteAssayData")) return(list(data = "metabolite_data", id = "metabolite_id_column"))
    abortOmicsRuvFactors("theObject must be a LipidomicsAssayData or MetaboliteAssayData object")
}

#' One assay as a feature-by-sample matrix on the design's samples, in the design's order; non-finite values missing.
#' @noRd
omicsRuvAssayValues <- function(assay, id_column, design_samples, assay_name) {
    samples <- design_samples[design_samples %in% colnames(assay)]
    ids <- as.character(assay[[id_column]])
    if (length(samples) < 2L || !length(ids) || anyNA(ids) || anyDuplicated(ids)) {
        abortOmicsRuvFactors(sprintf("assay '%s' needs unique feature ids and at least two samples of the design", assay_name))
    }
    values <- as.matrix(as.data.frame(lapply(assay[samples], as.numeric), check.names = FALSE))
    values[!is.finite(values)] <- NA_real_
    dimnames(values) <- list(ids, samples)
    values
}

#' The number of factors per assay (`ruv_k`): one whole number of 0 or more for every assay, or one per assay by name. `NA` is a
#' k that could not be estimated; the assay is then refused, not guessed. An infinite k or `NaN` is a wrong argument.
#' @noRd
resolveOmicsRuvFactorsK <- function(ruv_k, assay_names) {
    per_assay <- if (is.list(ruv_k)) ruv_k else stats::setNames(rep(list(ruv_k), length(assay_names)), assay_names)
    if (!all(assay_names %in% names(per_assay))) {
        abortOmicsRuvFactors(sprintf("ruv_k names no value for assay(s): %s", paste(setdiff(assay_names, names(per_assay)), collapse = ", ")))
    }
    lapply(per_assay[assay_names], function(k) {
        if (length(k) == 1L && is.na(k) && !(is.numeric(k) && is.nan(k))) return(NA_integer_)
        if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 0 || k != round(k)) {
            abortOmicsRuvFactors("ruv_k must be a whole number of 0 or more (one for every assay, or one per assay by name)")
        }
        as.integer(k)
    })
}

#' Today's estimate of k, bounded (KIT-127 D4): the automatic sweep's `best_k` per assay, limited to the sample-size
#' limit `calculateAdaptiveMaxK` (in the sweep that limit only weights the score, so a larger k can be returned).
#' @param theObject A `LipidomicsAssayData` or `MetaboliteAssayData` object.
#' @param ruv_results The per-assay list of `runLipidPerAssayRuvOptimization` or `runPerAssayRuvOptimization`.
#' @return `list(k, ruv_k_rule = "sweep_bounded", sweep_k, limit)`, each per assay; `k` is `NA` where the sweep failed.
#' @noRd
boundOmicsRuvSweepK <- function(theObject, ruv_results) {
    slots <- omicsRuvAssaySlots(theObject)
    assays <- methods::slot(theObject, slots$data)
    if (!is.list(ruv_results) || !all(names(assays) %in% names(ruv_results))) {
        abortOmicsRuvFactors("ruv_results must hold the sweep's result of every assay")
    }
    design_samples <- as.character(theObject@design_matrix[[theObject@sample_id]])
    sweep_k <- lapply(ruv_results[names(assays)], function(result) {
        if (isTRUE(result$success) && length(result$best_k) == 1L && !is.na(result$best_k)) as.integer(result$best_k) else NA_integer_
    })
    limit <- lapply(assays, function(assay) calculateAdaptiveMaxK(length(intersect(colnames(assay), design_samples))))
    list(k = Map(function(k, cap) if (is.na(k)) NA_integer_ else min(k, cap), sweep_k, limit), ruv_k_rule = "sweep_bounded",
        sweep_k = sweep_k, limit = limit)
}

#' The RUV-III factors of one assay (KIT-124 S4, `ruvW`): W from the features observed in every sample.
#' @param values Feature-by-sample matrix of the assay.
#' @param replicate_matrix Sample-by-replicate-set matrix M, rows in the order of `values`' columns.
#' @param k Number of factors, or `NA` when it could not be estimated.
#' @param control_ids The supplied controls, or NULL for every feature observed in every sample.
#' @param units The biological unit of each sample on a design with declared technical replicates, or NULL.
#' @param groups The number of groups of the test (the object's group column). With the groups as the replicate sets it
#'   is the number of replicate sets; on the technical-replicate design (the biological units as the replicate sets,
#'   arms 15 and 16 of simulator S4) it is smaller, and the limits count the groups across the units (decision M28).
#' @return `list(status, reason, reason_class, W, k, control_ids, complete_features)`; `status` is `ok`, `none` (k = 0)
#'   or `refused`; `reason_class` names a refusal by the class the other omics raise for it (KIT-129 item 2):
#'   `ruv_k_above_limit` or `ruv_k_not_estimated`.
#' @noRd
computeOmicsRuvAssayFactors <- function(values, replicate_matrix, k, control_ids = NULL, units = NULL,
                                        groups = ncol(replicate_matrix)) {
    complete <- rowSums(is.na(values)) == 0L
    y <- t(values[complete, , drop = FALSE])
    ctl <- if (is.null(control_ids)) rep(TRUE, ncol(y)) else colnames(y) %in% control_ids
    result <- function(status, reason = NA_character_, w = NULL, k_used = 0L, reason_class = if (identical(status, "refused")) "ruv_k_above_limit" else NA_character_) {
        list(status = status, reason = reason, reason_class = reason_class, W = w, k = k_used, control_ids = colnames(y)[ctl],
            complete_features = ncol(y))
    }
    if (is.na(k)) return(result("refused", "k could not be estimated", reason_class = "ruv_k_not_estimated"))
    if (k == 0L) return(result("none"))
    if (k > sum(ctl)) {
        return(result("refused", sprintf("k = %d is above the %d controls observed in every sample", k, sum(ctl))))
    }
    if (k > nrow(y) - groups - 1L) {
        return(result("refused", sprintf("k = %d leaves no residual degree of freedom (%d samples in %d groups)",
            k, nrow(y), groups)))
    }
    # The rank RUV-III's replicate sets leave (binding only when the replicate sets are not the groups).
    if (k > nrow(y) - ncol(replicate_matrix)) {
        return(result("refused", sprintf("k = %d is above the %d factors RUV-III can estimate from %d samples in %d replicate sets",
            k, nrow(y) - ncol(replicate_matrix), nrow(y), ncol(replicate_matrix))))
    }
    # Decision M28: with declared technical replicates the biological units, not the injections, limit k.
    if (!is.null(units)) {
        limit <- ruvFactorCountLimit(nrow(y), ncol(replicate_matrix), groups, sum(ctl), units = length(unique(units)))
        if (k > limit$k_max_units) {
            return(result("refused", sprintf(paste("k = %d is above the limit of %d: the design has technical replicates, so its",
                "%d biological units, not its %d injections, carry the test and need a residual degree of freedom beside the",
                "%d groups"), k, limit$k_max_units, limit$units, nrow(y), groups)))
        }
    }
    # W as `computeRuvCovariates` forms it (KIT-123 R2): the controls centred on their means, times ac' (ac ac')^-1.
    ac <- ruv::RUVIII(Y = y, M = replicate_matrix, ctl = ctl, k = k, return.info = TRUE)$fullalpha[seq_len(k), ctl, drop = FALSE]
    w <- sweep(y[, ctl, drop = FALSE], 2L, colMeans(y[, ctl, drop = FALSE])) %*% t(ac) %*% solve(ac %*% t(ac))
    w <- sweep(w, 2L, colMeans(w))
    dimnames(w) <- list(rownames(y), paste0("W", seq_len(k)))
    result("ok", w = w, k_used = k)
}

#' One corrected assay with each feature's level restored (KIT-123 R2), against the uncorrected values.
#' @noRd
restoreOmicsRuvAssayLevels <- function(corrected_assay, original_values, id_column) {
    samples <- colnames(original_values)
    ids <- as.character(corrected_assay[[id_column]])
    corrected_values <- as.matrix(as.data.frame(corrected_assay[samples], check.names = FALSE))
    rownames(corrected_values) <- ids
    restored <- restoreRuvFeatureLevels(corrected_values, original_values[ids, samples, drop = FALSE])
    corrected_assay[samples] <- as.data.frame(restored, check.names = FALSE)
    corrected_assay
}

#' RUV with the factors in the model, for lipidomics and metabolomics (KIT-127)
#'
#' Per assay: the RUV-III factors W from the features observed in every sample, with the levels of
#' `ruv_grouping_variable` as replicate sets and the design's controls (KIT-124 S4, arms 7 and 8); and the corrected
#' matrix, `ruvIII_C_Varying` with the same replicate sets, k and controls, each feature's level restored. The
#' differential test reads the input object, not the corrected one, with `factors` (`runLipidsDA`, `runMetabolitesDA`):
#' fold changes and p-values come from the uncorrected values with W in the model, and the corrected matrix is for
#' display and export (intensities, heatmaps, PCA, the integration input; KIT-124 design review, M10). The two are
#' different corrections of the same unwanted variation, RUV-III-C's own against the model's fit of W per feature, so
#' a group difference read from the corrected matrix need not equal the model's fold change (up to 0.18 log2 on the
#' audit's data for the features observed in every sample).
#' An assay whose factors cannot be computed gets no RUV and a warning of class `omics_ruv_factors_refused` (today's
#' behaviour, kept: the other omics' new designs stop instead); its record's `reason_class` and the warning name the
#' refusal by the class the other omics raise (`ruv_k_above_limit`, `ruv_k_not_estimated`; KIT-129 item 2).
#'
#' @param theObject A `LipidomicsAssayData` or `MetaboliteAssayData` object (log2 values, missing values kept).
#' @param ruv_design `"groups_all_features"`: every feature observed in every sample is a control.
#'   `"supplied_controls"`: the features of `ruv_control_ids` are the controls, and are not tested (the names of
#'   KIT-129 section 1e).
#' @param ruv_k The number of factors: a whole number of 0 or more, or one per assay by name. No default under
#'   `"user"` and `"sweep_bounded"`; not given (or NULL) under `"num_sv_be"`. 0 is no RUV; `NA` (a k that could not be
#'   estimated) refuses the assay.
#' @param ruv_k_rule Where k comes from, recorded with it: `"num_sv_be"` (the estimated count per assay, by the shared
#'   resolver: `estimateRuvFactorCount` with 200 permutations on the assay's features observed in every sample, the
#'   groups removed, or the biological units on a design with declared technical replicates; KIT-129 item 3, decisions
#'   M28, M30), `"user"` or `"sweep_bounded"` (`boundOmicsRuvSweepK`). `"sweep_bounded"` is refused with
#'   `"groups_all_features"`: that k was not measured with all features as controls. An assay whose count cannot be
#'   estimated is refused as any other (no RUV, the warning), with the count's reason. This kernel keeps its own limit
#'   on k and its refusal mode (an assay refused, not the run): the resolver is given no limit.
#' @param ruv_k_seed The estimated count's seed (recorded with it).
#' @param ruv_control_ids For `"supplied_controls"`: the control feature ids per assay, by assay name. A
#'   supplied control with a missing value is used neither for W nor for the correction; it is named in a warning of
#'   class `omics_ruv_controls_incomplete` and in `controls_left_out`, and is still not tested.
#' @param ruv_grouping_variable The design column whose levels are the replicate sets; the object's group column
#'   when NULL. The declared technical-replicate column gives the technical-replicate design (the biological units as
#'   the replicate sets, arms 15 and 16 of simulator S4): k is then limited by the units and the groups across them
#'   (decision M28), not by the units counted as groups.
#' @return `list(corrected, factors, parameters)`: `corrected` the object with the corrected assays, its `args`
#'   stamped with `ruv_design`, `ruv_grouping_variable`, `ruv_number_k` (the key `ruvIII_C_Varying` writes) and `ctrl` when an assay was corrected (the DA
#'   refuses a stamped object together with factors, class `omics_ruv_object_already_corrected`); `factors` per
#'   assay `list(assay, status, reason, reason_class, W, k, k_requested, ruv_k_rule, ruv_k_estimate, k_estimate, control_ids, control_count, controls_left_out,
#'   complete_features, not_tested, ruv_design, ruv_grouping_variable)`, `W` a sample-by-k matrix or NULL,
#'   `k_estimate` the count's record for `plotRuvFactorCount` (NULL under the other rules);
#'   `parameters` the arguments.
#' @noRd
runOmicsRuvFactors <- function(theObject, ruv_design, ruv_k, ruv_k_rule, ruv_control_ids = NULL,
                               ruv_grouping_variable = NULL, ruv_k_seed = RUV_FACTOR_COUNT_SEED) {
    slots <- omicsRuvAssaySlots(theObject)
    if (missing(ruv_design) || !is.character(ruv_design) || length(ruv_design) != 1L || !ruv_design %in% names(omicsRuvFactorsDesigns)) {
        abortOmicsRuvFactors(sprintf("ruv_design must be one of: %s", paste(names(omicsRuvFactorsDesigns), collapse = ", ")))
    }
    if (missing(ruv_k_rule) || !is.character(ruv_k_rule) || length(ruv_k_rule) != 1L || !ruv_k_rule %in% RUV_K_RULES) {
        abortOmicsRuvFactors("ruv_k_rule must be \"num_sv_be\", \"user\" or \"sweep_bounded\"")
    }
    counted <- identical(ruv_k_rule, "num_sv_be")
    if (counted && !missing(ruv_k) && !is.null(ruv_k)) abortOmicsRuvFactors("ruv_k is not given with ruv_k_rule = \"num_sv_be\": the count gives k")
    if (!counted && missing(ruv_k)) abortOmicsRuvFactors("ruv_k has no default: give the number of factors")
    if (!checkmate::test_int(ruv_k_seed)) abortOmicsRuvFactors("ruv_k_seed must be a whole number")
    assays <- methods::slot(theObject, slots$data)
    id_column <- methods::slot(theObject, slots$id)
    design <- as.data.frame(theObject@design_matrix)
    sample_id <- theObject@sample_id
    grouping <- ruv_grouping_variable %||% theObject@group_id
    if (!is.character(grouping) || length(grouping) != 1L || !grouping %in% colnames(design)) {
        abortOmicsRuvFactors("ruv_grouping_variable must name a column of the design")
    }
    if (anyNA(design[[grouping]])) {
        abortOmicsRuvFactors(sprintf("the design column '%s' has missing values: every sample needs a replicate set", grouping))
    }
    if (!length(assays) || is.null(names(assays)) || !all(nzchar(names(assays)))) abortOmicsRuvFactors("the object needs named assays")
    supplied <- identical(ruv_design, "supplied_controls")
    if (!supplied && !is.null(ruv_control_ids)) abortOmicsRuvFactors("ruv_control_ids is for ruv_design = \"supplied_controls\" only")
    if (supplied && !(is.list(ruv_control_ids) && all(names(assays) %in% names(ruv_control_ids)))) {
        abortOmicsRuvFactors("ruv_design = \"supplied_controls\" needs ruv_control_ids for every assay, by assay name")
    }
    # Audit F1: the bounded sweep's k was never measured with all features as controls (KIT-129 brings the estimate).
    if (!supplied && identical(ruv_k_rule, "sweep_bounded")) {
        abortOmicsRuvFactors("ruv_k_rule = \"sweep_bounded\" is refused with ruv_design = \"groups_all_features\": that k was not validated for it")
    }
    k_requested <- if (counted) stats::setNames(rep(list(NA_integer_), length(assays)), names(assays)) else
        resolveOmicsRuvFactorsK(ruv_k, names(assays))
    design_samples <- as.character(design[[sample_id]])

    values <- list()
    factors <- list()
    for (assay_name in names(assays)) {
        values[[assay_name]] <- omicsRuvAssayValues(assays[[assay_name]], id_column, design_samples, assay_name)
        control_ids <- if (supplied) as.character(ruv_control_ids[[assay_name]])
        absent <- setdiff(control_ids, rownames(values[[assay_name]]))
        if (supplied && (!length(control_ids) || length(absent))) {
            abortOmicsRuvFactors(sprintf("assay '%s': %s", assay_name, if (!length(control_ids)) "no control ids supplied" else
                sprintf("%d supplied control(s) are not in the assay: %s", length(absent), paste(utils::head(absent, 3L), collapse = ", "))))
        }
        samples <- colnames(values[[assay_name]])
        assay_design <- design[match(samples, design_samples), , drop = FALSE]
        replicate_matrix <- getRuvIIIReplicateMatrixHelper(assay_design, !!rlang::sym(sample_id), !!rlang::sym(grouping))
        technical <- resolveProtTechnicalReplicateBlock(assay_design, theObject@group_id, theObject@technical_replicate_id)
        units <- if (isTRUE(technical$technical)) technical$block
        # The groups of the test, counted apart from the replicate sets (the technical-replicate design: the units).
        groups <- if (isTRUE(theObject@group_id %in% colnames(assay_design))) {
            length(unique(as.character(assay_design[[theObject@group_id]])))
        } else ncol(replicate_matrix)
        # KIT-129 item 3: the estimated count per assay by the shared resolver (no limit there: this kernel's own limit
        # refuses the assay below). An assay the count cannot be made on is refused with the count's reason.
        counted_k <- if (counted) tryCatch(resolveRuvFactorCount("num_sv_be", if (supplied) "supplied_list" else "all_features",
            list(k_max = .Machine$integer.max), estimate_count = \() estimateRuvFactorCount(values[[assay_name]],
                assay_design[[grouping]], seed = ruv_k_seed, technical_replicates = isTRUE(technical$technical), units = units)),
            ruv_factor_count_unavailable = \(e) e)
        k_assay <- if (!counted) k_requested[[assay_name]] else if (inherits(counted_k, "condition")) NA_integer_ else counted_k$k
        record <- computeOmicsRuvAssayFactors(values[[assay_name]], replicate_matrix[samples, , drop = FALSE],
            k_assay, control_ids, units = units, groups = groups)
        if (inherits(counted_k, "condition")) {
            record$reason <- conditionMessage(counted_k)
            record$reason_class <- class(counted_k)[[1L]]
        }
        if (counted && identical(record$status, "none")) record$reason <- counted_k$reason
        factors[[assay_name]] <- list(assay = assay_name, status = record$status, reason = record$reason,
            reason_class = record$reason_class, W = record$W, k = record$k,
            k_requested = k_requested[[assay_name]], ruv_k_rule = ruv_k_rule,
            ruv_k_estimate = if (counted && !inherits(counted_k, "condition")) counted_k$k_estimate else NA_integer_,
            k_estimate = if (counted && !inherits(counted_k, "condition")) counted_k$estimate, control_ids = record$control_ids,
            control_count = length(record$control_ids), controls_left_out = setdiff(control_ids, record$control_ids),
            complete_features = record$complete_features,
            not_tested = if (supplied) control_ids else character(), ruv_design = ruv_design, ruv_grouping_variable = grouping)
    }
    # Audit F9: a supplied control with a missing value is no control, for W and for the correction alike.
    left_out <- Filter(length, lapply(factors, `[[`, "controls_left_out"))
    if (length(left_out)) {
        rlang::warn(c("supplied controls with missing values are not used as controls (they are still not tested)",
            stats::setNames(sprintf("%s: %s", names(left_out), vapply(left_out, paste, "", collapse = ", ")), rep("x", length(left_out)))),
            class = "omics_ruv_controls_incomplete")
    }
    refused <- Filter(function(record) identical(record$status, "refused"), factors)
    if (length(refused)) {
        rlang::warn(c("RUV factors could not be computed for some assays; they are analysed without RUV",
            stats::setNames(sprintf("%s: %s (%s)", names(refused), vapply(refused, `[[`, "", "reason"),
                vapply(refused, `[[`, "", "reason_class")), rep("x", length(refused)))),
            class = "omics_ruv_factors_refused")
    }

    # The corrected matrix: today's method on the assays that have factors, with the same replicate sets, k and controls.
    corrected <- theObject
    to_correct <- names(Filter(function(record) identical(record$status, "ok"), factors))
    if (length(to_correct)) {
        subset_object <- theObject
        methods::slot(subset_object, slots$data) <- assays[to_correct]
        correction_controls <- lapply(stats::setNames(to_correct, to_correct), function(assay_name) {
            if (supplied) factors[[assay_name]]$control_ids else rownames(values[[assay_name]])
        })
        result <- ruvIII_C_Varying(subset_object, ruv_grouping_variable = grouping, ruv_number_k = lapply(factors[to_correct], `[[`, "k"),
            ctrl = correction_controls)
        result_assays <- methods::slot(result, slots$data)
        for (assay_name in to_correct) {
            corrected_assay <- result_assays[[assay_name]]
            if (is.null(corrected_assay) || identical(corrected_assay, assays[[assay_name]]) ||
                !all(colnames(values[[assay_name]]) %in% colnames(corrected_assay))) {
                rlang::abort(sprintf("the RUV-III correction returned no corrected values for assay '%s' or for one of its samples", assay_name),
                    class = "omics_ruv_correction_failed", call = NULL)
            }
            assays[[assay_name]] <- restoreOmicsRuvAssayLevels(corrected_assay, values[[assay_name]], id_column)
        }
        methods::slot(corrected, slots$data) <- assays
        # Audit F3: the corrected object records its correction, under the keys `ruvIII_C_Varying` writes, and the design.
        corrected@args$ruv_design <- ruv_design
        corrected@args$ruv_grouping_variable <- grouping
        corrected@args$ruv_number_k <- lapply(factors, `[[`, "k")
        corrected@args$ctrl <- correction_controls
    }
    list(corrected = corrected, factors = factors,
        parameters = list(ruv_design = ruv_design, ruv_k = k_requested, ruv_k_rule = ruv_k_rule,
            ruv_k_seed = if (counted) as.integer(ruv_k_seed) else NA_integer_,
            ruv_control_ids = ruv_control_ids, ruv_grouping_variable = grouping))
}

#' Stop when an object this kernel has corrected is given to the DA together with factors: W in the model of
#' corrected values removes the unwanted variation twice (audit F3).
#' @noRd
refuseOmicsRuvCorrectedObject <- function(theObject) {
    if (is.null(theObject@args$ruv_design)) return(invisible(NULL))
    rlang::abort(paste("theObject is the corrected object of runOmicsRuvFactors(); with ruv_factors the DA reads the",
        "uncorrected object"), class = c("omics_ruv_object_already_corrected", "omics_ruv_factors_invalid"), call = NULL)
}

#' The tested features with no value in one whole group, by feature and group (audit F5).
#' @return A data frame with the columns `feature` and `group`, in the matrix's row order; no row when none.
#' @noRd
findOmicsRuvUnobservedGroups <- function(expr_matrix, design_matrix, sample_id_col, group_col) {
    samples <- as.character(design_matrix[[sample_id_col]])
    common <- intersect(colnames(expr_matrix), samples)
    groups <- if (is.null(group_col)) character() else as.character(design_matrix[[group_col]])[match(common, samples)]
    levels <- unique(groups[!is.na(groups)])
    if (!nrow(expr_matrix) || !length(levels)) return(data.frame(feature = character(), group = character(), stringsAsFactors = FALSE))
    seen <- is.finite(expr_matrix[, common, drop = FALSE])
    at <- which(do.call(cbind, lapply(levels, function(level) rowSums(seen[, groups %in% level, drop = FALSE]) == 0L)), arr.ind = TRUE)
    at <- at[order(at[, 1L], at[, 2L]), , drop = FALSE]
    data.frame(feature = rownames(expr_matrix)[at[, 1L]], group = levels[at[, 2L]], stringsAsFactors = FALSE)
}

#' One assay's DA input with the RUV factors in the model (KIT-124 S4, `testWithW` and `daKernel`)
#'
#' The supplied controls are left out of the test; the W columns are added to the design and to the formula. A
#' feature is left out only when its observed samples cannot estimate a W coefficient (the fit would then compare
#' its groups without the adjustment); those features are counted in a warning of class
#' `omics_ruv_features_not_tested` and named in `not_tested`. A feature with no value in one whole group is kept:
#' limma leaves that group's coefficient missing for it and tests it in the comparisons between the groups where it
#' is observed, and it is listed in `not_observed` (the user's decision of 2026-11-18 on audit F5; KIT-124 design
#' review, M9). This departs from the kernel of KIT-124's simulator, which left out every feature that could not
#' estimate all coefficients of the model; the two rules agree where every feature is observed in every group, which
#' is the data the simulator's arms 7 and 8 were scored on. With W in the model the DA forms each comparison's standard
#' error exactly per feature (`runOmicsRuvFactorsDaKernel`, decision M29); limma's contrast step alone would form it
#' from the full design's correlations, approximately for a feature with missing values. An assay without factors (k = 0 or refused) keeps the design and the formula. A model that leaves no residual degree
#' of freedom on the assay's samples stops with class `omics_ruv_no_residual_df`.
#' @param expr_matrix Feature-by-sample matrix of uncorrected values.
#' @param design_matrix,formula_string,sample_id_col As the DA kernel reads them.
#' @param group_col The design's group column, for `not_observed`; NULL lists nothing.
#' @param assay_factors The assay's entry of `runOmicsRuvFactors()$factors`.
#' @param assay_name The assay's name; the factors of another assay are refused.
#' @return `list(expr_matrix, design_matrix, formula_string, not_estimable, not_tested, not_observed)`: `not_tested`
#'   is `list(supplied_controls, not_estimable)`, feature ids; `not_observed` a data frame of `feature` and `group`.
#' @noRd
applyOmicsRuvFactorsToDaInput <- function(expr_matrix, design_matrix, formula_string, sample_id_col, assay_factors,
                                          assay_name, group_col = NULL, technical_replicate_id = NULL) {
    if (!is.list(assay_factors) || !all(c("assay", "status", "W", "not_tested") %in% names(assay_factors))) {
        abortOmicsRuvFactors(sprintf("ruv_factors holds no RUV factors for assay '%s'", assay_name))
    }
    if (!identical(assay_factors$assay, assay_name)) {
        abortOmicsRuvFactors(sprintf("the RUV factors of assay '%s' were given for assay '%s'", assay_factors$assay, assay_name))
    }
    supplied_controls <- intersect(rownames(expr_matrix), assay_factors$not_tested)
    expr_matrix <- expr_matrix[!rownames(expr_matrix) %in% assay_factors$not_tested, , drop = FALSE]
    w <- assay_factors$W
    if (is.null(w)) {
        return(list(expr_matrix = expr_matrix, design_matrix = design_matrix, formula_string = formula_string,
            not_estimable = character(), not_tested = list(supplied_controls = supplied_controls, not_estimable = character()),
            not_observed = findOmicsRuvUnobservedGroups(expr_matrix, design_matrix, sample_id_col, group_col)))
    }
    samples <- as.character(design_matrix[[sample_id_col]])
    absent <- setdiff(intersect(colnames(expr_matrix), samples), rownames(w))
    if (length(absent)) {
        abortOmicsRuvFactors(sprintf("assay '%s': the RUV factors have no row for sample(s): %s", assay_name,
            paste(utils::head(absent, 3L), collapse = ", ")))
    }
    if (any(colnames(w) %in% colnames(design_matrix))) {
        abortOmicsRuvFactors(sprintf("the design already has a column named %s", paste(intersect(colnames(w), colnames(design_matrix)), collapse = ", ")))
    }
    for (name in colnames(w)) design_matrix[[name]] <- unname(w[match(samples, rownames(w)), name])
    formula_string <- paste(formula_string, "+", paste(colnames(w), collapse = " + "))
    not_estimable <- character()
    model_design <- as.data.frame(design_matrix)
    rownames(model_design) <- samples
    model_formula <- stats::as.formula(formula_string)
    x <- stats::model.matrix(model_formula, stats::model.frame(model_formula, model_design))
    common <- intersect(colnames(expr_matrix), rownames(x))
    x <- x[common, , drop = FALSE]
    # Audit F2: the kernel's bound on k counts the replicate sets; the fitted model can hold more (covariates).
    residual_df <- nrow(x) - qr(x)$rank
    if (residual_df < 1L) {
        rlang::abort(sprintf(paste("assay '%s': the model '%s' with the RUV factors leaves no residual degree of freedom",
            "(%d samples, rank %d); use fewer factors or a smaller model"), assay_name, formula_string, nrow(x), nrow(x) - residual_df),
            class = "omics_ruv_no_residual_df", call = NULL)
    }
    # Decision M28: with declared technical replicates the biological units carry the test (the block), so the model needs
    # a residual degree of freedom between the units too.
    technical <- if (!is.null(group_col) && !is.null(technical_replicate_id))
        resolveProtTechnicalReplicateBlock(model_design[common, , drop = FALSE], group_col, technical_replicate_id)
    unit_df <- if (isTRUE(technical$technical)) ruvUnitResidualDf(x, technical$block) else NA_integer_
    if (!is.na(unit_df) && unit_df < 1L) {
        rlang::abort(sprintf(paste("assay '%s': the model '%s' with the RUV factors leaves no residual degree of freedom between",
            "the biological units: the design has technical replicates, so its %d units, not its %d injections, carry the test;",
            "use fewer factors (at most %d with this formula) or a smaller model"), assay_name, formula_string,
            length(unique(technical$block)), nrow(x), ncol(w) + unit_df - 1L),
            class = "omics_ruv_no_residual_df", call = NULL)
    }
    if (anyNA(expr_matrix)) {
        seen <- is.finite(expr_matrix[, common, drop = FALSE])
        pattern <- apply(seen, 1L, function(row) paste(which(row), collapse = ","))
        # A W coefficient is estimable when its column is not in the span of the model's other columns on the
        # feature's observed samples (audit F5: a group without values does not leave the feature out).
        estimable <- vapply(split(seq_len(nrow(seen)), pattern), function(rows) {
            x_seen <- x[seen[rows[[1L]], ], , drop = FALSE]
            rank_all <- qr(x_seen)$rank
            all(vapply(colnames(w), function(name) qr(x_seen[, colnames(x_seen) != name, drop = FALSE])$rank < rank_all, logical(1L)))
        }, logical(1L))
        keep <- unname(estimable[pattern])
        not_estimable <- rownames(expr_matrix)[!keep]
        expr_matrix <- expr_matrix[keep, , drop = FALSE]
        if (length(not_estimable)) {
            rlang::warn(sprintf(paste("assay '%s': %d feature(s) have too few values to estimate the RUV factors' coefficients",
                "and are not tested"), assay_name, length(not_estimable)), class = "omics_ruv_features_not_tested")
        }
    }
    list(expr_matrix = expr_matrix, design_matrix = design_matrix, formula_string = formula_string, not_estimable = not_estimable,
        not_tested = list(supplied_controls = supplied_controls, not_estimable = not_estimable),
        not_observed = findOmicsRuvUnobservedGroups(expr_matrix, design_matrix, sample_id_col, group_col))
}

#' The DA kernel's call with the assay's RUV factors in the model (KIT-129 B1): when the technical replicates'
#' correlation cannot be estimated, the refusal names k. The kernel's bound on k counts the replicate sets over the
#' injections, not the biological units, so a k within it can leave the units' correlation inestimable; the run's
#' handler (KIT-128) then records the assay as failed as for any other failure, with this message.
#' @param expr The DA kernel's call (evaluated here).
#' @param assay_factors The assay's entry of `runOmicsRuvFactors()$factors`, or NULL without RUV factors.
#' @return The value of `expr`. A failure of class `omics_da_block_correlation_failed` with k > 0 is raised again with
#'   k in its message and the class `omics_ruv_k_block_correlation_failed` added; every other condition passes unchanged.
#' @noRd
withOmicsRuvKOnBlockFailure <- function(expr, assay_factors) {
    k <- if (is.list(assay_factors) && !is.null(assay_factors$W)) ncol(assay_factors$W) else 0L
    if (k < 1L) return(expr)
    tryCatch(expr, omics_da_block_correlation_failed = function(e) {
        rlang::abort(sprintf("with k = %d RUV factors in the model, %s; use fewer factors", k, conditionMessage(e)),
            class = c("omics_ruv_k_block_correlation_failed", "omics_da_block_correlation_failed"), call = NULL)
    })
}

#' The DA kernel of an assay with RUV factors in the model (decision M29)
#'
#' The kernel's two halves as `runTestsContrastsLipidDA` and `runTestsContrastsMetabDA` call them, with each
#' comparison's standard error exact per feature (`ruvExactContrastsFit`, the exactness round's arm 25): with W in the
#' model, `limma::contrasts.fit` only approximates it for a feature with missing values. The technical-replicate block
#' (when the fit used one) is the fit half's, with its consensus correlation. An assay without factors keeps the
#' exported kernel.
#' @param fitFn,applyFn The omic's fit and contrast halves (`fitLipidDaAssayModel` and `applyLipidDaContrasts`, or the
#'   metabolomics pair).
#' @param data,contrast_strings,design_matrix,formula_string,sample_id_col,treat_lfc_cutoff,eBayes_trend,eBayes_robust,group_id,technical_replicate_id
#'   As `runTestsContrastsLipidDA`.
#' @return The value of `applyFn`.
#' @noRd
runOmicsRuvFactorsDaKernel <- function(fitFn, applyFn, data, contrast_strings, design_matrix, formula_string, sample_id_col,
                                       treat_lfc_cutoff, eBayes_trend, eBayes_robust, group_id, technical_replicate_id) {
    fit <- fitFn(data, design_matrix, formula_string, sample_id_col, group_id, technical_replicate_id)
    applyFn(fit, contrast_strings, treat_lfc_cutoff, eBayes_trend, eBayes_robust,
        contrastsFitFn = omicsRuvExactContrastsFitFn(fit, data, design_matrix, sample_id_col, group_id, technical_replicate_id))
}

#' The exact contrast step of an assay fit with RUV factors in the model (decision M29): `ruvExactContrastsFit` with the
#' refit of the fit's own model, on the fit's samples, with the fit half's technical-replicate block and consensus
#' correlation when it used one. Shared by the kernel (`runOmicsRuvFactorsDaKernel`) and the contrast declaration of
#' the run path (`omicsDaContrastOperation(binding, "declared")`, KIT-129 integration item 1).
#' @param fit The fit half's fit (`fitLipidDaAssayModel` or `fitMetabDaAssayModel`).
#' @param data,design_matrix,sample_id_col,group_id,technical_replicate_id What the fit half was given.
#' @return `function(fit, contrasts)` for the contrast half's `contrastsFitFn`.
#' @noRd
omicsRuvExactContrastsFitFn <- function(fit, data, design_matrix, sample_id_col, group_id, technical_replicate_id) {
    samples <- rownames(fit$design)
    data_subset <- data[, samples, drop = FALSE]
    block <- NULL
    if (isTRUE(fit$uses_duplicate_correlation)) {
        # the fit half's block: the same rule on the same design rows
        rows <- design_matrix
        if (sample_id_col %in% colnames(rows)) rownames(rows) <- rows[[sample_id_col]]
        block <- resolveProtTechnicalReplicateBlock(as.data.frame(rows)[samples, , drop = FALSE], group_id, technical_replicate_id)$block
    }
    correlation <- fit$consensus_correlation
    refit <- function(design) {
        if (is.null(block)) return(limma::lmFit(data_subset, design = design))
        limma::lmFit(data_subset, design = design, block = block, correlation = correlation)
    }
    \(fit, contrasts) ruvExactContrastsFit(fit, contrasts, fit$design, refit)
}
