# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ============================================================================
# func_metab_norm.R
# ============================================================================
# Purpose: Metabolomics normalization standalone functions (non-S4)
# 
# NOTE: S4 methods for metabolomics normalization (setMethod calls) remain in
# their source files (metabolite_normalization.R, metaboliteVsSamplesS4Objects.R)
# to maintain proper class-method coupling.
#
# This file is a placeholder for any future standalone (non-S4) metabolomics
# normalization helper functions that may be extracted.
#
# S4 methods for metabolomics normalization:
# - logTransformAssays() -> metabolite_normalization.R
# - normaliseUntransformedData() -> metabolite_normalization.R
# - normaliseBetweenSamples() -> metaboliteVsSamplesS4Objects.R
# - getNegCtrlMetabAnova() -> metaboliteVsSamplesS4Objects.R
#
# Dependencies:
# - limma
# - func_general_helpers.R (for utility functions)
# ============================================================================

# ============================================================================
# Helper Functions for Metabolomics Normalization Module
# ============================================================================



#' Generate QC Plots for Metabolomics Data and Save to Disk
#'
#' @title Generate QC plots for metabolomics data and save to disk
#' @description Generates PCA, RLE, Density, and Correlation plots for each assay
#'              in a MetaboliteAssayData object and saves them as PNG files.
#'
#' @param theObject MetaboliteAssayData S4 object
#' @param experiment_paths List of experiment directory paths (must contain metabolite_qc_dir)
#' @param stage One of "post_filter", "post_norm", "ruv_corrected"
#' @param grouping_variable Column name from design matrix for plot coloring
#' @param shape_variable Column name from design matrix for plot shapes (optional)
#' @param font_size Font size for plots
#'
#' @return Invisible list of file paths saved (named by assay and plot type)
#'
#' @importFrom ggplot2 ggsave
#' @importFrom purrr imap set_names
#' @importFrom logger log_info log_warn log_error
#' @noRd
generateMetabQcPlots <- function(
    theObject
    , experiment_paths
    , stage = c("post_filter", "post_norm", "ruv_corrected")
    , grouping_variable = NULL
    , shape_variable = NULL
    , font_size = 8
) {
    # --- Input validation ---
    stage <- match.arg(stage)

    stopifnot(
        inherits(theObject, "MetaboliteAssayData")
        , !is.null(experiment_paths$metabolite_qc_dir)
    )

    qc_dir <- experiment_paths$metabolite_qc_dir
    if (!dir.exists(qc_dir)) {
        dir.create(qc_dir, recursive = TRUE)
    }

    # --- Get assay names ---
    assay_list <- theObject@metabolite_data
    assay_names <- names(assay_list)

    if (length(assay_names) == 0) {
        logger::log_warn("No assays found in MetaboliteAssayData object")
        return(invisible(list()))
    }

    # --- Determine grouping variable ---
    if (is.null(grouping_variable)) {
        grouping_variable <- theObject@group_id
    }

    # --- Stage prefix mapping ---
    stage_prefix <- switch(
        stage
        , "post_filter" = "pre_norm"
        , "post_norm" = "post_norm"
        , "ruv_corrected" = "ruv_corrected"
    )

    # --- Generate plots for each assay ---
    saved_paths <- list()

    for (assay_name in assay_names) {
        safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
        # KIT-077: each file comes from the store; on a hit nothing below is plotted.
        write_figure <- function(path, type, width, height, render) {
            writeWorkflowOmicsFigure(path, render, "metabolomics",
                content = list(assay = assay_list[[assay_name]], design = theObject@design_matrix),
                figure = tools::file_path_sans_ext(basename(path)), assay = assay_name, type = type,
                parameters = list(stage = stage, grouping_variable = grouping_variable, shape_variable = shape_variable,
                    font_size = font_size), width = width, height = height, dpi = 150)
        }
        logger::log_info(paste("Generating QC plots for assay:", assay_name))

        # --- PCA Plot ---
        tryCatch({
            pca_path <- write_figure(file.path(qc_dir, sprintf("%s_%s_pca.png", safe_name, stage_prefix)), "pca", 8, 6, function(path) {
                pca_plots <- plotPca(
                    theObject
                    , grouping_variable = grouping_variable
                    , label_column = ""
                    , shape_variable = shape_variable
                    , title = ""
                    , font_size = font_size
                )

                # Select plot for this assay
                pca_plot <- if (is.list(pca_plots) && assay_name %in% names(pca_plots)) {
                    pca_plots[[assay_name]]
                } else if (inherits(pca_plots, "ggplot")) {
                    pca_plots
                } else if (is.list(pca_plots) && length(pca_plots) > 0) {
                    pca_plots[[1]]
                } else {
                    NULL
                }

                if (is.null(pca_plot)) return(NULL)
                ggplot2::ggsave(path, pca_plot, width = 8, height = 6, dpi = 150)
                omicsFigurePlotData(pca_plot)
            })
            if (!is.null(pca_path)) saved_paths[[paste0(assay_name, "_pca")]] <- pca_path
        }, error = function(e) {
            logger::log_warn(paste("Error generating PCA for", assay_name, ":", e$message))
            # KIT-131: a figure drawn but not written stops the step; any other failure costs this figure only.
            if (inherits(e, "workflow_figure_unwritten")) stop(e)
        })

        # --- RLE Plot ---
        tryCatch({
            rle_path <- write_figure(file.path(qc_dir, sprintf("%s_%s_rle.png", safe_name, stage_prefix)), "rle", 10, 6, function(path) {
                rle_plots <- plotRle(
                    theObject
                    , grouping_variable = grouping_variable
                    , yaxis_limit = c(-6, 6)
                )

                rle_plot <- if (is.list(rle_plots) && assay_name %in% names(rle_plots)) {
                    rle_plots[[assay_name]]
                } else if (inherits(rle_plots, "ggplot")) {
                    rle_plots
                } else if (is.list(rle_plots) && length(rle_plots) > 0) {
                    rle_plots[[1]]
                } else {
                    NULL
                }

                if (is.null(rle_plot)) return(NULL)
                ggplot2::ggsave(path, rle_plot, width = 10, height = 6, dpi = 150)
                omicsFigurePlotData(rle_plot)
            })
            if (!is.null(rle_path)) saved_paths[[paste0(assay_name, "_rle")]] <- rle_path
        }, error = function(e) {
            logger::log_warn(paste("Error generating RLE for", assay_name, ":", e$message))
            if (inherits(e, "workflow_figure_unwritten")) stop(e)
        })

        # --- Density Plot ---
        tryCatch({
            density_path <- write_figure(file.path(qc_dir, sprintf("%s_%s_density.png", safe_name, stage_prefix)), "density", 8, 6, function(path) {
                # Generate PCA first for density plot
                pca_for_density <- plotPca(
                    theObject
                    , grouping_variable = grouping_variable
                    , label_column = ""
                    , shape_variable = shape_variable
                    , title = ""
                    , font_size = font_size
                )

                density_plots <- plotDensity(
                    theObject
                    , grouping_variable = grouping_variable
                    , title = ""
                    , font_size = font_size
                )

                density_plot <- if (is.list(density_plots) && assay_name %in% names(density_plots)) {
                    density_plots[[assay_name]]
                } else if (inherits(density_plots, "ggplot")) {
                    density_plots
                } else if (is.list(density_plots) && length(density_plots) > 0) {
                    density_plots[[1]]
                } else {
                    NULL
                }

                if (is.null(density_plot)) return(NULL)
                ggplot2::ggsave(path, density_plot, width = 8, height = 6, dpi = 150)
                omicsFigurePlotData(density_plot)
            })
            if (!is.null(density_path)) saved_paths[[paste0(assay_name, "_density")]] <- density_path
        }, error = function(e) {
            logger::log_warn(paste("Error generating Density for", assay_name, ":", e$message))
            if (inherits(e, "workflow_figure_unwritten")) stop(e)
        })

        # --- Pearson Correlation Plot ---
        tryCatch({
            corr_path <- write_figure(file.path(qc_dir, sprintf("%s_%s_correlation.png", safe_name, stage_prefix)), "correlation", 10, 8, function(path) {
                corr_plots <- plotPearson(
                    theObject
                    , correlation_group = grouping_variable
                )

                corr_plot <- if (is.list(corr_plots) && assay_name %in% names(corr_plots)) {
                    corr_plots[[assay_name]]
                } else if (inherits(corr_plots, "ggplot")) {
                    corr_plots
                } else if (is.list(corr_plots) && length(corr_plots) > 0) {
                    corr_plots[[1]]
                } else {
                    NULL
                }

                if (is.null(corr_plot)) return(NULL)
                ggplot2::ggsave(path, corr_plot, width = 10, height = 8, dpi = 150)
                omicsFigurePlotData(corr_plot)
            })
            if (!is.null(corr_path)) saved_paths[[paste0(assay_name, "_correlation")]] <- corr_path
        }, error = function(e) {
            logger::log_warn(paste("Error generating Correlation for", assay_name, ":", e$message))
            if (inherits(e, "workflow_figure_unwritten")) stop(e)
        })
    }

    logger::log_info(paste("Generated", length(saved_paths), "QC plots for stage:", stage))
    return(invisible(saved_paths))
}


#' Which features of an omics assay are observed (finite) in every design sample: the controls RUV-III can use
#' for every sample, and so the upper bound on the factors it can estimate (`ruvCompleteControlRows`, KIT-005).
#' @return Logical vector, one entry per row of the assay, named by feature id.
#' @noRd
omicsRuvCompleteRows <- function(assay_data, id_column, design_samples) {
    sample_cols <- intersect(colnames(assay_data), design_samples)
    values <- suppressWarnings(vapply(assay_data[sample_cols], as.numeric, numeric(nrow(assay_data))))
    values <- matrix(values, nrow = nrow(assay_data))
    stats::setNames(rowSums(!is.finite(values)) == 0L, as.character(assay_data[[id_column]]))
}

#' How many of a control set are complete: by feature id when the set is named, else by position.
#' @return The count, or `NA_integer_` when the set cannot be aligned with the assay.
#' @noRd
omicsRuvCompleteControlCount <- function(complete_rows, ctrl_indices) {
    if (!is.null(names(ctrl_indices))) {
        return(as.integer(sum(complete_rows[names(ctrl_indices)[which(ctrl_indices)]], na.rm = TRUE)))
    }
    if (length(ctrl_indices) != length(complete_rows)) {
        return(NA_integer_)
    }
    as.integer(sum(complete_rows[which(ctrl_indices)], na.rm = TRUE))
}

#' The object a lipid or metabolite RUV sweep hands to `ruvCancor`: assays with missing values carry their NIPALS
#' imputation of the matrix `ruvCancor` would build, so the method finds nothing left to impute. The imputation does
#' not depend on the percentage, so a sweep prepares it once per assay (KIT-005's rule for proteomics). An assay this
#' cannot prepare stays as it is, and `ruvCancor` handles and reports it as before. Only the canonical-correlation
#' plots read this object: the control selection and the correction keep the data with its missing values.
#' @param theObject A `MetaboliteAssayData` or `LipidomicsAssayData`.
#' @param data_slot,id_column_slot The names of the class's assay list and feature id slots.
#' @return `theObject`, with imputed values in the assays that had missing values.
#' @noRd
prepareOmicsRuvCancorObject <- function(theObject, data_slot, id_column_slot) {
    assays <- methods::slot(theObject, data_slot)
    id_column <- methods::slot(theObject, id_column_slot)
    design_samples <- as.character(theObject@design_matrix[[theObject@sample_id]])
    # The value `ruvCancor` resolves when a sweep calls it: the object's stored argument, else its default.
    num_components <- theObject@args[["ruvCancor"]][["num_components_to_impute"]] %||% 2
    for (assay_name in names(assays)) {
        assay <- assays[[assay_name]]
        sample_cols <- intersect(colnames(assay), design_samples)
        if (!id_column %in% colnames(assay) || length(sample_cols) < 2L) next
        values <- as.matrix(as.data.frame(lapply(assay[sample_cols], as.numeric), check.names = FALSE))
        values[!is.finite(values)] <- NA_real_
        if (!anyNA(values)) next
        ordered_samples <- design_samples[design_samples %in% sample_cols]
        imputed <- tryCatch(
            imputeNipalsMatrix(t(values[, ordered_samples, drop = FALSE]), ncomp = num_components),
            error = function(e) NULL
        )
        if (is.null(imputed)) next
        assay[ordered_samples] <- as.data.frame(t(imputed), check.names = FALSE)
        assays[[assay_name]] <- assay
    }
    methods::slot(theObject, data_slot) <- assays
    theObject
}

#' Whether lipid and metabolite RUV handle a sample, a feature or an assay that cannot be corrected by its cause
#' (KIT-242 Part B; VAL-242 OD-1 (c), Q6.1, Q6.2 (c); VAL-290): TRUE under the current rule, FALSE only while a record
#' made before it replays (`workflowReplayRuleMarkers$ruv_fallback`), when an assay whose correction failed was dropped.
#' @return One logical.
#' @noRd
omicsRuvFallbackCurrent <- function() {
    marker <- workflowReplayRuleMarkers$ruv_fallback
    !identical(getOption(marker$option), marker$earlier)
}

#' Stop a lipid or metabolite RUV step (KIT-242, Q6.2 (c)): a failed computation or an invalid argument, with the assays
#' and the step that failed (`reason`: `invalid_k`, `invalid_controls`, `malformed_assay`, `too_few_samples`,
#' `correction_error`, `search_error`, or the search's `error_*` status).
#' @noRd
abortOmicsRuvFailed <- function(assay, reason, text, parent = NULL) {
    rlang::abort(text, class = "omics_ruv_failed", assay = assay, reason = reason, parent = parent, call = NULL)
}

#' The samples and features with no finite value in each assay of a lipid or metabolite object, over the design's
#' samples (KIT-242, VAL-242 Q6.1 and OD-1 (c)), and the object without them: RUV computes each assay without them and
#' puts them back as missing. An assay with neither is unchanged.
#' @param theObject A `LipidomicsAssayData` or `MetaboliteAssayData` object.
#' @return `list(object, samples, features)`; `samples` and `features` per assay (character, empty where none).
#' @noRd
omicsRuvWithoutEmpty <- function(theObject) {
    slots <- omicsRuvAssaySlots(theObject)
    assays <- methods::slot(theObject, slots$data)
    id_column <- methods::slot(theObject, slots$id)
    design_samples <- as.character(theObject@design_matrix[[theObject@sample_id]])
    samples <- stats::setNames(rep(list(character()), length(assays)), names(assays))
    features <- samples
    for (assay_name in names(assays)) {
        assay <- assays[[assay_name]]
        cols <- intersect(colnames(assay), design_samples)
        if (!length(cols) || !nrow(assay) || !id_column %in% colnames(assay)) next
        seen <- matrix(vapply(assay[cols], function(x) is.finite(suppressWarnings(as.numeric(x))), logical(nrow(assay))),
            nrow = nrow(assay))
        empty_cols <- cols[colSums(seen) == 0L]
        empty_rows <- rowSums(seen) == 0L
        if (!length(empty_cols) && !any(empty_rows)) next
        samples[[assay_name]] <- empty_cols
        features[[assay_name]] <- as.character(assay[[id_column]][empty_rows])
        assays[[assay_name]] <- assay[!empty_rows, setdiff(colnames(assay), empty_cols), drop = FALSE]
    }
    methods::slot(theObject, slots$data) <- assays
    list(object = theObject, samples = samples, features = features)
}

#' One assay after RUV with its left-out samples and features put back as missing (KIT-242): the corrected table, plus
#' an all-missing column for each left-out sample and an all-missing row (the input's annotation) for each left-out
#' feature, in the input's row and sample-column order.
#' @noRd
omicsRuvRestoreEmpty <- function(corrected, input, id_column, sample_cols, samples, features) {
    if (!length(samples) && !length(features)) return(corrected)
    for (sample in samples) corrected[[sample]] <- NA_real_
    if (length(features)) {
        extra <- input[as.character(input[[id_column]]) %in% features, intersect(colnames(corrected), colnames(input)), drop = FALSE]
        for (column in intersect(colnames(extra), sample_cols)) extra[[column]] <- NA_real_
        corrected <- dplyr::bind_rows(corrected, extra)
        corrected <- corrected[order(match(as.character(corrected[[id_column]]), as.character(input[[id_column]]))), , drop = FALSE]
    }
    present <- sample_cols[sample_cols %in% colnames(corrected)]
    corrected[, c(setdiff(colnames(corrected), sample_cols), present), drop = FALSE]
}

#' Legacy RUV-III-C on a lipid or metabolite object by cause (KIT-242 Part B; VAL-242 OD-1 (c), OD-3 (d), Q6.1,
#' Q6.2 (c)): the `ruvIII_C_Varying` methods' body under the current rule. Per assay: an invalid k or control argument,
#' or a malformed assay, stops (`omics_ruv_failed`); the samples and features with no value are left out and put back
#' as missing after the correction (`omics_ruv_sample_excluded`, `omics_ruv_feature_excluded`); a feature the control
#' vector does not name is not a control; an assay with no control feature keeps its values without RUV
#' (`omics_ruv_not_applied`, reason `no_usable_controls`); the correction runs under the earlier rule on what remains,
#' and an assay it cannot correct (an error, no matrix, no finite value) stops (`correction_error`). No assay is
#' dropped. With nothing left out, every assay corrected and the controls aligned, the result is the earlier rule's.
#' A NULL k or ctrl (read from the object's arguments) is handed to the earlier rule unchanged.
#' @return The corrected object.
#' @noRd
omicsRuvCorrectByCause <- function(theObject, ruv_grouping_variable, ruv_number_k, ctrl) {
    marker <- workflowReplayRuleMarkers$ruv_fallback
    earlier <- stats::setNames(list(marker$earlier), marker$option)
    inner <- function(object, k, controls) withr::with_options(earlier,
        ruvIII_C_Varying(object, ruv_grouping_variable = ruv_grouping_variable, ruv_number_k = k, ctrl = controls))
    slots <- omicsRuvAssaySlots(theObject)
    assays <- methods::slot(theObject, slots$data)
    id_column <- methods::slot(theObject, slots$id)
    is_k_list <- is.list(ruv_number_k)
    is_ctrl_list <- is.list(ctrl) && !is.data.frame(ctrl)
    assay_names <- names(assays)
    if (is.null(ruv_number_k) || is.null(ctrl) || !length(assays) || is.null(assay_names) || any(!nzchar(assay_names)) ||
        (is_k_list && !all(assay_names %in% names(ruv_number_k))) || (is_ctrl_list && !all(assay_names %in% names(ctrl)))) {
        return(inner(theObject, ruv_number_k, ctrl))
    }
    design_samples <- as.character(theObject@design_matrix[[theObject@sample_id]])
    plan <- list()
    for (assay_name in assay_names) {
        k <- if (is_k_list) ruv_number_k[[assay_name]] else ruv_number_k
        if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 0) {
            abortOmicsRuvFailed(assay_name, "invalid_k", sprintf(paste("RUV-III was not run: the number of unwanted factors for",
                "assay '%s' must be one whole number of 0 or more (given: %s)."), assay_name, paste(format(k), collapse = ", ")))
        }
        assay <- assays[[assay_name]]
        sample_cols <- if (is.data.frame(assay)) intersect(colnames(assay), design_samples) else character()
        if (!is.data.frame(assay) || !id_column %in% colnames(assay) || !nrow(assay) || length(sample_cols) < 2L) {
            abortOmicsRuvFailed(assay_name, "malformed_assay", sprintf(paste("RUV-III was not run: assay '%s' needs the feature",
                "id column '%s', at least one feature and at least two samples of the design."), assay_name, id_column))
        }
        if (!tibble::is_tibble(assay)) assay <- tibble::as_tibble(assay)
        ids <- as.character(assay[[id_column]])
        seen <- matrix(vapply(assay[sample_cols], function(x) is.finite(suppressWarnings(as.numeric(x))), logical(nrow(assay))),
            nrow = nrow(assay))
        empty_samples <- sample_cols[colSums(seen) == 0L]
        empty_features <- unique(ids[rowSums(seen) == 0L])
        if (length(sample_cols) - length(empty_samples) < 2L || length(empty_features) >= length(unique(ids))) {
            abortOmicsRuvFailed(assay_name, "too_few_samples", sprintf(paste("RUV-III was not run: assay '%s' has values in fewer",
                "than two samples or in no feature."), assay_name))
        }
        if (length(empty_samples)) {
            rlang::warn(sprintf(paste("RUV-III for assay '%s' was computed without samples %s, which have no value there; they",
                "stay empty in that assay."), assay_name, paste(empty_samples, collapse = ", ")),
                class = "omics_ruv_sample_excluded", assay = assay_name, samples = empty_samples)
        }
        if (length(empty_features)) {
            rlang::warn(sprintf(paste("RUV-III for assay '%s' was computed without %d features with no value (%s); they stay",
                "empty."), assay_name, length(empty_features), paste(utils::head(empty_features, 3L), collapse = ", ")),
                class = "omics_ruv_feature_excluded", assay = assay_name, features = empty_features)
        }
        kept_ids <- setdiff(ids, empty_features)
        # VAL-242 Q6.1: the controls by feature id; a feature the control vector does not name is not a control.
        given <- if (is_ctrl_list) ctrl[[assay_name]] else ctrl
        invalid_controls <- function(what) abortOmicsRuvFailed(assay_name, "invalid_controls", sprintf(paste("RUV-III was not run:",
            "the control features for assay '%s' %s."), assay_name, what))
        aligned <- TRUE
        control_ids <- if (is.null(given)) character() else if (is.numeric(given)) {
            if (any(is.na(given)) || any(given < 1) || any(given > length(ids))) invalid_controls("are positions outside the assay")
            ids[given]
        } else if (is.logical(given)) {
            if (length(given) == length(ids)) {
                ids[given %in% TRUE]
            } else if (!is.null(names(given))) {
                aligned <- all(ids %in% names(given))
                names(given)[given %in% TRUE]
            } else {
                invalid_controls("are an unnamed logical vector of another length than the assay")
            }
        } else if (is.character(given)) {
            given
        } else invalid_controls("must be feature ids, positions or a logical vector")
        control_ids <- intersect(kept_ids, control_ids)
        changed <- length(empty_samples) || length(empty_features) || !aligned || !tibble::is_tibble(assays[[assay_name]])
        plan[[assay_name]] <- list(k = k, assay = assay, sample_cols = sample_cols, samples = empty_samples,
            features = empty_features, controls = control_ids, changed = changed)
        if (!length(control_ids)) {
            rlang::warn(sprintf(paste("RUV-III was not applied to assay %s: none of its features with values is a control",
                "feature. The assay keeps its values without RUV."), assay_name),
                class = "omics_ruv_not_applied", assay = assay_name, reason = "no_usable_controls")
        }
    }
    to_correct <- names(Filter(function(p) length(p$controls) > 0L, plan))
    if (!length(to_correct)) {
        theObject@args$ruv_grouping_variable <- ruv_grouping_variable
        theObject@args$ruv_number_k <- ruv_number_k
        theObject@args$ctrl <- ctrl
        return(theObject)
    }
    direct <- identical(to_correct, assay_names) && !any(vapply(plan, `[[`, logical(1L), "changed"))
    if (direct) {
        result <- inner(theObject, ruv_number_k, ctrl)
        inputs <- assays
    } else {
        inputs <- lapply(plan[to_correct], function(p) p$assay[!as.character(p$assay[[id_column]]) %in% p$features,
            setdiff(colnames(p$assay), p$samples), drop = FALSE])
        sub <- theObject
        methods::slot(sub, slots$data) <- inputs
        result <- inner(sub, lapply(plan[to_correct], `[[`, "k"), lapply(plan[to_correct], `[[`, "controls"))
    }
    result_assays <- methods::slot(result, slots$data)
    failed <- Filter(function(assay_name) {
        out <- result_assays[[assay_name]]
        is.null(out) || (plan[[assay_name]]$k > 0 && identical(out, inputs[[assay_name]]))
    }, to_correct)
    if (length(failed)) {
        abortOmicsRuvFailed(failed, "correction_error", sprintf(paste("RUV-III-C could not correct assay %s (an error, no",
            "corrected matrix, or no finite value; see the log). Choose another k or control set, or run without RUV."),
            paste(failed, collapse = ", ")))
    }
    if (direct) return(result)
    final <- result
    methods::slot(final, slots$data) <- stats::setNames(lapply(assay_names, function(assay_name) {
        p <- plan[[assay_name]]
        if (!assay_name %in% to_correct) return(assays[[assay_name]])
        omicsRuvRestoreEmpty(result_assays[[assay_name]], p$assay, id_column, p$sample_cols, p$samples, p$features)
    }), assay_names)
    final@design_matrix <- theObject@design_matrix
    final@args$ruv_grouping_variable <- ruv_grouping_variable
    final@args$ruv_number_k <- ruv_number_k
    final@args$ctrl <- ctrl
    tryCatch(cleanDesignMatrix(final), error = function(e) final)
}

#' Run Per-Assay RUV Optimization
#'
#' @title Run RUV optimization per assay
#' @description Automatic mode uses a percentage-first whole-object search loop:
#'   \code{getNegCtrlMetabAnova()} and \code{ruvCancor()} are called once per
#'   tested percentage (not once per assay), and per-assay result rows are then
#'   derived from the shared whole-object results.  Manual mode calls each helper
#'   once and uses the caller-supplied k and percentage directly.
#'
#' @param theObject MetaboliteAssayData S4 object
#' @param ruv_mode One of "automatic", "manual"
#' @param params List of RUV parameters:
#'   - percentage_min: Minimum percentage for auto mode (default 1)
#'   - percentage_max: Maximum percentage for auto mode (default 20)
#'   - max_acceptable_k: Hard cap on acceptable k (default 3, overridden by
#'     adaptive_k_penalty)
#'   - ruv_grouping_variable: Column name for grouping
#'   - separation_metric: One of "max_difference", "mean_difference", "auc",
#'     "weighted_difference" (deprecated)
#'   - k_penalty_weight: Penalty weight from 0 to 1 (default 0.5)
#'   - adaptive_k_penalty: Use sample-size-adaptive max_acceptable_k (default TRUE)
#'   - manual_k: k value for manual mode (default 2)
#'   - manual_percentage: Percentage for manual mode (default 10)
#' @param experiment_paths List of experiment directory paths (currently unused,
#'   kept for interface stability)
#'
#' @return Named list (one entry per assay).  Automatic-mode entries include:
#'   \itemize{
#'     \item success, best_k, best_percentage, best_realized_num_controls,
#'       best_realized_percentage, control_genes_index, cancor_plot,
#'       separation_score, composite_score, optimization_results, error
#'   }
#'   Manual-mode entries include: success, best_k, best_percentage,
#'   control_genes_index, cancor_plot, separation_score, optimization_results, error.
#'
#' @importFrom purrr map imap map_dfr
#' @importFrom dplyr filter arrange desc
#' @importFrom logger log_info log_warn
#' @noRd
runPerAssayRuvOptimization <- function(
    theObject
    , ruv_mode = c("automatic", "manual")
    , params = list()
    , experiment_paths = NULL
) {

    ruv_mode <- match.arg(ruv_mode)

    stopifnot(inherits(theObject, "MetaboliteAssayData"))

    # --- Extract parameters ---
    percentage_min       <- params$percentage_min       %||% 1
    percentage_max       <- params$percentage_max       %||% 20
    max_acceptable_k     <- params$max_acceptable_k     %||% 3L
    ruv_grouping_variable <- params$ruv_grouping_variable
    separation_metric    <- params$separation_metric    %||% "max_difference"
    k_penalty_weight     <- params$k_penalty_weight     %||% 0.5
    adaptive_k_penalty   <- params$adaptive_k_penalty   %||% TRUE
    manual_k             <- params$manual_k             %||% 2
    manual_percentage    <- params$manual_percentage    %||% 10

    # Emit deprecation warning once at public entry boundary
    if (ruv_mode == "automatic" && identical(separation_metric, "weighted_difference")) {
        warning(
            paste0(
                "'weighted_difference' is deprecated because it up-weights high-K values "
                , "while the composite score penalizes high K. "
                , "Use 'max_difference' or 'mean_difference' instead."
            )
            , call. = FALSE
        )
    }


    # --- Validate grouping variable ---
    if (is.null(ruv_grouping_variable) || !ruv_grouping_variable %in% colnames(theObject@design_matrix)) {
        ruv_grouping_variable <- theObject@group_id
    }

    assay_list  <- theObject@metabolite_data
    assay_names <- names(assay_list)
    logger::log_info(paste("Running RUV optimization for", length(assay_names), "assays in", ruv_mode, "mode"))

    # =========================================================================
    # AUTOMATIC MODE
    # =========================================================================
    if (ruv_mode == "automatic") {

        percentage_range <- seq(percentage_min, percentage_max, by = 1)

        # --- Percentage-first whole-object loop ---
        # getNegCtrlMetabAnova() and ruvCancor() are called ONCE PER PERCENTAGE,
        # not once per assay.  Both return whole-object named lists.
        cancor_object <- tryCatch(
            prepareOmicsRuvCancorObject(theObject, "metabolite_data", "metabolite_id_column")
            , error = function(e) theObject
        )
        percentage_results <- purrr::map(percentage_range, function(pct) {

            ctrl_list <- tryCatch(
                getNegCtrlMetabAnova(
                    theObject
                    , ruv_grouping_variable = ruv_grouping_variable
                    , percentage_as_neg_ctrl = pct
                )
                , error = function(e) e
            )

            if (inherits(ctrl_list, "error")) {
                return(list(
                    percentage_requested = pct
                    , object_error = paste0("error_neg_ctrl_selection: ", ctrl_list$message)
                    , ctrl_list   = NULL
                    , cancor_list = NULL
                    , refusal     = if (inherits(ctrl_list, "multischolar_ruv_grouping_unusable")) ctrl_list
                ))
            }

            cancor_list <- tryCatch(
                ruvCancor(
                    cancor_object
                    , ctrl = ctrl_list
                    , ruv_grouping_variable = ruv_grouping_variable
                )
                , error = function(e) e
            )

            if (inherits(cancor_list, "error")) {
                return(list(
                    percentage_requested = pct
                    , object_error = paste0("error_cancor: ", cancor_list$message)
                    , ctrl_list   = ctrl_list
                    , cancor_list = NULL
                ))
            }

            list(
                percentage_requested = pct
                , object_error = NULL
                , ctrl_list   = ctrl_list
                , cancor_list = cancor_list
            )
        })

        # A grouping that cannot be used is refused to the caller, not recorded as failed rows (KIT-190). Signalled
        # here, outside purrr::map, so the refusal reaches the caller with its own class and message.
        for (pct_res in percentage_results) {
            if (!is.null(pct_res$refusal)) {
                stop(pct_res$refusal)
            }
        }

        # --- Per-assay derivation from shared percentage results ---
        design_samples <- as.character(theObject@design_matrix[[theObject@sample_id]])

        results <- purrr::imap(assay_list, function(assay_data, assay_name) {

            # Sample size: assay columns that appear in the design matrix
            sample_size <- length(intersect(colnames(assay_data), design_samples))
            candidate_feature_count <- nrow(assay_data)

            effective_max_k <- if (adaptive_k_penalty) {
                calculateAdaptiveMaxK(sample_size)
            } else {
                as.integer(max_acceptable_k)
            }

            # Build one result row per tested percentage for this assay
            complete_rows <- omicsRuvCompleteRows(assay_data, theObject@metabolite_id_column, design_samples)

            assay_pct_rows <- purrr::map(percentage_results, function(pct_res) {
                pct <- pct_res$percentage_requested

                # Whole-object failure → failed row for this assay, run continues
                if (!is.null(pct_res$object_error)) {
                    err_status <- if (startsWith(pct_res$object_error, "error_neg_ctrl_selection")) {
                        "error_neg_ctrl_selection"
                    } else {
                        "error_cancor"
                    }
                    return(list(
                        percentage_requested    = pct
                        , candidate_feature_count = candidate_feature_count
                        , realized_num_controls   = NA_integer_
                        , realized_percentage     = NA_real_
                        , sample_size             = sample_size
                        , best_k                  = NA_integer_
                        , separation_score        = NA_real_
                        , composite_score         = NA_real_
                        , status                  = err_status
                        , error_reason            = pct_res$object_error
                        , ctrl_indices            = NULL
                        , cancor_plot             = NULL
                    ))
                }

                # KIT-242 (VAL-373 Q-A): under the current rule an assay the control selection skipped has no control set
                # of its own (the earlier rule took the whole list, on which sum() below fails); the reason travels on.
                if (omicsRuvFallbackCurrent() && is.list(pct_res$ctrl_list) && !is.null(names(pct_res$ctrl_list)) &&
                    !assay_name %in% names(pct_res$ctrl_list)) {
                    skip <- attr(pct_res$ctrl_list, "skipped")[[assay_name]] %||% list(reason = "not_selected")
                    return(list(
                        percentage_requested    = pct
                        , candidate_feature_count = candidate_feature_count
                        , realized_num_controls   = 0L
                        , realized_percentage     = 0
                        , sample_size             = sample_size
                        , best_k                  = NA_integer_
                        , separation_score        = NA_real_
                        , composite_score         = NA_real_
                        , status                  = "skipped_no_controls_selected"
                        , error_reason            = skip$reason
                        , ctrl_indices            = NULL
                        , cancor_plot             = NULL
                        , skip_reason             = skip
                    ))
                }

                # Extract per-assay slices from whole-object results
                ctrl_indices <- if (is.list(pct_res$ctrl_list) && assay_name %in% names(pct_res$ctrl_list)) {
                    pct_res$ctrl_list[[assay_name]]
                } else {
                    pct_res$ctrl_list
                }

                cancor_plot <- if (is.list(pct_res$cancor_list) && assay_name %in% names(pct_res$cancor_list)) {
                    pct_res$cancor_list[[assay_name]]
                } else {
                    pct_res$cancor_list
                }

                realized_num_controls <- sum(ctrl_indices, na.rm = TRUE)
                realized_percentage <- if (candidate_feature_count > 0L) {
                    100 * realized_num_controls / candidate_feature_count
                } else {
                    NA_real_
                }

                if (realized_num_controls < 5L) {
                    return(list(
                        percentage_requested    = pct
                        , candidate_feature_count = candidate_feature_count
                        , realized_num_controls   = realized_num_controls
                        , realized_percentage     = realized_percentage
                        , sample_size             = sample_size
                        , best_k                  = NA_integer_
                        , separation_score        = NA_real_
                        , composite_score         = NA_real_
                        , status                  = "skipped_insufficient_controls"
                        , error_reason            = ""
                        , ctrl_indices            = ctrl_indices
                        , cancor_plot             = NULL
                    ))
                }

                # KIT-242 (VAL-242 Q6.2 (c)): under the current rule an assay the canonical-correlation step returned no
                # plot for (it could not be computed) is a failed computation, not a curve without an elbow.
                if (omicsRuvFallbackCurrent() && is.list(pct_res$cancor_list) && !assay_name %in% names(pct_res$cancor_list)) {
                    return(list(
                        percentage_requested    = pct
                        , candidate_feature_count = candidate_feature_count
                        , realized_num_controls   = realized_num_controls
                        , realized_percentage     = realized_percentage
                        , sample_size             = sample_size
                        , best_k                  = NA_integer_
                        , separation_score        = NA_real_
                        , composite_score         = NA_real_
                        , status                  = "error_cancor"
                        , error_reason            = "error_cancor: the canonical correlations could not be computed for this assay"
                        , ctrl_indices            = ctrl_indices
                        , cancor_plot             = NULL
                    ))
                }

                separation_score <- tryCatch(
                    calculateSeparationScore(cancor_plot, metric = separation_metric)
                    , error = function(e) NA_real_
                )

                best_k <- tryCatch(
                    findBestKElbow(cancor_plot)
                    , error = function(e) NA_integer_
                )

                if (is.na(best_k)) {
                    return(list(
                        percentage_requested    = pct
                        , candidate_feature_count = candidate_feature_count
                        , realized_num_controls   = realized_num_controls
                        , realized_percentage     = realized_percentage
                        , sample_size             = sample_size
                        , best_k                  = NA_integer_
                        , separation_score        = separation_score
                        , composite_score         = NA_real_
                        , status                  = "invalid_cancor_plot"
                        , error_reason            = "findBestKElbow returned NA"
                        , ctrl_indices            = ctrl_indices
                        , cancor_plot             = cancor_plot
                    ))
                }

                best_k <- as.integer(best_k)
                complete_controls <- omicsRuvCompleteControlCount(complete_rows, ctrl_indices)
                # RUV-III cannot estimate more factors than it has controls observed in every sample (KIT-005).
                if (isTRUE(best_k > complete_controls)) {
                    return(list(
                        percentage_requested    = pct
                        , candidate_feature_count = candidate_feature_count
                        , realized_num_controls   = realized_num_controls
                        , realized_percentage     = realized_percentage
                        , sample_size             = sample_size
                        , best_k                  = best_k
                        , separation_score        = separation_score
                        , composite_score         = NA_real_
                        , status                  = "skipped_k_exceeds_complete_controls"
                        , error_reason            = sprintf("k = %d exceeds %d controls complete in every sample", best_k, complete_controls)
                        , ctrl_indices            = ctrl_indices
                        , cancor_plot             = cancor_plot
                    ))
                }
                composite_score <- tryCatch(
                    calculateCompositeScore(separation_score, best_k, k_penalty_weight, effective_max_k)
                    , error = function(e) NA_real_
                )

                if (!is.finite(composite_score)) {
                    return(list(
                        percentage_requested    = pct
                        , candidate_feature_count = candidate_feature_count
                        , realized_num_controls   = realized_num_controls
                        , realized_percentage     = realized_percentage
                        , sample_size             = sample_size
                        , best_k                  = best_k
                        , separation_score        = separation_score
                        , composite_score         = NA_real_
                        , status                  = "error_scoring"
                        , error_reason            = "composite score is non-finite"
                        , ctrl_indices            = ctrl_indices
                        , cancor_plot             = cancor_plot
                    ))
                }

                list(
                    percentage_requested    = pct
                    , candidate_feature_count = candidate_feature_count
                    , realized_num_controls   = realized_num_controls
                    , realized_percentage     = realized_percentage
                    , sample_size             = sample_size
                    , best_k                  = best_k
                    , separation_score        = separation_score
                    , composite_score         = composite_score
                    , status                  = "ok"
                    , error_reason            = ""
                    , ctrl_indices            = ctrl_indices
                    , cancor_plot             = cancor_plot
                )
            })

            # Build optimization_results trace (one row per tested percentage)
            optimization_results <- purrr::map_dfr(assay_pct_rows, function(row) {
                data.frame(
                    percentage_requested    = row$percentage_requested
                    , candidate_feature_count = row$candidate_feature_count
                    , realized_num_controls   = row$realized_num_controls %||% NA_integer_
                    , realized_percentage     = row$realized_percentage   %||% NA_real_
                    , sample_size             = row$sample_size
                    , best_k                  = row$best_k                %||% NA_integer_
                    , separation_score        = row$separation_score      %||% NA_real_
                    , composite_score         = row$composite_score       %||% NA_real_
                    , status                  = row$status
                    , error_reason            = row$error_reason
                    , stringsAsFactors = FALSE
                )
            })

            # Deterministic winner selection: highest composite → lowest k →
            # highest separation → lowest percentage_requested
            valid_df <- optimization_results |>
                dplyr::filter(is.finite(composite_score)) |>
                dplyr::arrange(
                    dplyr::desc(composite_score)
                    , best_k
                    , dplyr::desc(separation_score)
                    , percentage_requested
                )

            if (nrow(valid_df) == 0L) {
                logger::log_warn(paste("RUV: no valid percentage for assay:", assay_name))
                refusal <- ruvNoValidPercentageMessage(optimization_results$percentage_requested
                    , optimization_results$status
                    , optimization_results$best_k
                    , vapply(assay_pct_rows, function(row) {
                        if (identical(row$status, "skipped_k_exceeds_complete_controls")) omicsRuvCompleteControlCount(complete_rows, row$ctrl_indices) else NA_integer_
                    }, integer(1L))
                    , "metabolites")
                failed <- list(
                    success                  = FALSE
                    , best_k                 = NA_integer_
                    , best_percentage        = NA_real_
                    , best_realized_num_controls = NA_integer_
                    , best_realized_percentage   = NA_real_
                    , control_genes_index    = NULL
                    , cancor_plot            = NULL
                    , separation_score       = NA_real_
                    , composite_score        = NA_real_
                    , optimization_results   = optimization_results
                    , error                  = refusal %||% "No valid percentage found"
                )
                # KIT-242 (VAL-373 Q-A): why the control selection skipped the assay, when it did (else no field).
                skip <- Find(Negate(is.null), lapply(assay_pct_rows, `[[`, "skip_reason"))
                if (!is.null(skip)) failed$skip_reason <- skip
                return(failed)
            }

            best_pct <- valid_df$percentage_requested[[1L]]
            winner_idx <- which(
                vapply(assay_pct_rows, function(r) r$percentage_requested == best_pct, logical(1L))
            )[[1L]]
            winner <- assay_pct_rows[[winner_idx]]


            list(
                success                      = TRUE
                , best_k                     = as.integer(winner$best_k)
                , best_percentage            = winner$percentage_requested
                , best_realized_num_controls = as.integer(winner$realized_num_controls)
                , best_realized_percentage   = winner$realized_percentage
                , control_genes_index        = winner$ctrl_indices
                , cancor_plot                = winner$cancor_plot
                , separation_score           = winner$separation_score
                , composite_score            = winner$composite_score
                , optimization_results       = optimization_results
                , error                      = NULL
            )
        })

        names(results) <- assay_names
        return(results)
    }

    # =========================================================================
    # MANUAL MODE (unchanged semantics)
    # =========================================================================
    results <- purrr::imap(assay_list, function(assay_data, assay_name) {
        logger::log_info(paste("Processing assay:", assay_name))

        tryCatch({
            ctrl_result <- getNegCtrlMetabAnova(
                theObject
                , ruv_grouping_variable = ruv_grouping_variable
                , percentage_as_neg_ctrl = manual_percentage
            )

            # KIT-242 (VAL-373 Q-A): under the current rule an assay the control selection skipped fails the Manual
            # search with its reason (the earlier rule took the whole list as that assay's controls).
            if (omicsRuvFallbackCurrent() && is.list(ctrl_result) && !is.null(names(ctrl_result)) &&
                !assay_name %in% names(ctrl_result)) {
                skip <- attr(ctrl_result, "skipped")[[assay_name]] %||% list(reason = "not_selected")
                return(list(
                    success               = FALSE
                    , best_k              = NA_integer_
                    , best_percentage     = NA_real_
                    , control_genes_index = NULL
                    , cancor_plot         = NULL
                    , separation_score    = NA_real_
                    , optimization_results = NULL
                    , error               = paste("the control selection skipped the assay:", skip$reason)
                    , skip_reason         = skip
                ))
            }

            ctrl_indices <- if (is.list(ctrl_result) && assay_name %in% names(ctrl_result)) {
                ctrl_result[[assay_name]]
            } else {
                ctrl_result
            }

            cancor_result <- ruvCancor(
                theObject
                , ctrl = ctrl_indices
                , ruv_grouping_variable = ruv_grouping_variable
            )

            cancor_plot <- if (is.list(cancor_result) && assay_name %in% names(cancor_result)) {
                cancor_result[[assay_name]]
            } else {
                cancor_result
            }


            list(
                success               = TRUE
                , best_k              = manual_k
                , best_percentage     = manual_percentage
                , control_genes_index = ctrl_indices
                , cancor_plot         = cancor_plot
                , separation_score    = NA_real_
                , optimization_results = NULL
                , error               = NULL
            )

        }, error = function(e) {
            if (inherits(e, "multischolar_ruv_grouping_unusable")) {
                return(list(refusal = e))
            }
            logger::log_warn(paste("RUV manual mode failed for", assay_name, ":", e$message))
            list(
                success               = FALSE
                , best_k              = NA_integer_
                , best_percentage     = NA_real_
                , control_genes_index = NULL
                , cancor_plot         = NULL
                , separation_score    = NA_real_
                , optimization_results = NULL
                , error               = e$message
            )
        })
    })

    # A grouping that cannot be used is refused to the caller (KIT-190); signalled outside purrr::imap so the
    # refusal keeps its own class and message.
    for (res in results) {
        if (!is.null(res$refusal)) {
            stop(res$refusal)
        }
    }

    names(results) <- assay_names
    return(results)
}






