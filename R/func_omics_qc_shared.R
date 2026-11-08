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
# func_omics_qc_shared.R
# ============================================================================
# Purpose: Shared Quality Control functions for multiple omics types
#
# This file contains QC functions that are common across different omics
# workflows (e.g., lipidomics, metabolomics) or have been unified to handle
# multiple object types to avoid code duplication and namespace conflicts.
#
# Dependencies:
# - dplyr
# - purrr
# - methods
# ============================================================================

# ----------------------------------------------------------------------------
# findDuplicateFeatureIDs (Unified)
# ----------------------------------------------------------------------------
#' @title Find Duplicate Feature IDs in Omics Assay Data
#'
#' @description
#' Identifies duplicate feature IDs within each assay of an omics data object.
#' Supports both LipidomicsAssayData and MetaboliteAssayData objects.
#'
#' @param theObject A LipidomicsAssayData or MetaboliteAssayData object.
#'
#' @return A named list where each element corresponds to an assay and contains
#'         a tibble of duplicate feature IDs (if any), or NULL.
#'
#' @importFrom dplyr count filter pull %>%
#' @importFrom purrr map set_names
#' @importFrom methods slot
#' @importFrom rlang sym !!
#'
#' @examples
#' \dontrun{
#' # For Lipidomics
#' duplicates_lipid <- findDuplicateFeatureIDs(lipid_assay_obj)
#'
#' # For Metabolomics
#' duplicates_metab <- findDuplicateFeatureIDs(metab_assay_obj)
#' }
#' @export
findDuplicateFeatureIDs <- function(theObject) {
    if (inherits(theObject, "LipidomicsAssayData")) {
        return(findLipidDuplicateFeatureIDs(theObject))
    } else if (inherits(theObject, "MetaboliteAssayData")) {
        return(findMetabDuplicateFeatureIDs(theObject))
    } else {
        stop("findDuplicateFeatureIDs: Unsupported object type: ", class(theObject))
    }
}


# ----------------------------------------------------------------------------
# validateColumnMapping (Unified)
# ----------------------------------------------------------------------------
#' @title Validate Column Mapping for Omics Data
#' @description Checks that required columns exist in the data and returns validation status.
#'              Supports both lipidomics and metabolomics validation.
#'
#' @param data Data frame to validate.
#' @param id_column Name of the primary feature ID column (e.g., lipid_id or metabolite_id).
#' @param sample_columns Character vector of sample column names.
#' @param omics_type Character specifying "lipidomics" or "metabolomics". 
#'                   If NULL, attempts to guess based on context.
#'
#' @return A list containing validation status and summary statistics.
#' @keywords internal
#' @noRd
validateOmicsColumnMapping <- function(data, id_column, sample_columns, omics_type = NULL) {
    # If omics_type is not provided, we check which function to call.
    # We default to metabolomics if unknown, but prefer explicit type.
    if (!is.null(omics_type) && tolower(omics_type) == "lipidomics") {
        return(validateLipidColumnMapping(data, id_column, sample_columns))
    } else {
        # Default or explicit metabolomics
        return(validateMetabColumnMapping(data, id_column, sample_columns))
    }
}
# ----------------------------------------------------------------------------
# resolveDuplicateFeaturesByIntensity
# ----------------------------------------------------------------------------
#' Resolve Duplicate Features by Keeping Highest Average Intensity
#'
#' Within an assay tibble, identifies features with duplicate IDs and keeps only
#' the one with the highest average intensity across sample columns.
#'
#' @param assay_tibble A data frame or tibble representing one assay.
#' @param id_col Character string. The name of the column containing the feature IDs.
#' @param sample_cols Character vector. The names of the columns containing quantitative sample data.
#'
#' @return A tibble with duplicate features resolved based on highest average intensity.
#' @keywords internal
#' @importFrom dplyr group_by summarise ungroup filter slice_max select rowwise mutate c_across any_of
#' @importFrom rlang sym !!
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tibble column_to_rownames rownames_to_column
#' @export
resolveDuplicateFeaturesByIntensity <- function(assay_tibble, id_col, sample_cols) {
    if (!id_col %in% colnames(assay_tibble)) {
        warning(sprintf("ID column '%s' not found in assay tibble. Returning original tibble.", id_col))
        return(assay_tibble)
    }

    if (length(sample_cols) == 0) {
        warning("No sample columns provided. Returning original tibble.")
        return(assay_tibble)
    }

    # Check for duplicates first
    id_counts <- assay_tibble %>% dplyr::count(!!rlang::sym(id_col), name = "feature_count")
    duplicates_exist <- any(id_counts$feature_count > 1)

    if (!duplicates_exist) {
        # message(sprintf("No duplicates found in ID column '%s'. Returning original tibble.", id_col))
        return(assay_tibble)
    }

    message(sprintf("Resolving duplicates in ID column '%s' by keeping highest average intensity feature...", id_col))

    # Ensure sample columns are numeric for mean calculation
    assay_tibble_numeric <- assay_tibble %>%
        dplyr::mutate(dplyr::across(dplyr::any_of(sample_cols), as.numeric))

    # Calculate average intensity (handle NAs)
    # Using rowwise is more robust to non-numeric columns than converting to matrix first
    resolved_tibble <- assay_tibble_numeric %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            avg_intensity = mean(dplyr::c_across(dplyr::any_of(sample_cols)), na.rm = TRUE)
        ) %>%
        dplyr::ungroup() %>%
        # Handle cases where avg_intensity might be NaN (if all samples are NA)
        dplyr::mutate(avg_intensity = ifelse(is.nan(avg_intensity), -Inf, avg_intensity)) %>%
        # Group by the ID and keep the one with the highest average intensity
        dplyr::group_by(!!rlang::sym(id_col)) %>%
        # slice_max keeps ties by default; with_ties = FALSE ensures only one row per ID
        dplyr::slice_max(order_by = avg_intensity, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        # Remove the temporary average intensity column
        dplyr::select(-avg_intensity)

    # Report how many rows were removed
    rows_removed <- nrow(assay_tibble) - nrow(resolved_tibble)
    if (rows_removed > 0) {
        message(sprintf("Removed %d lower-intensity duplicate feature row(s).", rows_removed))
    }

    return(resolved_tibble)
}

# ============================================================================
# KIT-104 W2: the QC tab's computation for lipid, metab and trx, out of the modules (design 4.6, 5). Each function is
# the lipid and metab code at b423b764, which was the same code twice; where the two differed the difference is named.
# The shared QC sub-tabs (R/mod_omic_qc_tabs.R) call these with the omic's binding and its `omicGuiConfig(<omic>)$qc`.
# The module names and file:line notes below are provenance at b423b764: KIT-104 W5 and W6 deleted those functions and
# several of the files (mod_lipid_qc_intensity_server_helpers.R, mod_metab_qc_itsd.R and the metab helper files).
# ============================================================================

#' The intensity filter's per-assay table from the declaration's counts (lipid's inline table,
#' mod_lipid_qc_intensity_server_helpers.R:174-189; metab's `buildMetabIntensityFilterStats`, mod_metab_qc_intensity.R:173).
#' @param counts The declaration's `counts`: `before` and `after`, named by assay.
#' @return A data frame: `Assay`, `Original`, `Filtered`, `Removed`, `Percent_Retained` (rounded to one decimal).
#' @noRd
omicsQcCountsTable <- function(counts) {
    table <- data.frame(Assay = names(counts$before), Original = as.numeric(counts$before),
        Filtered = as.numeric(counts$after), Removed = as.numeric(counts$before) - as.numeric(counts$after),
        stringsAsFactors = FALSE)
    table$Percent_Retained <- round((table$Filtered / table$Original) * 100, 1)
    table
}

#' Duplicate feature ids in the head state, read from its id column only (KIT-082; lipid's
#' `handleLipidDuplicateDetection`, metab's `detectMetabDuplicateFeatures`).
#' @param binding The omic's declaration binding.
#' @param manager The session owner.
#' @return `list(duplicates, total)`: `findDuplicateFeatureIDs`' list per assay and the number of duplicated ids.
#' @noRd
detectWorkflowOmicsDuplicates <- function(binding, manager) {
    shiny::req(manager)
    ids <- readWorkflowOmicsIds(binding, manager)
    shiny::req(ids)
    if (!inherits(ids, binding$class)) stop(sprintf("Current state is not a %s object", binding$class))
    duplicates <- findDuplicateFeatureIDs(ids)
    total <- sum(vapply(duplicates, function(table) if (is.null(table)) 0L else nrow(table), integer(1)))
    list(duplicates = duplicates, total = total)
}

#' The duplicate resolution's per-assay statistics from the declaration's value: metab's product holds them
#' (`stats`), lipid's holds the row counts (`counts`), read per assay as `handleLipidDuplicateResolution` did.
#' @param value The declaration's value.
#' @param assay_names The state's assay names (lipid's product only).
#' @return A list per assay of `original`, `resolved` and `removed`.
#' @noRd
omicsDuplicateResolutionStats <- function(value, assay_names) {
    if (!is.null(value$stats)) return(value$stats)
    counts <- value$counts
    lapply(stats::setNames(assay_names, assay_names), function(assay) list(original = unname(counts$before[[assay]]),
        resolved = unname(counts$after[[assay]]), removed = unname(counts$before[[assay]]) - unname(counts$after[[assay]])))
}

#' The internal standards of the head state: the pattern given, else the state's own; per assay the first column of
#' the id, the annotation and `itsd$name_columns` in which the metrics kernel finds standards; their CVs and intensities.
#' Only the rows the analysis can use are read (KIT-082). Generalised from metab's `readMetabQcItsdState` and
#' `analyzeMetabQcItsdData` (mod_metab_qc_itsd.R:48, :57); lipid's inline analysis (mod_lipid_qc_itsd.R:218-380) gives
#' the same values, text and refusals. A NULL pattern is an empty one (metab's handling; lipid's errored).
#' @param binding The omic's declaration binding.
#' @param manager The session owner.
#' @param pattern The pattern typed in the tab; NULL or empty uses the state's `internal_standard_regex`.
#' @param itsd The omic's `qc$itsd`: `name_columns`, `metrics(assay_data, is_pattern, id_col, sample_id_col)`,
#'   `quant(assay_data)` and the refusal's `tip`.
#' @return `list(metrics, long_data, result_text, total, pattern, metrics_by_assay)`; `long_data` is NULL when no
#'   standard's intensities could be read.
#' @noRd
analyzeWorkflowOmicsItsd <- function(binding, manager, pattern, itsd) {
    shell <- readWorkflowOmicsShell(binding, manager)
    shiny::req(shell)
    if (!inherits(shell, binding$class)) stop(sprintf("Current state is not a %s object", binding$class))
    regex <- shell@internal_standard_regex
    pattern <- if (!is.null(pattern) && nzchar(pattern)) pattern else if (!is.na(regex) && nzchar(regex)) regex
    if (is.null(pattern)) stop("No internal standard pattern provided. Please enter a regex pattern to identify IS.")
    logger::log_info(paste("Analyzing internal standards with pattern:", pattern))
    id_column <- methods::slot(shell, binding$id_slot)
    columns <- unique(c(id_column, shell@annotation_id_column, itsd$name_columns))
    sample_id <- shell@sample_id
    rm(shell)
    assays <- methods::slot(readWorkflowOmicsRows(binding, manager, NULL,
        function(page, id_column) omicsInternalStandardRows(page, pattern, columns)), binding$data_slot)
    metrics_by_assay <- list()
    long_by_assay <- list()
    for (i in seq_along(assays)) {
        assay_data <- assays[[i]]
        assay_name <- if (!is.null(names(assays))) names(assays)[[i]] else paste0("Assay_", i)
        found <- data.frame(is_id = character(), mean_intensity = numeric(), cv = numeric())
        used <- NULL
        for (column in columns) {
            if (!is.null(column) && column %in% colnames(assay_data)) {
                found <- itsd$metrics(assay_data = assay_data, is_pattern = pattern, id_col = column, sample_id_col = sample_id)
                if (nrow(found) > 0) {
                    used <- column
                    logger::log_info(paste("Found IS in column:", column))
                    break
                }
            }
        }
        if (nrow(found) == 0) next
        found$assay <- assay_name
        metrics_by_assay[[assay_name]] <- found
        filter_column <- if (!is.null(used)) used else id_column
        sample_columns <- itsd$quant(assay_data)$sample_names
        if (length(sample_columns) == 0 || !filter_column %in% colnames(assay_data)) next
        rows <- assay_data[[filter_column]] %in% found$is_id
        if (!any(rows)) next
        long <- tidyr::pivot_longer(assay_data[rows, c(filter_column, sample_columns), drop = FALSE],
            cols = dplyr::all_of(sample_columns), names_to = "Sample", values_to = "Intensity")
        long$assay <- assay_name
        names(long)[1] <- "IS_ID"
        long_by_assay[[assay_name]] <- long
    }
    if (length(metrics_by_assay) == 0) {
        first_columns <- if (length(assays) > 0) colnames(assays[[1]]) else character()
        stop(paste0("No internal standards found. Pattern: ", pattern,
            "\nColumns searched: ", paste(columns[columns %in% first_columns], collapse = ", "), "\n", itsd$tip))
    }
    metrics <- do.call(rbind, metrics_by_assay)
    total <- nrow(metrics)
    result_text <- paste("Internal Standard Analysis Complete", "====================================",
        sprintf("Pattern used: %s", pattern), sprintf("Total IS detected: %d", total), "", "CV Statistics:",
        sprintf("  Median CV: %.1f%%", stats::median(metrics$cv, na.rm = TRUE)),
        sprintf("  Max CV: %.1f%%", max(metrics$cv, na.rm = TRUE)),
        sprintf("  IS with CV > 30%%: %d", sum(metrics$cv > 30, na.rm = TRUE)), "", "Per-Assay IS Counts:",
        paste(vapply(names(metrics_by_assay), function(name) {
            sprintf("  %s: %d internal standards", name, nrow(metrics_by_assay[[name]]))
        }, character(1)), collapse = "\n"), sep = "\n")
    list(metrics = metrics, long_data = if (length(long_by_assay) > 0) do.call(rbind, long_by_assay),
        result_text = result_text, total = total, pattern = pattern, metrics_by_assay = metrics_by_assay)
}

#' The ITSD CV lollipop figure (metab's `buildMetabQcItsdCvPlot`, mod_metab_qc_itsd.R:279; lipid's inline copy,
#' mod_lipid_qc_itsd.R:155-201).
#' @param metrics The analysis's `metrics`.
#' @return A ggplot.
#' @noRd
omicsItsdCvPlot <- function(metrics) {
    metrics$is_id <- stats::reorder(metrics$is_id, metrics$cv)
    metrics$cv_status <- dplyr::case_when(
        metrics$cv <= 15 ~ "Good (<15%)",
        metrics$cv <= 30 ~ "Acceptable (15-30%)",
        TRUE ~ "Review (>30%)"
    )
    metrics$cv_status <- factor(metrics$cv_status, levels = c("Good (<15%)", "Acceptable (15-30%)", "Review (>30%)"))
    ggplot2::ggplot(metrics, ggplot2::aes(x = is_id, y = cv)) +
        ggplot2::geom_hline(yintercept = 15, linetype = "dashed", color = "#2ecc71", linewidth = 0.8, alpha = 0.7) +
        ggplot2::geom_hline(yintercept = 30, linetype = "dashed", color = "#e67e22", linewidth = 0.8, alpha = 0.7) +
        ggplot2::geom_segment(ggplot2::aes(xend = is_id, yend = 0), color = "grey60", linewidth = 0.8) +
        ggplot2::geom_point(ggplot2::aes(color = cv_status), size = 4) +
        ggplot2::facet_wrap(~assay, scales = "free_y") +
        ggplot2::scale_color_manual(
            values = c("Good (<15%)" = "#2ecc71", "Acceptable (15-30%)" = "#e67e22", "Review (>30%)" = "#e74c3c"),
            name = "CV Status"
        ) +
        ggplot2::labs(title = "Internal Standard CV", subtitle = "Dashed lines: 15% (green) and 30% (orange) thresholds",
            x = NULL, y = "CV (%)") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
            legend.position = "bottom",
            panel.grid.major.x = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(face = "bold")
        ) +
        ggplot2::coord_flip()
}

#' The ITSD intensity-trend figure (metab's `buildMetabQcItsdIntensityPlot`, mod_metab_qc_itsd.R:340; lipid's inline
#' copy, mod_lipid_qc_itsd.R:404-422).
#' @param long The analysis's `long_data`.
#' @return A ggplot.
#' @noRd
omicsItsdIntensityPlot <- function(long) {
    ggplot2::ggplot(long, ggplot2::aes(x = Sample, y = Intensity, color = IS_ID, group = IS_ID)) +
        ggplot2::geom_line(alpha = 0.7) +
        ggplot2::geom_point(size = 1, alpha = 0.5) +
        ggplot2::facet_wrap(~assay, scales = "free_y") +
        ggplot2::labs(title = "Internal Standard Intensity Across Samples", x = "Sample", y = "Intensity",
            color = "Internal Standard") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
            legend.position = "bottom")
}

#' The QC finalise's checks on the head state's shell (lipid's `validateLipidQcS4FinalizeState`,
#' mod_lipid_qc_s4.R:347-404; metab's `validateMetabQcS4FinalizeState`; trx's checks in `trxQcFinalizeCollaborators`):
#' the class; the assays (any zero-row assay refused when `finalize_zero_rows`, otherwise only an empty list, D2); the
#' design rows; then the omic's alignment (`finalize_alignment`, D10/D11), for every omic.
#' @param shell The head state's shell (`readWorkflowOmicsShell`, `measures = "features"`).
#' @param config The omic's `omicGuiConfig(<omic>)$qc`: `finalize_zero_rows`, `texts` and `finalize_alignment`.
#' @param binding The omic's declaration binding (the class and its assay slot).
#' @return The shell, when it passes. A refusal is an error with today's text.
#' @noRd
validateWorkflowOmicsQcFinalize <- function(shell, config, binding) {
    shiny::req(shell)
    if (!inherits(shell, binding$class)) stop(sprintf("Current state is not a %s object", binding$class))
    slots <- methods::slotNames(shell)
    assays <- if (binding$data_slot %in% slots) methods::slot(shell, binding$data_slot)
    rows <- attr(shell, "assay_rows") %||% vapply(assays, nrow, integer(1))
    if (length(assays) == 0L || (isTRUE(config$finalize_zero_rows) && any(rows == 0L))) {
        stop(config$texts[["finalize_no_assays"]])
    }
    if (!"design_matrix" %in% slots || !is.data.frame(shell@design_matrix) || nrow(shell@design_matrix) == 0L) {
        stop(config$texts[["finalize_no_design"]])
    }
    # Every omic's finalise checks its alignment; a configuration without one is refused, never skipped.
    if (!is.function(config$finalize_alignment)) {
        cli::cli_abort("the qc configuration of {.val {binding$omic_type}} has no {.field finalize_alignment}.",
            class = "omic_gui_config_invalid")
    }
    alignment <- config$finalize_alignment(shell)
    if (!isTRUE(alignment$valid)) {
        stop(sprintf("Current %s object has assay/design mismatch: %s", binding$class, paste(alignment$errors, collapse = "; ")))
    }
    shell
}

#' D11. Lipid's and metab's finalise alignment: the shared design alignment with the omic's rules and lipid's arguments
#' at mod_lipid_qc_s4.R:374-394 (the design's runs as the sample columns, no replicate requirement), on the shell's
#' zero-row prototypes.
#' @param shell The head state's shell.
#' @param omic_type `lipidomics` or `metabolomics`.
#' @return `list(valid, errors)`.
#' @noRd
omicsQcFinalizeAlignment <- function(shell, omic_type) {
    binding <- omicsDeclarationBindings()[[omic_type]]
    rules <- omicsDesignRules(omic_type)
    design <- shell@design_matrix
    runs <- if (shell@sample_id %in% names(design)) as.character(design[[shell@sample_id]]) else character()
    mapping <- c(stats::setNames(list(methods::slot(shell, binding$id_slot)), rules$id_key),
        list(annotation_col = shell@annotation_id_column, sample_columns = runs))
    alignment <- validateOmicsDesignAlignment(design, methods::slot(shell, binding$data_slot), rules,
        columnMapping = mapping, sampleIdCol = shell@sample_id, groupCol = shell@group_id, requireReplicates = FALSE)
    alignment[c("valid", "errors")]
}
