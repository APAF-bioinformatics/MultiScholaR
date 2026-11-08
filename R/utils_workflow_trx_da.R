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

# KIT-106 D2: the transcriptomics DA tab's core functions, moved unchanged with their names from R/mod_trx_da.R (at
# 2269a45f), which W3 deletes: the model text, the RUVs head's metadata, the contrast entries the binding reads, the static
# volcano and the export's report. The DA tab's configuration names them (`da$model$text`, `da$source$read_metadata`,
# `da$contrast_entries`, `da$contrast_entry`, `da$views$static_volcano$plot`, `da$run$export_report`).

#' The realised formula and model text of the head (metadata only).
#' @noRd
trxDaModelText <- function(workflowData) {
    manager <- workflowData$state_manager
    formula <- tryCatch(withWorkflowStateScope(manager, function() trxRealisedFormula(manager, getWorkflowStateReference(manager))),
        error = function(error) NULL)
    if (is.null(formula) || !grepl("W_", formula, fixed = TRUE)) return("Run the normalisation (RUVs) first; the formula is realised from its design.")
    paste(c(paste("Formula:", formula), "Fit: glmQLFit(robust = TRUE)", "Test: QL F-test (glmQLFTest)", "Adjustment: BH (FDR)"), collapse = "\n")
}

#' The head's metadata when it is a RUVs state (its design has W columns), the state the DA runs on; NULL otherwise.
#' @noRd
trxDaReadMetadata <- function(workflowData) {
    meta <- tryCatch(readWorkflowOmicsMetadata(transcriptomicsDeclarationBinding, workflowData$state_manager), error = function(error) NULL)
    if (is.null(meta) || !any(startsWith(names(meta$design_matrix), "W_"))) return(NULL)
    meta
}

#' The contrasts the binding reads: the builder's `name=expression` entries (its `full_format`).
#' @noRd
trxDaContrastEntries <- function(contrastsTbl) {
    entries <- if ("full_format" %in% names(contrastsTbl)) contrastsTbl$full_format else contrastsTbl$contrasts
    trxContrastsTable(as.character(entries))
}

#' The contrast entry of a contrast name (the tab's selectors show names; the results are keyed by entry).
#' @noRd
trxDaContrastEntry <- function(workflowData, name) {
    entries <- trxDaContrastEntries(workflowData$contrasts_tbl)
    hit <- entries$contrasts[entries$friendly_names == name]
    if (length(hit)) hit[[1L]] else name
}

#' The static volcano of one contrast: KIT-093's renderer (the notebook's): y is -log10 FDR (BH), significance FDR (BH) <
#' threshold; the +-1 logFC lines are decorative.
#' @noRd
trxDaStaticVolcano <- function(da_results_list, selected_contrast, selected_assay = NULL, da_q_val_thresh = 0.05, lfc_threshold = 1, ...) {
    shiny::req(da_results_list, selected_contrast)
    long <- da_results_list$da_genes_long
    trxVolcano(long[long$comparison == selected_contrast, , drop = FALSE], selected_contrast, fdr_thresh = da_q_val_thresh)$plot
}

#' The worker's DA publication in the words of the tab's report (`success`).
#' @noRd
trxDaExport <- function(export) {
    export$success <- isTRUE(export$saved)
    export
}
