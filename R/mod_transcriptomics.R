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

# KIT-096: the transcriptomics workflow module (import, design, QC). The tabs are the shared import (KIT-102) and the
# lipid design-builder and QC modules under the transcriptomics profile (`trxOmicProfile()`), which supplies the omic's collaborators at the seams those
# modules already take; nothing is copied. The tabs hold state references only (KIT-075) and run on the managed runner, or
# the binding's wrappers in process (KIT-078). Recorded divergence from conventions/R-shiny.md: the shinydashboard and
# wellPanel shell, the managed runner in place of ExtendedTask and the import's fileInput fallback are kept, to match the
# other omics (the user's instruction) under the existing-code scope rule of conventions/R.md.

# ---------------------------------------------------------------------------------------------------------------
# The profile.
# ---------------------------------------------------------------------------------------------------------------

#' The transcriptomics profile of the lipid design, QC, normalisation, DA and summary modules: labels, test-id prefix,
#' the binding's readers and wrappers, and the collaborators each module takes at its seams. The import is the shared
#' feature-table panel under `omicGuiConfig("transcriptomics")` (KIT-102), which the profile's `omic_type` selects.
#' @return A plain list.
#' @noRd
trxOmicProfile <- function() {
    profile <- list(omic_type = "transcriptomics", prefix = "trx", testid = "trx", label = "Transcriptomics", label_lower = "transcriptomics",
        class = "TranscriptomicsAssayData", workflow_type = "TRANSCRIPTOMICS",
        design_help = paste("Samples are matched to the count columns by sample name only, in the counts file's order; there",
            "is no positional mapping. Each group needs at least two samples. The order in which you add the factor levels is",
            "the group level order the analysis uses."),
        run_order = identity, group_levels = trxDesignGroupLevels, validate_design = trxDesignCheck, build_object = trxDesignObject,
        record_provenance = recordTrxDesignProvenance, state_name = "trx_raw_data_s4",
        state_description = "Initial TranscriptomicsAssayData S4 object created after design matrix",
        qc_title = "Transcriptomics Quality Control & Filtering", qc_tabs = trxQcTabsUi, qc_init = initializeTrxQcSubmodules,
        read_metadata = function(manager) readWorkflowOmicsMetadata(transcriptomicsDeclarationBinding, manager),
        read_state = readTrxWorkflowState, register_qc_summary = registerTrxQcSummaryOutput, finalize = trxQcFinalizeCollaborators(),
        binding = transcriptomicsDeclarationBinding, report_template = "transcriptomics_report.rmd", report_dir = "transcriptomics")
    c(profile, trxNormDaSummaryProfile())
}

#' The normalisation, DA and summary entries of the transcriptomics profile (KIT-097).
#' @noRd
trxNormDaSummaryProfile <- function() {
    list(norm_options = trxNormOptionsPanel, norm_server = mod_trx_norm_server, norm_tabs = c("ruv", "pca", "rle", "log"),
        norm_stage_labels = c("Raw logCPM", "TMM logCPM", "RUVs-corrected logCPM (display only)"),
        norm_ruv_heading = "RUVs Control Selection (canonical correlation)",
        norm_ruv_help = paste("The CCA panels and the grid of percentages: the controls are inferred negative controls (the least",
            "significant genes in a preliminary group GLM)."),
        da_server = mod_trx_da_server, da_model_ui = trxDaModelUi, da_threshold_label = "FDR (BH) threshold:",
        da_lfc_label = "Fold-change line (display only; no test):", long_table = "da_genes_long",
        da_columns = list(id = "GeneID", label = "GeneName", pvalue = "PValue", significance = "FDR", significance_label = "FDR (BH)"),
        da_static_volcano = trxDaStaticVolcano, da_volcano_interactive = FALSE, da_widget = trxDaWidget, da_with_source = trxDaWithSource,
        da_formula = trxRealisedFormula, da_contrasts = trxDaContrastEntries, da_export = trxDaExport, da_run_reference = runTrxReferenceDaAnalysis,
        summary_export_ui = trxSummaryExportUi, summary_server = trxSummaryServer, summary_context = trxSummaryContext,
        summary_state_object = trxSummaryStateObject, summary_placeholder = "e.g., my_transcriptomics_analysis",
        summary_description = "Full transcriptomics analysis workflow with normalization, DE and the RUVIII integration export",
        ruviii_m_choices = trxRuviiiMChoices, integration_filename = "transcriptomics_normalized_data.tsv",
        integration_write = trxSummaryIntegrationWrite, integration_ready = trxSummaryIntegrationReady,
        integration_missing_message = "The integration export needs the RUVIII export: run it in the panel above, then save again.")
}

# ---------------------------------------------------------------------------------------------------------------
# Design: the level order, KIT-093's validation at save, the design state.
# ---------------------------------------------------------------------------------------------------------------

#' The group level order the user entered: the assigned groups ordered by the positions of their factor levels in the
#' builder's factor list (the order they were added); a group of several factors by each factor in turn.
#' @noRd
trxDesignGroupLevels <- function(design, factors) {
    assigned <- !is.na(design$group) & design$group != ""
    groups <- unique(design$group[assigned])
    first <- design[match(groups, design$group), , drop = FALSE]
    position <- function(column) if (column %in% names(first)) match(first[[column]], factors) else rep(NA_integer_, length(groups))
    groups[order(position("factor1"), position("factor2"), position("factor3"), seq_along(groups))]
}

#' A refusal's text as the tab shows it: the message and its class.
#' @noRd
trxConditionText <- function(error) {
    text <- conditionMessage(error)
    class <- class(error)[[1L]]
    if (grepl(class, text, fixed = TRUE)) text else sprintf("%s (%s)", text, class)
}

#' KIT-093's design validation at save, by name only, against every imported count column except the runs the user
#' removed (an unassigned run is a count column outside the design), with at least two runs per group (RUVs and RUVIII
#' need replicates). Returns the lipid preflight's shape.
#' @noRd
trxDesignCheck <- function(results, columnMapping) {
    runs <- setdiff(columnMapping$sample_columns, results$excluded_runs)
    tryCatch({
        levels <- results$group_levels %||% trxDesignGroupLevels(results$design_matrix, character())
        checked <- validateTrxDesign(runs, results$design_matrix, levels)
        sizes <- table(checked$group)
        small <- names(sizes)[sizes < 2L]
        if (length(small)) {
            trxAbort(sprintf("group(s) with fewer than two samples: %s", paste(small, collapse = ", ")),
                c("trx_design_group_singleton", "trx_design_invalid", "trx_error"))
        }
        list(valid = TRUE, errors = character())
    }, trx_error = function(error) list(valid = FALSE, errors = trxConditionText(error)))
}

#' The design state's object: the builder's count table and design, the recorded level order, the chosen id and name
#' columns, and the builder's `tech_rep_group` as the technical replicates.
#' @noRd
trxDesignObject <- function(results, columnMapping) {
    createTranscriptomicsAssayData(results$data_cln, results$design_matrix, results$group_levels,
        roles = list(gene_id_column = columnMapping$lipid_id_col, annotation_id_column = columnMapping$annotation_col,
            technical_replicate_id = "tech_rep_group"))
}

# ---------------------------------------------------------------------------------------------------------------
# The module.
# ---------------------------------------------------------------------------------------------------------------

#' @title Transcriptomics Workflow Module
#' @description The transcriptomics workflow's tabs: import, design and QC (KIT-096), normalisation, DA and the summary
#'   (KIT-097). KIT-102: the shared shell (mod_omic_workflow.R) under the transcriptomics configuration; its tabs are the
#'   lipid modules under `trxOmicProfile()` (the component registry, `omicGuiComponents`).
#' @param id Module id (`transcriptomics_workflow` in the app).
#' @name mod_transcriptomics
#' @noRd
mod_transcriptomics_ui <- function(id) {
    mod_omic_workflow_ui(id, omicGuiConfig("transcriptomics"))
}

#' @rdname mod_transcriptomics
#' @param project_dirs,omic_type,experiment_label,volumes As `mod_lipidomics_server`.
#' @return The module's bus.
#' @noRd
mod_transcriptomics_server <- function(id, project_dirs, omic_type, experiment_label, volumes = NULL) {
    mod_omic_workflow_server(id, project_dirs, omic_type, experiment_label, volumes, omicGuiConfig("transcriptomics"))
}
