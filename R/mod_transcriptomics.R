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
#' feature-table panel under `omicGuiConfig("transcriptomics")` (KIT-102), which the profile's `omic_type` selects. So is
#' the design tab (KIT-103 W5): its eight design fields left the profile for `omicGuiConfig("transcriptomics")$design`,
#' and its design functions (`trxDesignGroupLevels`, `trxDesignCheck`, `trxDesignObject`) are in func_trx_design.R.
#' So is the QC tab (KIT-104 W6): its QC fields left the profile for `omicGuiConfig("transcriptomics")$qc`, the shared
#' shell and finalise (mod_omic_qc.R, mod_omic_qc_tabs.R); KIT-104 W8 removed `class` and `read_state`, which nothing
#' else read (the configuration's binding carries both).
#' @return A plain list.
#' @noRd
trxOmicProfile <- function() {
    profile <- list(omic_type = "transcriptomics", prefix = "trx", testid = "trx", label = "Transcriptomics", label_lower = "transcriptomics",
        workflow_type = "TRANSCRIPTOMICS",
        binding = transcriptomicsDeclarationBinding, report_template = "transcriptomics_report.rmd", report_dir = "transcriptomics")
    c(profile, trxNormDaSummaryProfile())
}

#' The summary entries of the transcriptomics profile (KIT-097). KIT-110 W3 removed the six normalisation entries: their
#' values are `omicGuiConfig("transcriptomics")$norm` (`omicGuiTrxNorm`). KIT-106 W3 removed the 14 DA entries: their
#' values are `omicGuiConfig("transcriptomics")$da` (`omicGuiTrxDa`), the long table is the binding's, and the tab is the
#' shared DA frame (mod_omic_da.R).
#' @noRd
trxNormDaSummaryProfile <- function() {
    list(summary_export_ui = trxSummaryExportUi, summary_server = trxSummaryServer, summary_context = trxSummaryContext,
        summary_state_object = trxSummaryStateObject, summary_placeholder = "e.g., my_transcriptomics_analysis",
        summary_description = "Full transcriptomics analysis workflow with normalization, DE and the RUVIII integration export",
        ruviii_m_choices = trxRuviiiMChoices, integration_filename = "transcriptomics_normalized_data.tsv",
        integration_write = trxSummaryIntegrationWrite, integration_ready = trxSummaryIntegrationReady,
        integration_missing_message = "The integration export needs the RUVIII export: run it in the panel above, then save again.")
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
