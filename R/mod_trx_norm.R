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

# KIT-097: the transcriptomics normalisation tab, the lipid norm module under the transcriptomics profile. The panels are
# the lipid tab's (RUV QC, PCA, RLE, Log; one assay row) with its renderers; the options and the sequence are the
# transcriptomics ones: TMM and logCPM (`trx_tmm_logcpm`), the RUV grid (`trx_ruv_grid`, automatic as the notebook or
# manual), RUVs (`trx_ruvs`, the DA head), then the PCA and RLE at raw, TMM and RUVs, the CCA panels and the stage grid.
# `runWorkflowOmicsManagedNormalisation` (ITSD, log2, between-sample, RUV-III) is not this sequence and is not reused.
# KIT-110 W2: the sequence, the reset and the in-process publications moved to R/utils_workflow_trx_norm.R (the core)
# and R/utils_workflow_omics_actions.R (`runWorkflowOmicsPublicationInProcess`); the options panel and the server stay
# here until W3.

#' The normalisation options of the transcriptomics profile (the lipid panel's layout and styling).
#' @noRd
trxNormOptionsPanel <- function(ns) {
    automatic <- shiny::conditionalPanel(condition = "input.ruv_mode == 'automatic'", ns = ns,
        shiny::h5("Automatic Optimization Parameters"),
        shiny::fluidRow(shiny::column(4, shiny::numericInput(ns("auto_percentage_min"), "Min %:", value = 5, min = 1, max = 50, width = "100%")),
            shiny::column(4, shiny::numericInput(ns("auto_percentage_max"), "Max %:", value = 30, min = 1, max = 50, width = "100%")),
            shiny::column(4, shiny::numericInput(ns("auto_percentage_step"), "Step %:", value = 5, min = 1, max = 25, width = "100%"))),
        shiny::numericInput(ns("k_penalty_weight"), "K Value Penalty Weight:", value = 0.05, min = 0, max = 0.9, step = 0.01, width = "100%"),
        shiny::helpText("The notebook's grid: 5 to 30 % by 5, k capped at 3 below 15 samples and 4 otherwise, penalty 0.05; the",
            "percentage and k with the maximum composite score are chosen."))
    manual <- shiny::conditionalPanel(condition = "input.ruv_mode == 'manual'", ns = ns,
        shiny::sliderInput(ns("ruv_percentage"), "% Genes as Inferred Negative Controls:", min = 1, max = 50, value = 10, step = 1, width = "100%"),
        shiny::numericInput(ns("ruv_k"), "Number of Factors (k):", value = 2, min = 1, max = 4, width = "100%"),
        shiny::helpText("k must not exceed the cap (3 below 15 samples, 4 otherwise); a larger k is refused."))
    shiny::column(3, shiny::wellPanel(
        shiny::h4("Normalization Options"), shiny::hr(),
        shiny::h4("TMM and logCPM"),
        shiny::helpText("edgeR TMM normalisation factors of the filtered counts; log2 counts per million with a prior count of 2."),
        shiny::hr(), shiny::h4("RUVs (RUVSeq)"),
        testid(shiny::radioButtons(ns("ruv_mode"), "RUV Parameter Tuning:", choices = list("Automatic (recommended)" = "automatic",
            "Manual (advanced users)" = "manual"), selected = "automatic", width = "100%"), "trx-norm-ruv-mode"),
        shiny::helpText("RUVs uses inferred negative controls (the least significant genes in a preliminary group GLM). Its factors",
            "W_1..W_k enter the DA model as covariates; the RUVs-corrected logCPM is for display only."),
        automatic, manual, shiny::hr(),
        testid(shiny::actionButton(ns("run_normalization"), "Run Normalization Pipeline", class = "btn-primary", width = "100%",
            icon = shiny::icon("play")), "trx-norm-run"),
        shiny::br(), shiny::br(),
        testid(shiny::actionButton(ns("reset_normalization"), "Reset to Pre-Normalization", class = "btn-warning", width = "100%",
            icon = shiny::icon("undo")), "trx-norm-reset")))
}


#' @title Transcriptomics Normalisation Tab
#' @description The lipid norm module's server under the transcriptomics profile (`mod_lipid_norm_server(profile =)`).
#' @param id,workflow_data,experiment_paths,omic_type,experiment_label,selected_tab As `mod_lipid_norm_server`.
#' @return The tab's reactive values.
#' @name mod_trx_norm
#' @noRd
mod_trx_norm_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        normData <- omicNormReactiveState()
        # One count table; its figures are the worker's trx_<stage>_<kind>.png in the gene QC directory.
        normData$assay_names <- "counts"
        addLog <- omicNormAddLog(normData)
        registerOmicNormQcImageOutputs(output, omicNormQcImageRenderer(normData, experiment_paths,
            sanitizeAssayNameFn = function(assay_name) "trx", qcDirKey = "gene_qc_dir"),
            stagePrefixes = c(pre_norm = "raw", post_norm = "tmm", ruv_corrected = "ruvs"))
        registerOmicNormAssayLabelOutputs(output, omicNormAssayLabelRenderer(normData))
        registerOmicNormLogOutput(output, omicNormLogRenderer(normData))
        registerOmicNormRuvQcOutput(output, omicNormRuvQcUiRenderer(normData, session$ns))
        registerOmicNormRuvCancorOutputs(output, normData)
        # KIT-110 W2: the run and reset are the core's handlers (utils_workflow_trx_norm.R), with the tab's notification and
        # status callbacks; the handlers are looked up by name when the button is pressed.
        shiny::observeEvent(input$run_normalization, {
            shiny::req(workflow_data$state_manager)
            runTrxNormalisation(list(managed = runTrxNormManaged, in_process = runTrxNormInProcess), transcriptomicsDeclarationBinding,
                workflow_data$state_manager, trxNormGridParameters(input), normData, workflow_data$config_list, experiment_paths,
                addLog, notify = function(...) shiny::showNotification(...), setStatus = function(tab, value) {
                    status <- workflow_data$tab_status
                    status[[tab]] <- value
                    workflow_data$tab_status <- status
                })
        })
        shiny::observeEvent(input$reset_normalization, resetTrxNormalisation(workflow_data$state_manager, normData,
            workflow_data$config_list, addLog, notify = function(...) shiny::showNotification(...), setStatus = function(tab, value) {
                status <- workflow_data$tab_status
                status[[tab]] <- value
                workflow_data$tab_status <- status
            }))
        normData
    })
}
