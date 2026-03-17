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
# mod_metab_norm.R
# ============================================================================
# Purpose: Metabolomics normalization orchestrator module (Proteomics UX)
#
# This module provides a unified normalization pipeline with per-assay
# visualization for metabolomics data. It orchestrates sub-modules for
# imputation, ITSD selection, and scaling/normalization.
# ============================================================================

#' @title Metabolomics Normalization Module
#' @description A Shiny module for multi-step normalization of metabolomics data.
#'              Harmonized with proteomics UX pattern.
#' @name mod_metab_norm
NULL

#' @rdname mod_metab_norm
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h4 h5 p hr br radioButtons selectInput textInput actionButton uiOutput verbatimTextOutput icon tags checkboxInput numericInput plotOutput conditionalPanel tabsetPanel tabPanel sliderInput helpText
#' @importFrom shinyjqui jqui_resizable
#' @importFrom DT dataTableOutput
mod_metab_norm_ui <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::wellPanel(
            shiny::fluidRow(
                # ================================================================
                # LEFT PANEL (width = 3): All normalization controls
                # ================================================================
                shiny::column(3
                    , shiny::wellPanel(
                        shiny::h4("Normalization Options")

                        # --- Plot Aesthetics ---
                        , shiny::hr()
                        , shiny::h5("Plot Aesthetics")
                        , shiny::selectInput(
                            ns("color_variable")
                            , "Color by:"
                            , choices = c("group", "factor1", "factor2", "batch")
                            , selected = "group"
                            , width = "100%"
                        )
                        , shiny::selectInput(
                            ns("shape_variable")
                            , "Shape by:"
                            , choices = c("group", "factor1", "factor2", "batch")
                            , selected = "group"
                            , width = "100%"
                        )

                        # --- ITSD Normalization ---
                        , shiny::hr()
                        , shiny::h4("ITSD Normalization")
                        , shiny::checkboxInput(
                            ns("apply_itsd")
                            , "Apply ITSD normalization"
                            , value = TRUE
                            , width = "100%"
                        )
                        , shiny::conditionalPanel(
                            condition = "input.apply_itsd == true"
                            , ns = ns
                            , shiny::selectInput(
                                ns("itsd_aggregation")
                                , "ITSD Aggregation:"
                                , choices = c("Median (recommended)" = "median", "Mean" = "mean", "Sum" = "sum")
                                , selected = "median"
                                , width = "100%"
                            )
                        )

                        # --- Log2 Transform ---
                        , shiny::hr()
                        , shiny::h4("Log2 Transform")
                        , shiny::numericInput(
                            ns("log_offset")
                            , "Offset for zeros:"
                            , value = 1
                            , min = 0
                            , step = 0.1
                            , width = "100%"
                        )

                        # --- Normalization Method ---
                        , shiny::hr()
                        , shiny::h4("Between-Sample Normalization")
                        , shiny::selectInput(
                            ns("norm_method")
                            , "Normalization Method:"
                            , choices = list(
                                "Cyclic Loess" = "cyclicloess"
                                , "Quantile" = "quantile"
                                , "Scale (Median Absolute Values)" = "scale"
                                , "None" = "none"
                            )
                            , selected = "cyclicloess"
                            , width = "100%"
                        )

                        # --- RUV-III Options ---
                        , shiny::hr()
                        , shiny::h4("RUV-III Batch Correction")
                        , shiny::radioButtons(
                            ns("ruv_mode")
                            , "RUV Parameter Tuning:"
                            , choices = list(
                                "Automatic (recommended)" = "automatic"
                                , "Manual (advanced users)" = "manual"
                                , "Skip RUV" = "skip"
                            )
                            , selected = "automatic"
                            , width = "100%"
                        )

                        , shiny::selectInput(
                            ns("ruv_grouping_variable")
                            , "RUV Grouping Variable:"
                            , choices = c("group", "factor1", "factor2", "batch")
                            , selected = "group"
                            , width = "100%"
                        )

                        # --- Automatic RUV Parameters ---
                        , shiny::conditionalPanel(
                            condition = "input.ruv_mode == 'automatic'"
                            , ns = ns
                            , shiny::h5("Automatic Optimization Parameters")
                            , shiny::fluidRow(
                                shiny::column(6, shiny::numericInput(ns("auto_percentage_min"), "Min %:", value = 1, min = 1, max = 50))
                                , shiny::column(6, shiny::numericInput(ns("auto_percentage_max"), "Max %:", value = 20, min = 1, max = 50))
                            )
                            , shiny::selectInput(ns("separation_metric"), "Metric:", 
                                               choices = list("Max Diff" = "max_difference", "Mean Diff" = "mean_difference", "AUC" = "auc"))
                            , shiny::sliderInput(ns("k_penalty_weight"), "K Penalty:", min = 0.1, max = 0.9, value = 0.5, step = 0.1)
                            , shiny::checkboxInput(ns("adaptive_k_penalty"), "Adaptive K Penalty", value = TRUE)
                        )

                        # --- Manual RUV Parameters ---
                        , shiny::conditionalPanel(
                            condition = "input.ruv_mode == 'manual'"
                            , ns = ns
                            , shiny::sliderInput(ns("ruv_percentage"), "% Neg Ctrl:", min = 1, max = 20, value = 5)
                            , shiny::numericInput(ns("ruv_k"), "Factors (k):", value = 2, min = 1, max = 10)
                        )

                        # --- Action Buttons ---
                        , shiny::hr()
                        , shiny::actionButton(
                            ns("run_normalization")
                            , "Run Normalization Pipeline"
                            , class = "btn-primary"
                            , width = "100%"
                            , icon = shiny::icon("play")
                        )
                        , shiny::br(), shiny::br()
                        , shiny::actionButton(
                            ns("reset_normalization")
                            , "Reset to Pre-Normalization"
                            , class = "btn-warning"
                            , width = "100%"
                            , icon = shiny::icon("undo")
                        )
                    )
                )

                # ================================================================
                # RIGHT PANEL (width = 9): Sub-modules in tabs
                # ================================================================
                , shiny::column(9
                    , shiny::tabsetPanel(
                        id = ns("norm_qc_tabs")
                        , mod_metab_norm_impute_ui(ns("impute"))
                        , mod_metab_norm_itsd_ui(ns("itsd"))
                        , shiny::tabPanel("RUV QC", icon = shiny::icon("chart-line"), shiny::uiOutput(ns("ruv_qc_ui")))
                        , mod_metab_norm_scale_ui(ns("pca"), "pca")
                        , mod_metab_norm_scale_ui(ns("density"), "density")
                        , mod_metab_norm_scale_ui(ns("rle"), "rle")
                        , mod_metab_norm_scale_ui(ns("correlation"), "correlation")
                        , shiny::tabPanel("Log", icon = shiny::icon("list-alt"), shiny::verbatimTextOutput(ns("norm_log")))
                    )
                )
            )
        )
    )
}

#' @rdname mod_metab_norm
#' @export
mod_metab_norm_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Reactive State
        norm_data <- shiny::reactiveValues(
            assay_names = NULL
            , plot_refresh_trigger = 0
            , normalization_log = character(0)
            , ruv_results = list()
        )

        # Initialize Sub-modules
        mod_metab_norm_impute_server("impute", workflow_data, omic_type, experiment_label)
        itsd_module <- mod_metab_norm_itsd_server("itsd", workflow_data, omic_type, experiment_label)
        mod_metab_norm_scale_server("pca", workflow_data, experiment_paths, norm_data, "pca")
        mod_metab_norm_scale_server("density", workflow_data, experiment_paths, norm_data, "density")
        mod_metab_norm_scale_server("rle", workflow_data, experiment_paths, norm_data, "rle")
        mod_metab_norm_scale_server("correlation", workflow_data, experiment_paths, norm_data, "correlation")

        # Logging helper
        add_log <- function(message) {
            timestamp <- format(Sys.time(), "%H:%M:%S")
            norm_data$normalization_log <- c(norm_data$normalization_log, sprintf("[%s] %s", timestamp, message))
        }

        # Initialize Assay Names
        shiny::observe({
            shiny::req(workflow_data$state_manager)
            current_s4 <- workflow_data$state_manager$getState()
            if (!is.null(current_s4) && inherits(current_s4, "MetaboliteAssayData")) {
                norm_data$assay_names <- names(current_s4@metabolite_data)
            }
        })

        # Run Normalization Logic
        shiny::observeEvent(input$run_normalization, {
            shiny::req(workflow_data$state_manager)
            add_log("Starting normalization pipeline...")
            
            shiny::withProgress(message = "Processing...", value = 0, {
                tryCatch({
                    current_s4 <- workflow_data$state_manager$getState()
                    
                    # 1. ITSD
                    if (isTRUE(input$apply_itsd)) {
                        add_log("Applying ITSD...")
                        ids <- itsd_module$get_feature_ids(current_s4)
                        current_s4 <- normaliseUntransformedData(current_s4, method = "ITSD", 
                                                               itsd_aggregation = input$itsd_aggregation, 
                                                               itsd_feature_ids = ids)
                    }
                    
                    # 2. Log2
                    add_log("Applying Log2...")
                    current_s4 <- logTransformAssays(current_s4, offset = input$log_offset)
                    
                    # 3. Between-sample
                    if (input$norm_method != "none") {
                        add_log(paste("Applying", input$norm_method))
                        current_s4 <- normaliseBetweenSamples(current_s4, normalisation_method = input$norm_method)
                    }
                    
                    # Capture Checkpoint cp05: Normalised S4 (Pre-RUV)
                    .capture_checkpoint(current_s4, "cp05", "normalised")

                    # 4. RUV-III
                    if (input$ruv_mode != "skip") {
                        add_log(paste("Running RUV-III:", input$ruv_mode))
                        ruv_params <- list(
                            percentage_min = input$auto_percentage_min,
                            percentage_max = input$auto_percentage_max,
                            ruv_grouping_variable = input$ruv_grouping_variable,
                            separation_metric = input$separation_metric,
                            k_penalty_weight = input$k_penalty_weight,
                            adaptive_k_penalty = input$adaptive_k_penalty,
                            manual_k = input$ruv_k,
                            manual_percentage = input$ruv_percentage
                        )
                        
                        ruv_results <- runPerAssayRuvOptimization(current_s4, input$ruv_mode, ruv_params, experiment_paths)
                        norm_data$ruv_results <- ruv_results
                        
                        best_k <- extractBestKPerAssay(ruv_results)
                        ctrl <- extractCtrlPerAssay(ruv_results)
                        
                        current_s4 <- ruvIII_C_Varying(current_s4, ruv_grouping_variable = input$ruv_grouping_variable, 
                                                     ruv_number_k = best_k, ctrl = ctrl)
                        
                        # Capture Checkpoint cp06: RUV Corrected S4
                        .capture_checkpoint(current_s4, "cp06", "ruv_corrected")
                        
                        add_log("RUV complete.")
                    }
                    
                    # Save State
                    workflow_data$state_manager$saveState("metab_normalized", current_s4, description = "Normalized metabolites")
                    
                    # Trigger Plots
                    generateMetabQcPlots(current_s4, experiment_paths, stage = "ruv_corrected", 
                                        grouping_variable = input$color_variable)
                    
                    norm_data$plot_refresh_trigger <- norm_data$plot_refresh_trigger + 1
                    add_log("Pipeline complete.")
                    shiny::showNotification("Normalization complete")
                    
                }, error = function(e) {
                    add_log(paste("Error:", e$message))
                    shiny::showNotification(e$message, type = "error")
                })
            })
        })

        output$norm_log <- shiny::renderText(paste(norm_data$normalization_log, collapse = "\n"))

        # Render RUV QC Plots
        output$ruv_qc_ui <- shiny::renderUI({
            shiny::req(norm_data$assay_names)
            assay_uis <- purrr::map(norm_data$assay_names, \(assay_name) {
                safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
                shiny::tagList(
                    shiny::h5(paste("Assay:", assay_name))
                    , shiny::plotOutput(ns(paste0("cancor_plot_", safe_name)), height = "300px")
                    , shiny::verbatimTextOutput(ns(paste0("ruv_summary_", safe_name)))
                    , shiny::hr()
                )
            })
            shiny::tagList(assay_uis)
        })

        shiny::observe({
            shiny::req(norm_data$ruv_results)
            purrr::iwalk(norm_data$ruv_results, \(res, assay_name) {
                safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
                output[[paste0("cancor_plot_", safe_name)]] <- shiny::renderPlot({ res$cancor_plot })
                output[[paste0("ruv_summary_", safe_name)]] <- shiny::renderText({
                    paste("Best K:", res$best_k, "\nControls:", sum(res$control_genes_index, na.rm=TRUE))
                })
            })
        })
    })
}

# <!-- APAF Bioinformatics | mod_metab_norm.R | Approved | 2026-03-17 -->
