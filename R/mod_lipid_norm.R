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
# mod_lipid_norm.R
# ============================================================================
# Purpose: Lipidomics normalization orchestrator module
# ============================================================================

#' @title Lipidomics Normalization Module
#' @description A Shiny module for multi-step normalization of lipidomics data.
#' @name mod_lipid_norm
NULL

#' @rdname mod_lipid_norm
#' @export
mod_lipid_norm_ui <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::wellPanel(
            shiny::fluidRow(
                # LEFT PANEL
                shiny::column(3
                    , shiny::wellPanel(
                        shiny::h4("Normalization Options")
                        , shiny::hr()
                        , shiny::selectInput(ns("color_variable"), "Color by:", choices = c("group"), width = "100%")
                        
                        , shiny::hr()
                        , shiny::h4("ITSD Normalization")
                        , shiny::checkboxInput(ns("apply_itsd"), "Apply ITSD normalization", value = TRUE)
                        , shiny::conditionalPanel(
                            condition = "input.apply_itsd == true", ns = ns
                            , shiny::selectInput(ns("itsd_aggregation"), "ITSD Aggregation:", 
                                               choices = c("Median" = "median", "Mean" = "mean", "Sum" = "sum"))
                        )

                        , shiny::hr()
                        , shiny::h4("Log2 Transform")
                        , shiny::numericInput(ns("log_offset"), "Offset:", value = 1, min = 0)

                        , shiny::hr()
                        , shiny::h4("Normalization")
                        , shiny::selectInput(ns("norm_method"), "Method:", 
                                           choices = list("Cyclic Loess" = "cyclicloess", "Quantile" = "quantile", "None" = "none"))

                        , shiny::hr()
                        , shiny::actionButton(ns("run_normalization"), "Run Pipeline", class = "btn-primary", width = "100%")
                    )
                )
                # RIGHT PANEL
                , shiny::column(9
                    , shiny::tabsetPanel(
                        id = ns("norm_qc_tabs")
                        , mod_lipid_norm_impute_ui(ns("impute"))
                        , mod_lipid_norm_itsd_ui(ns("itsd"))
                        , mod_lipid_norm_scale_ui(ns("pca"), "pca")
                        , mod_lipid_norm_scale_ui(ns("density"), "density")
                        , mod_lipid_norm_scale_ui(ns("rle"), "rle")
                        , shiny::tabPanel("Log", icon = shiny::icon("list-alt"), shiny::verbatimTextOutput(ns("norm_log")))
                    )
                )
            )
        )
    )
}

#' @rdname mod_lipid_norm
#' @export
mod_lipid_norm_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        norm_data <- shiny::reactiveValues(
            assay_names = NULL,
            plot_refresh_trigger = 0,
            normalization_log = character(0)
        )

        # Sub-modules
        mod_lipid_norm_impute_server("impute", workflow_data, omic_type, experiment_label)
        itsd_module <- mod_lipid_norm_itsd_server("itsd", workflow_data, omic_type, experiment_label)
        mod_lipid_norm_scale_server("pca", workflow_data, experiment_paths, norm_data, "pca")
        mod_lipid_norm_scale_server("density", workflow_data, experiment_paths, norm_data, "density")
        mod_lipid_norm_scale_server("rle", workflow_data, experiment_paths, norm_data, "rle")

        add_log <- function(msg) {
            norm_data$normalization_log <- c(norm_data$normalization_log, paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg))
        }

        shiny::observe({
            shiny::req(workflow_data$state_manager)
            current_s4 <- workflow_data$state_manager$getState()
            if (!is.null(current_s4) && inherits(current_s4, "LipidomicsAssayData")) {
                norm_data$assay_names <- names(current_s4@lipid_data)
            }
        })

        shiny::observe({
            if (!is.null(workflow_data$design_matrix)) {
                shiny::updateSelectInput(session, "color_variable", 
                                       choices = colnames(workflow_data$design_matrix),
                                       selected = "group")
            }
        })

        shiny::observeEvent(input$run_normalization, {
            shiny::req(workflow_data$state_manager)
            add_log("Starting lipidomics normalization...")
            
            shiny::withProgress(message = "Processing lipids...", value = 0, {
                tryCatch({
                    current_s4 <- workflow_data$state_manager$getState()
                    
                    if (isTRUE(input$apply_itsd)) {
                        add_log("ITSD normalization...")
                        ids <- itsd_module$get_feature_ids(current_s4)
                        current_s4 <- normaliseUntransformedData(current_s4, method = "ITSD", 
                                                               itsd_aggregation = input$itsd_aggregation, 
                                                               itsd_feature_ids = ids)
                    }
                    
                    add_log("Log2 transformation...")
                    current_s4 <- logTransformAssays(current_s4, offset = input$log_offset)
                    
                    if (input$norm_method != "none") {
                        add_log(paste("Method:", input$norm_method))
                        current_s4 <- normaliseBetweenSamples(current_s4, normalisation_method = input$norm_method)
                    }
                    
                    # Capture Checkpoint cp05: Normalised S4
                    .capture_checkpoint(current_s4, "cp05", "normalised")

                    workflow_data$state_manager$saveState("lipid_normalized", current_s4, description = "Normalized lipids")
                    
                    # NOTE: RUV-III logic for lipids could be added here if needed, 
                    # mirroring mod_metab_norm.R. For now, capturing cp06 as current state.
                    .capture_checkpoint(current_s4, "cp06", "ruv_corrected")
                    
                    generateLipidQcPlots(current_s4, experiment_paths, stage = "post_norm", 
                                       grouping_variable = input$color_variable)
                    
                    norm_data$plot_refresh_trigger <- norm_data$plot_refresh_trigger + 1
                    add_log("Lipid pipeline complete.")
                    shiny::showNotification("Lipidomics normalization complete")
                    
                }, error = function(e) {
                    add_log(paste("Error:", e$message))
                    shiny::showNotification(e$message, type = "error")
                })
            })
        })

        output$norm_log <- shiny::renderText(paste(norm_data$normalization_log, collapse = "\n"))
    })
}

# <!-- APAF Bioinformatics | mod_lipid_norm.R | Approved | 2026-03-17 -->
