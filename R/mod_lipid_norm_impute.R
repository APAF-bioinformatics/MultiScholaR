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
# mod_lipid_norm_impute.R
# ============================================================================
# Purpose: Lipidomics imputation Shiny module
# ============================================================================

#' @title Lipidomics Imputation Module
#' @description A Shiny module for applying missing value imputation to lipidomics data.
#' @name mod_lipid_norm_impute
NULL

#' @rdname mod_lipid_norm_impute
#' @export
mod_lipid_norm_impute_ui <- function(id) {
    ns <- shiny::NS(id)

    shiny::tabPanel(
        "Imputation",
        shiny::br(),
        shiny::tabsetPanel(
            id = ns("impute_internal_tabs"),
            shiny::tabPanel(
                "Controls",
                shiny::br(),
                shiny::fluidRow(
                    shiny::column(
                        4,
                        shiny::wellPanel(
                            shiny::h4("Missing Value Imputation"),
                            shiny::p("Handle missing values (zeros or NAs) in your lipidomics data."),
                            shiny::hr(),
                            shiny::selectInput(
                                ns("impute_method"),
                                "Imputation Method:",
                                choices = list(
                                    "None (skip)" = "none",
                                    "Limpa (Probabilistic)" = "limpa",
                                    "MissForest (Random Forest)" = "missforest",
                                    "Half-Minimum" = "half_min",
                                    "Zero" = "zero"
                                ),
                                selected = "none"
                            ),
                            shiny::hr(),
                            shiny::div(
                                shiny::actionButton(
                                    ns("apply_imputation"),
                                    "Apply Imputation",
                                    class = "btn-primary",
                                    width = "48%"
                                ),
                                shiny::actionButton(
                                    ns("revert_imputation"),
                                    "Revert",
                                    class = "btn-warning",
                                    width = "48%",
                                    style = "margin-left: 4%"
                                )
                            )
                        )
                    ),
                    shiny::column(8, shiny::verbatimTextOutput(ns("imputation_results")))
                )
            ),
            shiny::tabPanel(
                "Diagnostics",
                shiny::br(),
                shiny::uiOutput(ns("diagnostics_ui"))
            )
        )
    )
}

#' @rdname mod_lipid_norm_impute
#' @export
mod_lipid_norm_impute_server <- function(id, workflow_data, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Reactive values for diagnostics
        diag_data <- shiny::reactiveValues(
            dpc_plots = list(),
            missforest_logs = list(),
            current_method = "none"
        )

        shiny::observeEvent(input$apply_imputation, {
            shiny::req(workflow_data$state_manager)
            if (input$impute_method == "none") {
                return()
            }

            shiny::showNotification("Applying imputation...", id = "lipid_impute_working", duration = NULL)

            tryCatch(
                {
                    current_s4 <- workflow_data$state_manager$getState()
                    diag_data$current_method <- input$impute_method
                    imputed_s4 <- current_s4

                    if (input$impute_method == "limpa") {
                        imputed_s4 <- lipidMissingValueImputationLimpa(current_s4, verbose = TRUE)
                        # Capture DPC info
                        diag_data$dpc_plots <- purrr::map(names(imputed_s4@lipid_data), \(assay_name) {
                            assay_df <- imputed_s4@lipid_data[[assay_name]]
                            id_col <- imputed_s4@lipid_id_column
                            design_samples <- as.character(imputed_s4@design_matrix[[imputed_s4@sample_id]])
                            sample_cols <- intersect(colnames(assay_df), design_samples)
                            mat <- assay_df |>
                                dplyr::select(all_of(c(id_col, sample_cols))) |>
                                tibble::column_to_rownames(id_col) |>
                                as.matrix()
                            if (max(mat, na.rm = TRUE) > 50) mat <- log2(mat)
                            if (any(is.infinite(mat))) mat[is.infinite(mat)] <- NA

                            dpc_params <- tryCatch(
                                {
                                    limpa::dpc(mat)
                                },
                                error = function(e) NULL
                            )

                            if (!is.null(dpc_params) && !is.null(experiment_label)) {
                                tryCatch(
                                    {
                                        qc_dir <- file.path("results", experiment_label, "QC", "Imputation")
                                        if (!dir.exists(qc_dir)) dir.create(qc_dir, recursive = TRUE)
                                        safe_assay <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
                                        x_range <- seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out = 100)
                                        prob <- 1 / (1 + exp(-(dpc_params$alpha + dpc_params$beta * x_range)))
                                        df_plot <- data.frame(Intensity = x_range, Probability = prob)
                                        p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Intensity, y = Probability)) +
                                            ggplot2::geom_line(color = "darkred", size = 1) +
                                            ggplot2::ylim(0, 1) +
                                            ggplot2::labs(title = paste("DPC (Lipid):", assay_name), subtitle = sprintf("Alpha: %.2f, Beta: %.2f", dpc_params$alpha, dpc_params$beta)) +
                                            ggplot2::theme_minimal()
                                        ggplot2::ggsave(file.path(qc_dir, paste0("dpc_lipid_", safe_assay, ".png")), p, width = 6, height = 4)
                                        saveRDS(dpc_params, file.path(qc_dir, paste0("dpc_params_lipid_", safe_assay, ".rds")))
                                    },
                                    error = function(e) NULL
                                )
                            }
                            return(list(params = dpc_params))
                        }) |> purrr::set_names(names(imputed_s4@lipid_data))
                    } else if (input$impute_method == "missforest") {
                        imputed_s4 <- lipidMissingValueImputationMissForest(current_s4, verbose = TRUE)
                        diag_data$missforest_logs <- list(status = "Complete", timestamp = Sys.time())

                        if (!is.null(experiment_label)) {
                            tryCatch(
                                {
                                    qc_dir <- file.path("results", experiment_label, "QC", "Imputation")
                                    if (!dir.exists(qc_dir)) dir.create(qc_dir, recursive = TRUE)
                                    saveRDS(diag_data$missforest_logs, file.path(qc_dir, "missforest_diagnostics_lipid.rds"))
                                },
                                error = function(e) {
                                    logger::log_warn("Failed to save missForest diagnostics")
                                }
                            )
                        }
                    } else {
                        imputed_s4@lipid_data <- purrr::map(current_s4@lipid_data, function(assay_df) {
                            num_cols <- sapply(assay_df, is.numeric)
                            sample_cols <- names(assay_df)[num_cols]

                            assay_df[sample_cols] <- purrr::map_dfc(assay_df[sample_cols], function(col) {
                                if (input$impute_method == "zero") {
                                    col[is.na(col)] <- 0
                                } else if (input$impute_method == "half_min") {
                                    # FIX: Prevent warning
                                    pos_vals <- col[col > 0 & !is.na(col)]
                                    if (length(pos_vals) == 0) {
                                        min_val <- 1
                                    } else {
                                        min_val <- min(pos_vals, na.rm = TRUE)
                                    }
                                    col[is.na(col) | col == 0] <- min_val / 2
                                }
                                return(col)
                            })
                            return(assay_df)
                        })
                    }

                    # Capture Checkpoint cp02_imputed
                    if (exists(".capture_checkpoint")) {
                        .capture_checkpoint(imputed_s4, "cp02", "imputed")
                    }

                    workflow_data$state_manager$saveState("lipid_imputed", imputed_s4, description = paste("Imputed lipids:", input$impute_method))

                    output$imputation_results <- shiny::renderText("Imputation applied successfully.")
                    shiny::removeNotification("lipid_impute_working")
                    shiny::showNotification("Imputation complete")
                },
                error = function(e) {
                    shiny::removeNotification("lipid_impute_working")
                    shiny::showNotification(e$message, type = "error")
                }
            )
        })

        # Revert
        shiny::observeEvent(input$revert_imputation, {
            tryCatch(
                {
                    history <- workflow_data$state_manager$getHistory()
                    if (length(history) > 1) {
                        prev_state <- history[length(history) - 1]
                        workflow_data$state_manager$revertToState(prev_state)
                        output$imputation_results <- shiny::renderText(paste("Reverted to:", prev_state))
                    }
                },
                error = function(e) {
                    shiny::showNotification(e$message, type = "error")
                }
            )
        })

        # Diagnostics UI
        output$diagnostics_ui <- shiny::renderUI({
            if (diag_data$current_method == "limpa") {
                shiny::tagList(
                    shiny::h4("Limpa DPC Diagnostics"),
                    shiny::verbatimTextOutput(ns("limpa_diag_text"))
                )
            } else if (diag_data$current_method == "missforest") {
                shiny::tagList(
                    shiny::h4("MissForest OOB Report"),
                    shiny::verbatimTextOutput(ns("missforest_diag_text"))
                )
            } else {
                shiny::p("No diagnostics available.")
            }
        })

        output$limpa_diag_text <- shiny::renderPrint({
            diag_data$dpc_plots
        })
        output$missforest_diag_text <- shiny::renderPrint({
            diag_data$missforest_logs
        })
    })
}

# <!-- APAF Bioinformatics | mod_lipid_norm_impute.R | Approved | 2026-03-17 -->
