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
# mod_metab_da.R
# ============================================================================
# Purpose: Metabolomics differential abundance analysis Shiny module
#
# This module provides per-assay differential abundance analysis using limma,
# with interactive volcano plots (Glimma), heatmaps, and result tables.
# UI and functionality matches the proteomics DA module (mod_prot_da.R).
# ============================================================================

#' @title Metabolomics Differential Analysis Module
#' @description A Shiny module for performing per-assay differential abundance
#'              analysis on metabolomics data using limma.
#' @name mod_metab_da
NULL

#' @rdname mod_metab_da
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 h4 h5 p hr br selectInput numericInput actionButton uiOutput verbatimTextOutput plotOutput tags downloadButton checkboxInput textAreaInput tabsetPanel tabPanel conditionalPanel helpText
#' @importFrom DT DTOutput
#' @importFrom shinyjqui jqui_resizable
mod_metab_da_ui <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::fluidRow(
            # ================================================================
            # LEFT SIDEBAR (3/12)
            # ================================================================
            shiny::column(
                3,
                shiny::wellPanel(
                    shiny::h4("Differential Abundance Analysis")

                    # Formula input
                    , shiny::textAreaInput(
                        ns("formula_string"),
                        "Model Formula:",
                        value = "~ 0 + group",
                        height = "60px"
                    ),
                    shiny::helpText("Formula fetched from S4 object @args if available"),
                    shiny::hr()

                    # Analysis Parameters
                    , shiny::h5("Analysis Parameters"),
                    shiny::numericInput(
                        ns("da_q_val_thresh"),
                        "Q-value Threshold:",
                        value = 0.05,
                        min = 0.001,
                        max = 0.2,
                        step = 0.005
                    ),
                    shiny::numericInput(
                        ns("treat_lfc_cutoff"),
                        "Log Fold-Change Cutoff:",
                        value = 0,
                        min = 0,
                        max = 2,
                        step = 0.1
                    ),
                    shiny::hr()

                    # Available contrasts display
                    , shiny::h5("Available Contrasts"),
                    shiny::verbatimTextOutput(ns("contrasts_display")),
                    shiny::helpText("Contrasts defined in design matrix"),
                    shiny::hr()

                    # Load Filtered Session
                    , shiny::actionButton(
                        ns("load_filtered_session"),
                        "Load Filtered Session",
                        class = "btn-info",
                        width = "100%",
                        icon = shiny::icon("upload")
                    ),
                    shiny::helpText("Load previously exported filtered data for DA analysis"),
                    shiny::br(),
                    shiny::br()

                    # Run DA Analysis
                    , shiny::actionButton(
                        ns("run_da_analysis"),
                        "Run DA Analysis",
                        class = "btn-primary",
                        width = "100%",
                        icon = shiny::icon("play")
                    ),
                    shiny::br(),
                    shiny::br()

                    # Download Results
                    , shiny::downloadButton(
                        ns("download_da_results"),
                        "Download All Results",
                        class = "btn-success",
                        style = "width: 100%;"
                    ),
                    shiny::hr()

                    # Status display
                    , shiny::h5("Analysis Status"),
                    shiny::verbatimTextOutput(ns("da_status"))
                )
            )

            # ================================================================
            # RIGHT PANEL - TABSET (9/12)
            # ================================================================
            , shiny::column(
                9,
                shiny::uiOutput(ns("heatmap_manual_save_warning")),
                shiny::tabsetPanel(
                    id = ns("da_results_tabs")

                    # ----------------------------------------------------------
                    # TAB 1: VOLCANO PLOT
                    # ----------------------------------------------------------
                    , shiny::tabPanel(
                        "Volcano Plot",
                        shiny::br(),
                        shiny::fluidRow(
                            shiny::column(
                                4,
                                shiny::selectInput(
                                    ns("volcano_contrast"),
                                    "Select Contrast:",
                                    choices = NULL,
                                    width = "100%"
                                )
                            ),
                            shiny::column(
                                4,
                                shiny::selectInput(
                                    ns("volcano_assay"),
                                    "Select Assay:",
                                    choices = c("Combined"),
                                    width = "100%"
                                )
                            ),
                            shiny::column(
                                4,
                                shiny::checkboxInput(
                                    ns("volcano_interactive"),
                                    "Interactive Plot (Glimma)",
                                    value = TRUE
                                )
                            )
                        )

                        # Interactive volcano (Glimma)
                        , shiny::conditionalPanel(
                            condition = "input.volcano_interactive == true",
                            ns = ns,
                            shiny::div(
                                id = ns("volcano_glimma_container"),
                                style = "height: 600px; background-color: white; border-radius: 8px; padding: 10px;",
                                shiny::htmlOutput(ns("volcano_glimma"))
                            )
                        )

                        # Static volcano (ggplot2)
                        , shiny::conditionalPanel(
                            condition = "input.volcano_interactive == false",
                            ns = ns,
                            shinyjqui::jqui_resizable(
                                shiny::plotOutput(ns("volcano_static"), height = "600px")
                            )
                        )
                    )

                    # ----------------------------------------------------------
                    # TAB 2: HEATMAP
                    # ----------------------------------------------------------
                    , shiny::tabPanel(
                        "Heatmap",
                        shiny::br()

                        # Main controls (4 columns)
                        , shiny::fluidRow(
                            shiny::column(
                                3,
                                shiny::selectInput(
                                    ns("heatmap_contrast"),
                                    "Select Contrast:",
                                    choices = NULL,
                                    width = "100%"
                                )
                            ),
                            shiny::column(
                                3,
                                shiny::numericInput(
                                    ns("heatmap_top_n"),
                                    "Top N Metabolites:",
                                    value = 50,
                                    min = 10,
                                    max = 500,
                                    step = 10
                                )
                            ),
                            shiny::column(
                                3,
                                shiny::selectInput(
                                    ns("heatmap_clustering"),
                                    "Apply Clustering:",
                                    choices = list(
                                        "Both rows & columns" = "both",
                                        "Rows only" = "row",
                                        "Columns only" = "column",
                                        "None" = "none"
                                    ),
                                    selected = "both"
                                )
                            ),
                            shiny::column(
                                3,
                                shiny::selectInput(
                                    ns("heatmap_scaling"),
                                    "Data Scaling:",
                                    choices = list(
                                        "Row (metabolite) scaling" = "row",
                                        "Column (sample) scaling" = "column",
                                        "Both" = "both",
                                        "None" = "none"
                                    ),
                                    selected = "row"
                                )
                            )
                        )

                        # Advanced clustering options (conditional)
                        , shiny::conditionalPanel(
                            condition = "input.heatmap_clustering != 'none'",
                            ns = ns,
                            shiny::wellPanel(
                                shiny::h5("Clustering Options"),
                                shiny::fluidRow(
                                    shiny::column(
                                        4,
                                        shiny::selectInput(
                                            ns("heatmap_cluster_method"),
                                            "Clustering Method:",
                                            choices = list(
                                                "Ward (minimum variance)" = "ward.D2",
                                                "Ward (original)" = "ward.D",
                                                "Complete linkage" = "complete",
                                                "Single linkage" = "single",
                                                "Average linkage" = "average",
                                                "McQuitty (WPGMA)" = "mcquitty"
                                            ),
                                            selected = "ward.D2"
                                        )
                                    ),
                                    shiny::column(
                                        4,
                                        shiny::selectInput(
                                            ns("heatmap_distance_method"),
                                            "Distance Metric:",
                                            choices = list(
                                                "Euclidean" = "euclidean",
                                                "Manhattan" = "manhattan",
                                                "Pearson correlation" = "pearson",
                                                "Spearman correlation" = "spearman",
                                                "Maximum" = "maximum"
                                            ),
                                            selected = "euclidean"
                                        )
                                    ),
                                    shiny::column(
                                        4,
                                        shiny::selectInput(
                                            ns("heatmap_color_scheme"),
                                            "Color Scheme:",
                                            choices = list(
                                                "Red-Blue" = "RdBu",
                                                "Red-Yellow-Blue" = "RdYlBu",
                                                "Blue-White-Red" = "coolwarm",
                                                "Viridis" = "viridis",
                                                "Plasma" = "plasma",
                                                "Inferno" = "inferno"
                                            ),
                                            selected = "RdBu"
                                        )
                                    )
                                ),
                                shiny::fluidRow(
                                    shiny::column(
                                        4,
                                        shiny::selectInput(
                                            ns("heatmap_assay"),
                                            "Select Assay:",
                                            choices = c("Combined"),
                                            width = "100%"
                                        )
                                    ),
                                    shiny::column(
                                        4,
                                        shiny::checkboxInput(
                                            ns("heatmap_show_labels"),
                                            "Show Metabolite Labels",
                                            value = FALSE
                                        )
                                    )
                                )
                            )
                        )

                        # Heatmap output
                        , shiny::div(
                            id = ns("heatmap_container"),
                            style = "height: 650px;",
                            shinyjqui::jqui_resizable(
                                shiny::plotOutput(ns("heatmap_plot"), height = "600px")
                            )
                        )

                        # Save Button
                        , shiny::div(
                            style = "margin-top: 10px; margin-bottom: 20px;",
                            shiny::actionButton(
                                ns("save_heatmap"),
                                "Save Heatmap",
                                icon = shiny::icon("download"),
                                class = "btn-primary"
                            )
                        )

                        # Cluster summary (if tree cutting is applied)
                        , shiny::conditionalPanel(
                            condition = "input.heatmap_tree_cut_method != 'none'",
                            ns = ns,
                            shiny::wellPanel(
                                shiny::h4("Cluster Information"),
                                shiny::verbatimTextOutput(ns("cluster_summary"))
                            )
                        )
                    )

                    # ----------------------------------------------------------
                    # TAB 3: DA RESULTS TABLE
                    # ----------------------------------------------------------
                    , shiny::tabPanel(
                        "DA Results Table",
                        shiny::br()

                        # Controls
                        , shiny::fluidRow(
                            shiny::column(
                                3,
                                shiny::selectInput(
                                    ns("table_contrast"),
                                    "Select Contrast:",
                                    choices = NULL,
                                    width = "100%"
                                )
                            ),
                            shiny::column(
                                3,
                                shiny::selectInput(
                                    ns("table_assay"),
                                    "Select Assay:",
                                    choices = c("All"),
                                    width = "100%"
                                )
                            ),
                            shiny::column(
                                3,
                                shiny::selectInput(
                                    ns("table_significance"),
                                    "Show:",
                                    choices = list(
                                        "All Results" = "all",
                                        "Significant Only" = "significant",
                                        "Up-regulated" = "up",
                                        "Down-regulated" = "down"
                                    ),
                                    selected = "significant"
                                )
                            ),
                            shiny::column(
                                3,
                                shiny::numericInput(
                                    ns("table_max_rows"),
                                    "Max Rows:",
                                    value = 1000,
                                    min = 100,
                                    max = 10000,
                                    step = 100
                                )
                            )
                        )

                        # Summary statistics
                        , shiny::fluidRow(
                            shiny::column(
                                12,
                                shiny::wellPanel(
                                    shiny::h5("Summary Statistics"),
                                    shiny::verbatimTextOutput(ns("da_summary_stats"))
                                )
                            )
                        )

                        # Results table
                        , DT::DTOutput(ns("da_results_table"))
                    )
                )
            )
        )
    )
}





















































