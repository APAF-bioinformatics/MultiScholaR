# ============================================================================
# mod_lipid_norm.R
# ============================================================================
# Purpose: Lipidomics normalization Shiny module (Proteomics UX harmonized)
#
# This module provides a unified normalization pipeline with per-assay
# visualization for lipidomics data (pos/neg mode support).
#
# Key features:
# - 3/9 column layout matching proteomics pattern
# - All options visible (single "Run" button)
# - Per-assay ITSD selection with DT::datatable
# - Independent auto-RUV optimization per assay
# - 3-column QC comparison (Post-Filter / Post-Norm / RUV-Corrected)
# - Per-assay visualization (pos row on top, neg row below)
# - Correlation filtering as final step
# ============================================================================

#' @title Lipidomics Normalization Module
#' @description A Shiny module for multi-step normalization of lipidomics data.
#'              Harmonized with proteomics UX pattern.
#' @name mod_lipid_norm
NULL

#' @rdname mod_lipid_norm
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 h4 h5 p hr br radioButtons selectInput textInput actionButton uiOutput verbatimTextOutput icon tags checkboxInput numericInput plotOutput conditionalPanel tabsetPanel tabPanel sliderInput helpText
#' @importFrom shinyjqui jqui_resizable
#' @importFrom DT dataTableOutput
mod_lipid_norm_ui <- function(id) {
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
                        , shiny::helpText("Variable to use for plot coloring")

                        , shiny::selectInput(
                            ns("shape_variable")
                            , "Shape by:"
                            , choices = c("group", "factor1", "factor2", "batch")
                            , selected = "group"
                            , width = "100%"
                        )
                        , shiny::helpText("Variable to use for plot shapes")

                        # --- ITSD Normalization ---
                        , shiny::hr()
                        , shiny::h4("ITSD Normalization")
                        , shiny::checkboxInput(
                            ns("apply_itsd")
                            , "Apply ITSD normalization"
                            , value = TRUE
                            , width = "100%"
                        )
                        , shiny::helpText("Normalize using internal standards (per assay)")

                        , shiny::conditionalPanel(
                            condition = "input.apply_itsd == true"
                            , ns = ns
                            , shiny::selectInput(
                                ns("itsd_aggregation")
                                , "ITSD Aggregation:"
                                , choices = c(
                                    "Median (recommended)" = "median"
                                    , "Mean" = "mean"
                                    , "Sum" = "sum"
                                )
                                , selected = "median"
                                , width = "100%"
                            )
                            , shiny::helpText("Method to combine multiple internal standards per sample. Select standards in ITSD Selection tab.")
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
                        , shiny::helpText("Applied after ITSD normalization")

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
                        , shiny::helpText("Method for between-sample normalization")

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
                        , shiny::helpText("Experimental factor to preserve during batch correction (shared across assays)")

                        # --- Automatic RUV Parameters ---
                        , shiny::conditionalPanel(
                            condition = "input.ruv_mode == 'automatic'"
                            , ns = ns

                            , shiny::h5("Automatic Optimization Parameters")
                            , shiny::helpText("Auto-optimization runs independently for each assay (pos/neg)")

                            , shiny::fluidRow(
                                shiny::column(6
                                    , shiny::numericInput(
                                        ns("auto_percentage_min")
                                        , "Min %:"
                                        , value = 1
                                        , min = 1
                                        , max = 50
                                        , width = "100%"
                                    )
                                )
                                , shiny::column(6
                                    , shiny::numericInput(
                                        ns("auto_percentage_max")
                                        , "Max %:"
                                        , value = 20
                                        , min = 1
                                        , max = 50
                                        , width = "100%"
                                    )
                                )
                            )
                            , shiny::helpText("Range of percentages to test (1-50%)")

                            , shiny::selectInput(
                                ns("separation_metric")
                                , "Separation Quality Metric:"
                                , choices = list(
                                    "Maximum Difference (recommended)" = "max_difference"
                                    , "Mean Difference" = "mean_difference"
                                    , "Area Under Curve" = "auc"
                                    , "Weighted Difference" = "weighted_difference"
                                )
                                , selected = "max_difference"
                                , width = "100%"
                            )

                            , shiny::sliderInput(
                                ns("k_penalty_weight")
                                , "K Value Penalty Weight:"
                                , min = 0.1
                                , max = 0.9
                                , value = 0.5
                                , step = 0.1
                                , width = "100%"
                            )
                            , shiny::helpText("Controls penalty for high k values (0.1 = lenient, 0.9 = strict)")

                            , shiny::checkboxInput(
                                ns("adaptive_k_penalty")
                                , "Adaptive K Penalty (sample size aware)"
                                , value = TRUE
                                , width = "100%"
                            )
                            , shiny::br()
                        )

                        # --- Manual RUV Parameters ---
                        , shiny::conditionalPanel(
                            condition = "input.ruv_mode == 'manual'"
                            , ns = ns

                            , shiny::sliderInput(
                                ns("ruv_percentage")
                                , "% Lipids as Negative Controls:"
                                , min = 1
                                , max = 20
                                , value = 5
                                , step = 1
                                , width = "100%"
                            )
                            , shiny::helpText("Percentage of most stable lipids used as negative controls")

                            , shiny::numericInput(
                                ns("ruv_k")
                                , "Number of Factors (k):"
                                , value = 2
                                , min = 1
                                , max = 10
                                , width = "100%"
                            )
                            , shiny::helpText("Number of unwanted variation factors to remove")
                            , shiny::br()
                        )

                        # --- Skip RUV Info ---
                        , shiny::conditionalPanel(
                            condition = "input.ruv_mode == 'skip'"
                            , ns = ns

                            , shiny::tags$div(
                                style = "background-color: #fff3cd; border: 1px solid #ffc107; padding: 10px; border-radius: 4px;"
                                , shiny::h5("RUV-III Will Be Skipped", style = "color: #856404; margin-top: 0;")
                                , shiny::helpText(
                                    "Batch correction will not be applied. Consider including batch as a covariate in DE model if batch effects are present."
                                    , style = "margin-bottom: 0;"
                                )
                            )
                            , shiny::br()
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

                        , shiny::br()
                        , shiny::br()

                        , shiny::actionButton(
                            ns("reset_normalization")
                            , "Reset to Pre-Normalization"
                            , class = "btn-warning"
                            , width = "100%"
                            , icon = shiny::icon("undo")
                        )

                        , shiny::br()
                        , shiny::br()

                        , shiny::actionButton(
                            ns("export_session")
                            , "Export Normalized Data"
                            , class = "btn-success"
                            , width = "100%"
                            , icon = shiny::icon("download")
                        )
                        , shiny::helpText("Save normalized data for later analysis")
                    )
                )

                # ================================================================
                # RIGHT PANEL (width = 9): QC Visualization Tabs
                # ================================================================
                , shiny::column(9
                    , shiny::tabsetPanel(
                        id = ns("norm_qc_tabs")

                        # --------------------------------------------------------
                        # Tab 1: ITSD Selection (per-assay DT tables)
                        # --------------------------------------------------------
                        , shiny::tabPanel(
                            "ITSD Selection"
                            , icon = shiny::icon("flask")
                            , shiny::br()
                            , shiny::fluidRow(
                                shiny::column(12
                                    , shiny::h5("Internal Standard Selection per Assay", style = "text-align: center;")
                                    , shiny::helpText("Select internal standards for each assay. Pre-selected candidates are auto-detected from naming patterns (ITSD, IS, labeled compounds).")
                                    , shiny::hr()
                                    , shiny::uiOutput(ns("itsd_selection_ui"))
                                )
                            )
                        )

                        # --------------------------------------------------------
                        # Tab 2: RUV QC (per-assay stacked cancor plots)
                        # --------------------------------------------------------
                        , shiny::tabPanel(
                            "RUV QC"
                            , icon = shiny::icon("chart-line")
                            , shiny::br()
                            , shiny::fluidRow(
                                shiny::column(12
                                    , shiny::h5("RUV Parameter Optimization Results", style = "text-align: center;")
                                    , shiny::helpText("Canonical correlation plots and optimization summary per assay. Each assay is optimized independently.")
                                    , shiny::hr()
                                    , shiny::uiOutput(ns("ruv_qc_ui"))
                                )
                            )
                        )

                        # --------------------------------------------------------
                        # Tab 3: PCA (STATIC 2-row × 3-column grid)
                        # --------------------------------------------------------
                        , shiny::tabPanel(
                            "PCA"
                            , icon = shiny::icon("project-diagram")
                            , shiny::br()
                            # --- Row 1: Assay 1 ---
                            , shiny::fluidRow(
                                shiny::column(12
                                    , shiny::h5(shiny::textOutput(ns("assay1_label_pca")), style = "margin-bottom: 5px;")
                                )
                            )
                            , shiny::fluidRow(
                                shiny::column(4
                                    , shiny::tags$small("Post-Filtering", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("pca_post_filter_assay1"), height = "300px")
                                    )
                                )
                                , shiny::column(4
                                    , shiny::tags$small("Post-Normalization", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("pca_post_norm_assay1"), height = "300px")
                                    )
                                )
                                , shiny::column(4
                                    , shiny::tags$small("RUV-Corrected", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("pca_ruv_corrected_assay1"), height = "300px")
                                    )
                                )
                            )
                            , shiny::hr()
                            # --- Row 2: Assay 2 ---
                            , shiny::fluidRow(
                                shiny::column(12
                                    , shiny::h5(shiny::textOutput(ns("assay2_label_pca")), style = "margin-bottom: 5px;")
                                )
                            )
                            , shiny::fluidRow(
                                shiny::column(4
                                    , shiny::tags$small("Post-Filtering", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("pca_post_filter_assay2"), height = "300px")
                                    )
                                )
                                , shiny::column(4
                                    , shiny::tags$small("Post-Normalization", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("pca_post_norm_assay2"), height = "300px")
                                    )
                                )
                                , shiny::column(4
                                    , shiny::tags$small("RUV-Corrected", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("pca_ruv_corrected_assay2"), height = "300px")
                                    )
                                )
                            )
                        )

                        # --------------------------------------------------------
                        # Tab 4: Density (STATIC 2-row × 3-column grid)
                        # --------------------------------------------------------
                        , shiny::tabPanel(
                            "Density"
                            , icon = shiny::icon("chart-area")
                            , shiny::br()
                            # --- Row 1: Assay 1 ---
                            , shiny::fluidRow(
                                shiny::column(12
                                    , shiny::h5(shiny::textOutput(ns("assay1_label_density")), style = "margin-bottom: 5px;")
                                )
                            )
                            , shiny::fluidRow(
                                shiny::column(4
                                    , shiny::tags$small("Post-Filtering", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("density_post_filter_assay1"), height = "300px")
                                    )
                                )
                                , shiny::column(4
                                    , shiny::tags$small("Post-Normalization", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("density_post_norm_assay1"), height = "300px")
                                    )
                                )
                                , shiny::column(4
                                    , shiny::tags$small("RUV-Corrected", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("density_ruv_corrected_assay1"), height = "300px")
                                    )
                                )
                            )
                            , shiny::hr()
                            # --- Row 2: Assay 2 ---
                            , shiny::fluidRow(
                                shiny::column(12
                                    , shiny::h5(shiny::textOutput(ns("assay2_label_density")), style = "margin-bottom: 5px;")
                                )
                            )
                            , shiny::fluidRow(
                                shiny::column(4
                                    , shiny::tags$small("Post-Filtering", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("density_post_filter_assay2"), height = "300px")
                                    )
                                )
                                , shiny::column(4
                                    , shiny::tags$small("Post-Normalization", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("density_post_norm_assay2"), height = "300px")
                                    )
                                )
                                , shiny::column(4
                                    , shiny::tags$small("RUV-Corrected", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("density_ruv_corrected_assay2"), height = "300px")
                                    )
                                )
                            )
                        )

                        # --------------------------------------------------------
                        # Tab 5: RLE (STATIC 2-row × 3-column grid)
                        # --------------------------------------------------------
                        , shiny::tabPanel(
                            "RLE"
                            , icon = shiny::icon("chart-bar")
                            , shiny::br()
                            # --- Row 1: Assay 1 ---
                            , shiny::fluidRow(
                                shiny::column(12
                                    , shiny::h5(shiny::textOutput(ns("assay1_label_rle")), style = "margin-bottom: 5px;")
                                )
                            )
                            , shiny::fluidRow(
                                shiny::column(4
                                    , shiny::tags$small("Post-Filtering", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("rle_post_filter_assay1"), height = "300px")
                                    )
                                )
                                , shiny::column(4
                                    , shiny::tags$small("Post-Normalization", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("rle_post_norm_assay1"), height = "300px")
                                    )
                                )
                                , shiny::column(4
                                    , shiny::tags$small("RUV-Corrected", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("rle_ruv_corrected_assay1"), height = "300px")
                                    )
                                )
                            )
                            , shiny::hr()
                            # --- Row 2: Assay 2 ---
                            , shiny::fluidRow(
                                shiny::column(12
                                    , shiny::h5(shiny::textOutput(ns("assay2_label_rle")), style = "margin-bottom: 5px;")
                                )
                            )
                            , shiny::fluidRow(
                                shiny::column(4
                                    , shiny::tags$small("Post-Filtering", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("rle_post_filter_assay2"), height = "300px")
                                    )
                                )
                                , shiny::column(4
                                    , shiny::tags$small("Post-Normalization", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("rle_post_norm_assay2"), height = "300px")
                                    )
                                )
                                , shiny::column(4
                                    , shiny::tags$small("RUV-Corrected", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("rle_ruv_corrected_assay2"), height = "300px")
                                    )
                                )
                            )
                        )

                        # --------------------------------------------------------
                        # Tab 6: Correlation Heatmap (STATIC 2-row × 3-column grid)
                        # --------------------------------------------------------
                        , shiny::tabPanel(
                            "Correlation"
                            , icon = shiny::icon("th")
                            , shiny::br()
                            # --- Row 1: Assay 1 ---
                            , shiny::fluidRow(
                                shiny::column(12
                                    , shiny::h5(shiny::textOutput(ns("assay1_label_correlation")), style = "margin-bottom: 5px;")
                                )
                            )
                            , shiny::fluidRow(
                                shiny::column(4
                                    , shiny::tags$small("Post-Filtering", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("correlation_post_filter_assay1"), height = "300px")
                                    )
                                )
                                , shiny::column(4
                                    , shiny::tags$small("Post-Normalization", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("correlation_post_norm_assay1"), height = "300px")
                                    )
                                )
                                , shiny::column(4
                                    , shiny::tags$small("RUV-Corrected", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("correlation_ruv_corrected_assay1"), height = "300px")
                                    )
                                )
                            )
                            , shiny::hr()
                            # --- Row 2: Assay 2 ---
                            , shiny::fluidRow(
                                shiny::column(12
                                    , shiny::h5(shiny::textOutput(ns("assay2_label_correlation")), style = "margin-bottom: 5px;")
                                )
                            )
                            , shiny::fluidRow(
                                shiny::column(4
                                    , shiny::tags$small("Post-Filtering", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("correlation_post_filter_assay2"), height = "300px")
                                    )
                                )
                                , shiny::column(4
                                    , shiny::tags$small("Post-Normalization", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("correlation_post_norm_assay2"), height = "300px")
                                    )
                                )
                                , shiny::column(4
                                    , shiny::tags$small("RUV-Corrected", style = "display: block; text-align: center;")
                                    , shinyjqui::jqui_resizable(
                                        shiny::imageOutput(ns("correlation_ruv_corrected_assay2"), height = "300px")
                                    )
                                )
                            )
                        )

                        # --------------------------------------------------------
                        # Tab 7: Correlation Filtering (final step)
                        # --------------------------------------------------------
                        , shiny::tabPanel(
                            "Correlation Filtering"
                            , icon = shiny::icon("filter")
                            , shiny::br()
                            , shiny::fluidRow(
                                shiny::column(4
                                    , shiny::wellPanel(
                                        shiny::h5("Sample Correlation Filtering")
                                        , shiny::helpText("FINAL STEP: Filter samples based on their correlation with other samples in the same group. Low correlation indicates potential technical issues.")

                                        , shiny::sliderInput(
                                            ns("min_pearson_correlation_threshold")
                                            , "Min. Pearson Correlation Threshold:"
                                            , min = 0.0
                                            , max = 1.0
                                            , value = 0.5
                                            , step = 0.05
                                            , width = "100%"
                                        )
                                        , shiny::helpText("Samples with correlation below this threshold will be flagged")

                                        , shiny::br()

                                        , shiny::actionButton(
                                            ns("apply_correlation_filter")
                                            , "Apply Correlation Filter"
                                            , class = "btn-primary"
                                            , width = "100%"
                                            , icon = shiny::icon("play")
                                        )

                                        , shiny::br()
                                        , shiny::br()

                                        , shiny::actionButton(
                                            ns("skip_correlation_filter")
                                            , "Skip Filtering & Proceed to DE"
                                            , class = "btn-info"
                                            , width = "100%"
                                            , icon = shiny::icon("forward")
                                        )
                                        , shiny::helpText("Bypass sample filtering and proceed directly to DE")

                                        , shiny::br()
                                        , shiny::br()

                                        , shiny::h5("Filter Results")
                                        , shiny::verbatimTextOutput(ns("correlation_filter_summary"))
                                    )
                                )
                                , shiny::column(8
                                    , shiny::h5("Final Filtered Data QC", style = "text-align: center;")
                                    , shiny::helpText("QC plots for the final dataset after correlation filtering")
                                    , shinyjqui::jqui_resizable(
                                        shiny::plotOutput(ns("final_qc_plot"), height = "500px")
                                    )
                                )
                            )
                        )

                        # --------------------------------------------------------
                        # Tab 8: Normalization Log
                        # --------------------------------------------------------
                        , shiny::tabPanel(
                            "Log"
                            , icon = shiny::icon("list-alt")
                            , shiny::br()
                            , shiny::fluidRow(
                                shiny::column(12
                                    , shiny::h5("Normalization Pipeline Log", style = "text-align: center;")
                                    , shiny::verbatimTextOutput(ns("norm_log"))
                                )
                            )
                        )
                    )
                )
            )
        )
    )
}


#' @rdname mod_lipid_norm
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent req renderUI renderText renderPlot showNotification removeNotification tags renderImage updateSelectInput observe
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_density geom_point labs theme_minimal theme facet_wrap scale_color_brewer coord_flip
#' @importFrom logger log_info log_error log_warn
#' @importFrom purrr imap map walk set_names
#' @importFrom DT renderDataTable datatable formatStyle styleEqual
mod_lipid_norm_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        logger::log_info("=== METABOLOMICS NORMALIZATION MODULE STARTED ===")
        logger::log_info(paste("Module ID:", id))

        # ================================================================
        # Local Reactive Values
        # ================================================================
        norm_data <- shiny::reactiveValues(
            # --- Assay information ---
            assay_names = NULL
            , itsd_selections = list()  # Per-assay ITSD selections

            # --- Processing stages ---
            , pre_norm_qc_generated = FALSE
            , normalization_complete = FALSE
            , ruv_complete = FALSE
            , correlation_filtering_complete = FALSE

            # --- QC plot file paths (keyed by assay_name and stage) ---
            , qc_plot_paths = list()

            # --- S4 Objects at different stages ---
            , post_filter_obj = NULL
            , post_itsd_obj = NULL
            , post_log2_obj = NULL
            , post_norm_obj = NULL
            , ruv_corrected_obj = NULL
            , correlation_filtered_obj = NULL

            # --- RUV optimization results (per-assay) ---
            , ruv_optimization_results = list()

            # --- Correlation filtering ---
            , correlation_results = list()

            # --- Plot refresh trigger ---
            , plot_refresh_trigger = 0

            # --- Normalization log ---
            , normalization_log = character(0)
        )

        # ================================================================
        # Helper Functions
        # ================================================================

        # Add message to normalization log
        add_log <- function(message) {
            timestamp <- format(Sys.time(), "%H:%M:%S")
            norm_data$normalization_log <- c(
                norm_data$normalization_log
                , sprintf("[%s] %s", timestamp, message)
            )
        }

        # Get current plot aesthetics from input
        getPlotAesthetics <- function() {
            list(
                color_var = if (is.null(input$color_variable) || input$color_variable == "") {
                    "group"
                } else {
                    input$color_variable
                }
                , shape_var = if (is.null(input$shape_variable) || input$shape_variable == "") {
                    "group"
                } else {
                    input$shape_variable
                }
            )
        }

        # Generate composite QC figure from saved images
        # Mirrors createGridQC structure with row labels and patchwork layout
        # @param plot_files Vector of file paths to individual QC plot images
        # @param output_path Path where the composite figure will be saved
        # @param ncol Number of columns in the composite (default: 3)
        # @param row_labels Named list of row labels (e.g., list(pca = c("a)", "b)", "c)")))
        # @param column_labels Vector of column titles (e.g., c("Pre-Norm", "Post-Norm", "RUV"))
        # @return Path to saved composite or NULL on error
        generateCompositeFromFiles <- function(plot_files, output_path, ncol = 3, row_labels = NULL, column_labels = NULL) {
            logger::log_info(sprintf("[generateCompositeFromFiles] Generating composite from %d files...", length(plot_files)))

            if (!requireNamespace("patchwork", quietly = TRUE)) {
                warning("patchwork package required for composite generation")
                return(NULL)
            }
            if (!requireNamespace("ggplot2", quietly = TRUE)) {
                warning("ggplot2 package required for composite generation")
                return(NULL)
            }
            if (!requireNamespace("png", quietly = TRUE)) {
                warning("png package required for composite generation")
                return(NULL)
            }

            # --- Helper: Create label plot ---
            createLabelPlot <- function(title) {
                ggplot2::ggplot() +
                    ggplot2::annotate("text", x = 0, y = 0.5, label = title, size = 5, hjust = 0) +
                    ggplot2::xlim(0, 1) +
                    ggplot2::theme_void() +
                    ggplot2::theme(
                        plot.margin = ggplot2::margin(5, 5, 5, 5)
                        , panel.background = ggplot2::element_blank()
                    )
            }

            # --- Helper: Create column title plot ---
            createTitlePlot <- function(title) {
                ggplot2::ggplot() +
                    ggplot2::annotate("text", x = 0.5, y = 0.5, label = title, size = 6, fontface = "bold", hjust = 0.5) +
                    ggplot2::xlim(0, 1) +
                    ggplot2::theme_void() +
                    ggplot2::theme(
                        plot.margin = ggplot2::margin(5, 5, 10, 5)
                        , panel.background = ggplot2::element_blank()
                    )
            }

            # --- Helper: Load image as ggplot ---
            loadImageAsPlot <- function(file_path) {
                if (is.na(file_path) || !file.exists(file_path)) {
                    return(ggplot2::ggplot() + ggplot2::theme_void())
                }
                tryCatch({
                    img <- png::readPNG(file_path)
                    g <- grid::rasterGrob(img, interpolate = TRUE)
                    ggplot2::ggplot() +
                        ggplot2::annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
                        ggplot2::theme_void()
                }, error = function(e) {
                    logger::log_warn(sprintf("[generateCompositeFromFiles] Could not load image: %s", file_path))
                    ggplot2::ggplot() + ggplot2::theme_void()
                })
            }

            tryCatch({
                n_files <- length(plot_files)
                n_plot_types <- n_files / ncol

                # Default row labels if not provided
                if (is.null(row_labels)) {
                    all_labels <- letters[1:n_files]
                    row_labels <- split(paste0(all_labels, ")"), rep(1:n_plot_types, each = ncol))
                    names(row_labels) <- paste0("row", seq_len(n_plot_types))
                }

                plot_sections <- list()
                height_values <- c()

                # Add Column Titles if provided
                if (!is.null(column_labels)) {
                    if (length(column_labels) == ncol) {
                        title_plots <- lapply(column_labels, createTitlePlot)
                        plot_sections <- append(plot_sections, list(
                            patchwork::wrap_plots(title_plots, ncol = ncol)
                        ))
                        height_values <- c(height_values, 0.2)
                        logger::log_info("[generateCompositeFromFiles] Added column titles")
                    }
                }

                row_names <- names(row_labels)

                for (i in seq_along(row_names)) {
                    row_name <- row_names[i]
                    labels <- row_labels[[row_name]]

                    start_idx <- (i - 1) * ncol + 1
                    end_idx <- min(i * ncol, n_files)
                    row_files <- plot_files[start_idx:end_idx]

                    has_files <- any(!is.na(row_files) & sapply(row_files, function(f) !is.na(f) && file.exists(f)))

                    if (has_files) {
                        label_plots <- lapply(labels, createLabelPlot)
                        image_plots <- lapply(row_files, loadImageAsPlot)

                        plot_sections <- append(plot_sections, list(
                            patchwork::wrap_plots(label_plots, ncol = ncol)
                            , patchwork::wrap_plots(image_plots, ncol = ncol)
                        ))
                        height_values <- c(height_values, 0.1, 1)

                        logger::log_info(sprintf("[generateCompositeFromFiles] Added row: %s", row_name))
                    } else {
                        logger::log_info(sprintf("[generateCompositeFromFiles] Skipping empty row: %s", row_name))
                    }
                }

                if (length(plot_sections) == 0) {
                    warning("No valid plot sections to combine")
                    return(NULL)
                }

                logger::log_info("[generateCompositeFromFiles] Combining plot sections...")
                combined_plot <- patchwork::wrap_plots(plot_sections, ncol = 1) +
                    patchwork::plot_layout(heights = height_values)

                plot_width <- 4 + (ncol * 3)
                plot_height <- 4 + (length(height_values) * 2)

                ggplot2::ggsave(
                    output_path
                    , combined_plot
                    , width = plot_width
                    , height = plot_height
                    , dpi = 150
                    , limitsize = FALSE
                )

                rm(plot_sections, combined_plot)
                gc()

                logger::log_info(sprintf("[generateCompositeFromFiles] Composite saved to %s", output_path))
                return(output_path)

            }, error = function(e) {
                logger::log_error(paste("[generateCompositeFromFiles] Error:", e$message))
                return(NULL)
            })
        }

        # Helper function to generate pre-normalization QC plots
        generatePreNormalizationQc <- function() {
            logger::log_info("=== GENERATING PRE-NORMALIZATION QC PLOTS ===")
            
            # Get current S4 object
            shiny::req(workflow_data$state_manager)
            current_s4 <- workflow_data$state_manager$getState()
            
            if (is.null(current_s4)) {
                logger::log_warn("No S4 object available for QC plot generation")
                return()
            }
            
            # Ensure assay names are populated for renderImage functions
            if (inherits(current_s4, "LipidomicsAssayData")) {
                detected_assays <- names(current_s4@lipid_data)
                if (length(detected_assays) > 0 && is.null(norm_data$assay_names)) {
                    norm_data$assay_names <- detected_assays
                    logger::log_info(paste("Set assay names:", paste(detected_assays, collapse = ", ")))
                }
            }
            
            # Get aesthetics
            aesthetics <- getPlotAesthetics()
            
            # Generate plots
            tryCatch({
                generateMetabQcPlots(
                    theObject = current_s4
                    , experiment_paths = experiment_paths
                    , stage = "post_filter"
                    , grouping_variable = aesthetics$color_var
                    , shape_variable = aesthetics$shape_var
                )
                
                # Trigger plot refresh
                norm_data$plot_refresh_trigger <- norm_data$plot_refresh_trigger + 1
                norm_data$pre_norm_qc_generated <- TRUE
                logger::log_info("Pre-normalization QC plots generated successfully")
                
            }, error = function(e) {
                logger::log_error(paste("Error generating pre-normalization QC:", e$message))
                add_log(paste("Error generating Pre-QC:", e$message))
            })
        }

        # Render QC image from disk for a specific assay slot
        # Uses reactive filename resolution based on assay position (1 or 2)
        # @param assay_slot Integer (1 or 2) indicating which assay row
        # @param plot_type Character: "pca", "density", "rle", or "correlation"
        # @param stage_prefix Character: "pre_norm", "post_norm", or "ruv_corrected"
        render_qc_image_for_assay <- function(assay_slot, plot_type, stage_prefix) {
            shiny::renderImage({
                # Take dependency on refresh trigger AND assay names
                norm_data$plot_refresh_trigger
                assay_names <- norm_data$assay_names
                
                # Resolve actual assay name from slot
                if (is.null(assay_names) || length(assay_names) < assay_slot) {
                    return(list(src = "", alt = "Assay not detected yet"))
                }
                
                assay_name <- assay_names[[assay_slot]]
                safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
                filename <- paste0(safe_name, "_", stage_prefix, "_", plot_type, ".png")
                
                qc_dir <- experiment_paths$lipid_qc_dir
                if (is.null(qc_dir)) {
                    return(list(src = "", alt = "QC directory not configured"))
                }
                
                img_path <- file.path(qc_dir, filename)
                
                if (file.exists(img_path)) {
                    list(
                        src = img_path
                        , contentType = "image/png"
                        , width = "100%"
                        , height = "auto"
                        , alt = paste(plot_type, "-", assay_name)
                    )
                } else {
                    list(src = "", alt = paste("Plot not generated yet:", filename))
                }
            }, deleteFile = FALSE)
        }

        # ================================================================
        # Initialize Assay Names from S4 Object
        # ================================================================
        shiny::observe({
            shiny::req(workflow_data$state_manager)

            tryCatch({
                current_s4 <- workflow_data$state_manager$getState()
                
                if (inherits(current_s4, "LipidomicsAssayData")) {
                    detected_assays <- names(current_s4@lipid_data)
                    norm_data$assay_names <- detected_assays
                    logger::log_info(paste("Detected assays:", paste(detected_assays, collapse = ", ")))

                    # Initialize ITSD selections for each assay
                    for (assay_name in norm_data$assay_names) {
                        if (!assay_name %in% names(norm_data$itsd_selections)) {
                            norm_data$itsd_selections[[assay_name]] <- NULL
                        }
                    }
                }
            }, error = function(e) {
                logger::log_warn(paste("Could not detect assay names:", e$message))
            })
        })

        # ================================================================
        # Auto-trigger Pre-Normalization QC when tab is selected
        # ================================================================
        if (!is.null(selected_tab)) {
            shiny::observeEvent(selected_tab(), {
                # Only trigger if normalization tab is selected ("norm" in mod_lipidomics.R)
                if (!is.null(selected_tab()) && selected_tab() == "norm") {
                    logger::log_info("Normalization tab selected - checking if pre-QC needed")
                    
                    # Only trigger if we haven't run pre-norm QC yet and have valid data
                    if (!norm_data$pre_norm_qc_generated) {
                        shiny::req(workflow_data$state_manager)
                        current_s4 <- tryCatch(workflow_data$state_manager$getState(), error = function(e) NULL)
                        
                        if (!is.null(current_s4) && inherits(current_s4, "LipidomicsAssayData")) {
                            logger::log_info("Auto-triggering pre-normalization QC plots")
                            
                            shiny::withProgress(
                                message = "Generating Pre-Normalization QC..."
                                , value = 0.5
                                , {
                                    generatePreNormalizationQc()
                                }
                            )
                            norm_data$pre_norm_qc_generated <- TRUE
                        }
                    }
                }
            }, ignoreInit = FALSE)
        }

        # ================================================================
        # Update Plot Aesthetic Choices from Design Matrix
        # ================================================================
        shiny::observe({
            if (!is.null(workflow_data$design_matrix)) {
                design_cols <- colnames(workflow_data$design_matrix)

                plot_available_vars <- intersect(
                    design_cols
                    , c("group", "factor1", "factor2", "batch", "technical_replicate_id", "sample_id")
                )

                if (length(plot_available_vars) > 0) {
                    shiny::updateSelectInput(
                        session, "color_variable"
                        , choices = plot_available_vars
                        , selected = if ("group" %in% plot_available_vars) "group" else plot_available_vars[1]
                    )
                    shiny::updateSelectInput(
                        session, "shape_variable"
                        , choices = plot_available_vars
                        , selected = if ("group" %in% plot_available_vars) "group" else plot_available_vars[1]
                    )
                }

                ruv_available_vars <- intersect(design_cols, c("group", "factor1", "factor2", "batch"))
                if (length(ruv_available_vars) > 0) {
                    shiny::updateSelectInput(
                        session, "ruv_grouping_variable"
                        , choices = ruv_available_vars
                        , selected = if ("group" %in% ruv_available_vars) "group" else ruv_available_vars[1]
                    )
                }
            }
        })

        # ================================================================
        # Render Normalization Log
        # ================================================================
        output$norm_log <- shiny::renderText({
            if (length(norm_data$normalization_log) == 0) {
                return("Normalization log will appear here as you apply steps...")
            }
            paste(norm_data$normalization_log, collapse = "\n")
        })

        # ================================================================
        # Dynamic UI: ITSD Selection Tables (per-assay)
        # ================================================================
        output$itsd_selection_ui <- shiny::renderUI({
            shiny::req(norm_data$assay_names)

            assay_uis <- purrr::map(norm_data$assay_names, \(assay_name) {
                safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
                shiny::wellPanel(
                    shiny::h5(paste("Assay:", assay_name))
                    , DT::dataTableOutput(ns(paste0("itsd_table_", safe_name)))
                    , shiny::br()
                )
            })

            shiny::tagList(assay_uis)
        })

        # ================================================================
        # Dynamic UI: RUV QC Plots (per-assay stacked)
        # ================================================================
        output$ruv_qc_ui <- shiny::renderUI({
            shiny::req(norm_data$assay_names)

            assay_uis <- purrr::map(norm_data$assay_names, \(assay_name) {
                safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
                shiny::tagList(
                    shiny::fluidRow(
                        shiny::column(12
                            , shiny::h5(paste("Assay:", assay_name), style = "border-bottom: 1px solid #ddd; padding-bottom: 5px;")
                        )
                    )
                    , shiny::fluidRow(
                        shiny::column(8
                            , shinyjqui::jqui_resizable(
                                shiny::plotOutput(ns(paste0("cancor_plot_", safe_name)), height = "400px")
                            )
                        )
                        , shiny::column(4
                            , shiny::wellPanel(
                                shiny::h6("Optimization Summary")
                                , shiny::verbatimTextOutput(ns(paste0("ruv_summary_", safe_name)))
                                , shiny::br()
                                , shiny::h6("Results Table")
                                , DT::dataTableOutput(ns(paste0("ruv_table_", safe_name)))
                            )
                        )
                    )
                    , shiny::hr()
                )
            })

            shiny::tagList(assay_uis)
        })

        # ================================================================
        # STATIC Output Bindings: QC Plot Images (24 total)
        # Bound at server startup - no race condition
        # ================================================================
        
        # --- PCA outputs (6) ---
        output$pca_post_filter_assay1 <- render_qc_image_for_assay(1, "pca", "pre_norm")
        output$pca_post_norm_assay1 <- render_qc_image_for_assay(1, "pca", "post_norm")
        output$pca_ruv_corrected_assay1 <- render_qc_image_for_assay(1, "pca", "ruv_corrected")
        output$pca_post_filter_assay2 <- render_qc_image_for_assay(2, "pca", "pre_norm")
        output$pca_post_norm_assay2 <- render_qc_image_for_assay(2, "pca", "post_norm")
        output$pca_ruv_corrected_assay2 <- render_qc_image_for_assay(2, "pca", "ruv_corrected")
        
        # --- Density outputs (6) ---
        output$density_post_filter_assay1 <- render_qc_image_for_assay(1, "density", "pre_norm")
        output$density_post_norm_assay1 <- render_qc_image_for_assay(1, "density", "post_norm")
        output$density_ruv_corrected_assay1 <- render_qc_image_for_assay(1, "density", "ruv_corrected")
        output$density_post_filter_assay2 <- render_qc_image_for_assay(2, "density", "pre_norm")
        output$density_post_norm_assay2 <- render_qc_image_for_assay(2, "density", "post_norm")
        output$density_ruv_corrected_assay2 <- render_qc_image_for_assay(2, "density", "ruv_corrected")
        
        # --- RLE outputs (6) ---
        output$rle_post_filter_assay1 <- render_qc_image_for_assay(1, "rle", "pre_norm")
        output$rle_post_norm_assay1 <- render_qc_image_for_assay(1, "rle", "post_norm")
        output$rle_ruv_corrected_assay1 <- render_qc_image_for_assay(1, "rle", "ruv_corrected")
        output$rle_post_filter_assay2 <- render_qc_image_for_assay(2, "rle", "pre_norm")
        output$rle_post_norm_assay2 <- render_qc_image_for_assay(2, "rle", "post_norm")
        output$rle_ruv_corrected_assay2 <- render_qc_image_for_assay(2, "rle", "ruv_corrected")
        
        # --- Correlation outputs (6) ---
        output$correlation_post_filter_assay1 <- render_qc_image_for_assay(1, "correlation", "pre_norm")
        output$correlation_post_norm_assay1 <- render_qc_image_for_assay(1, "correlation", "post_norm")
        output$correlation_ruv_corrected_assay1 <- render_qc_image_for_assay(1, "correlation", "ruv_corrected")
        output$correlation_post_filter_assay2 <- render_qc_image_for_assay(2, "correlation", "pre_norm")
        output$correlation_post_norm_assay2 <- render_qc_image_for_assay(2, "correlation", "post_norm")
        output$correlation_ruv_corrected_assay2 <- render_qc_image_for_assay(2, "correlation", "ruv_corrected")
        
        # ================================================================
        # STATIC Output Bindings: Assay Labels (8 total)
        # ================================================================
        
        # Helper for assay label rendering
        render_assay_label <- function(assay_slot) {
            shiny::renderText({
                if (!is.null(norm_data$assay_names) && length(norm_data$assay_names) >= assay_slot) {
                    paste("Assay:", norm_data$assay_names[[assay_slot]])
                } else {
                    paste0("Assay ", assay_slot, ": (detecting...)")
                }
            })
        }
        
        # PCA tab labels
        output$assay1_label_pca <- render_assay_label(1)
        output$assay2_label_pca <- render_assay_label(2)
        
        # Density tab labels
        output$assay1_label_density <- render_assay_label(1)
        output$assay2_label_density <- render_assay_label(2)
        
        # RLE tab labels
        output$assay1_label_rle <- render_assay_label(1)
        output$assay2_label_rle <- render_assay_label(2)
        
        # Correlation tab labels
        output$assay1_label_correlation <- render_assay_label(1)
        output$assay2_label_correlation <- render_assay_label(2)

        # ================================================================
        # Render ITSD Selection Tables (per-assay)
        # ================================================================
        shiny::observe({
            shiny::req(norm_data$assay_names)
            shiny::req(workflow_data$state_manager)

            current_s4 <- tryCatch({
                workflow_data$state_manager$getState()
            }, error = function(e) NULL)

            if (is.null(current_s4) || !inherits(current_s4, "LipidomicsAssayData")) {
                return()
            }

            purrr::walk(norm_data$assay_names, \(assay_name) {
                safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
                output_id <- paste0("itsd_table_", safe_name)

                output[[output_id]] <- DT::renderDataTable({
                    assay_data <- current_s4@lipid_data[[assay_name]]
                    if (is.null(assay_data)) return(NULL)

                    lipid_id_col <- current_s4@lipid_id_column
                    annotation_col <- current_s4@annotation_id_column

                    selection_table <- buildItsdSelectionTable(
                        assay_data = assay_data
                        , lipid_id_col = lipid_id_col
                        , annotation_cols = annotation_col
                    )

                    # Pre-select candidates
                    preselected <- which(selection_table$is_candidate)

                    DT::datatable(
                        selection_table
                        , selection = list(mode = "multiple", selected = preselected)
                        , filter = "top"
                        , options = list(
                            pageLength = 10
                            , scrollX = TRUE
                            , order = list(list(4, "desc"), list(3, "asc"))
                        )
                        , rownames = FALSE
                    ) |>
                        DT::formatStyle(
                            "is_candidate"
                            , backgroundColor = DT::styleEqual(TRUE, "#d4edda")
                        ) |>
                        DT::formatRound(columns = c("mean_intensity", "cv_percent"), digits = 2)
                })
            })
        })

        # ================================================================
        # Track ITSD Selections from DT
        # ================================================================
        shiny::observe({
            shiny::req(norm_data$assay_names)

            purrr::walk(norm_data$assay_names, \(assay_name) {
                safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
                input_id <- paste0("itsd_table_", safe_name, "_rows_selected")

                shiny::observeEvent(input[[input_id]], {
                    norm_data$itsd_selections[[assay_name]] <- input[[input_id]]
                    logger::log_info(paste(
                        "ITSD selection updated for", assay_name, ":"
                        , length(input[[input_id]]), "features selected"
                    ))
                }, ignoreNULL = FALSE)
            })
        })

        # ================================================================
        # Main Normalization Pipeline
        # ================================================================
        shiny::observeEvent(input$run_normalization, {
            shiny::req(workflow_data$state_manager)

            add_log("Starting normalization pipeline...")

            shiny::withProgress(
                message = "Running normalization pipeline..."
                , value = 0
                , {
                    tryCatch({
                        current_s4 <- workflow_data$state_manager$getState()
                        shiny::req(current_s4)

                        aesthetics <- getPlotAesthetics()
                        total_steps <- 6  # Pre-QC, ITSD, Log2, Norm, RUV, Post-QC

                        # --- Store Post-Filtering state ---
                        shiny::incProgress(1/total_steps, detail = "Capturing pre-normalization state...")
                        norm_data$post_filter_obj <- current_s4
                        add_log("Post-filtering state captured")

                        # --- Generate Pre-Norm QC Plots ---
                        generateMetabQcPlots(
                            theObject = current_s4
                            , experiment_paths = experiment_paths
                            , stage = "post_filter"
                            , grouping_variable = aesthetics$color_var
                            , shape_variable = aesthetics$shape_var
                        )
                        add_log("Pre-normalization QC plots generated")

                        # --- Step 1: ITSD Normalization ---
                        shiny::incProgress(1/total_steps, detail = "Applying ITSD normalization...")
                        if (isTRUE(input$apply_itsd)) {
                            add_log(paste("Applying ITSD normalization (aggregation:", input$itsd_aggregation, ")"))

                            # Convert DT row selections to feature IDs per assay
                            itsd_feature_ids <- NULL
                            has_manual_selections <- any(sapply(norm_data$itsd_selections, \(x) length(x) > 0))

                            if (has_manual_selections) {
                                lipid_id_col <- current_s4@lipid_id_column
                                annotation_col <- current_s4@annotation_id_column
                                itsd_feature_ids <- purrr::imap(norm_data$itsd_selections, \(row_indices, assay_name) {
                                    if (is.null(row_indices) || length(row_indices) == 0) {
                                        return(NULL)  # No manual selection for this assay
                                    }
                                    assay_data <- current_s4@lipid_data[[assay_name]]
                                    if (is.null(assay_data)) return(NULL)

                                    # Rebuild selection table to get feature IDs
                                    selection_table <- buildItsdSelectionTable(
                                        assay_data = assay_data
                                        , lipid_id_col = lipid_id_col
                                        , annotation_cols = annotation_col
                                    )

                                    # Extract feature IDs for selected rows
                                    selected_ids <- selection_table$feature_id[row_indices]
                                    add_log(paste("Assay", assay_name, ":", length(selected_ids), "ITSD features selected"))
                                    return(selected_ids)
                                })
                                # Remove NULL entries (assays without manual selections)
                                itsd_feature_ids <- purrr::compact(itsd_feature_ids)
                                if (length(itsd_feature_ids) == 0) itsd_feature_ids <- NULL
                            }

                            current_s4 <- normaliseUntransformedData(
                                theObject = current_s4
                                , method = "ITSD"
                                , itsd_aggregation = input$itsd_aggregation
                                , itsd_feature_ids = itsd_feature_ids
                            )
                            norm_data$post_itsd_obj <- current_s4

                            workflow_data$state_manager$saveState(
                                state_name = "lipid_itsd_norm"
                                , s4_data_object = current_s4
                                , config_object = workflow_data$config_list
                                , description = paste("ITSD normalization (aggregation:", input$itsd_aggregation, ")")
                            )
                            add_log("ITSD normalization complete")
                        } else {
                            add_log("ITSD normalization skipped")
                        }

                        # --- Step 2: Log2 Transform ---
                        shiny::incProgress(1/total_steps, detail = "Applying log2 transformation...")
                        add_log(paste("Applying log2 transformation (offset:", input$log_offset, ")"))

                        current_s4 <- logTransformAssays(
                            theObject = current_s4
                            , offset = input$log_offset
                        )
                        norm_data$post_log2_obj <- current_s4

                        workflow_data$state_manager$saveState(
                            state_name = "lipid_log2"
                            , s4_data_object = current_s4
                            , config_object = workflow_data$config_list
                            , description = paste("Log2 transformation (offset:", input$log_offset, ")")
                        )
                        add_log("Log2 transformation complete")

                        # --- Step 3: Between-Sample Normalization ---
                        shiny::incProgress(1/total_steps, detail = "Applying between-sample normalization...")
                        if (input$norm_method != "none") {
                            add_log(paste("Applying between-sample normalization (method:", input$norm_method, ")"))

                            current_s4 <- normaliseBetweenSamples(
                                theObject = current_s4
                                , normalisation_method = input$norm_method
                            )
                        }
                        norm_data$post_norm_obj <- current_s4

                        workflow_data$state_manager$saveState(
                            state_name = "lipid_normalized"
                            , s4_data_object = current_s4
                            , config_object = workflow_data$config_list
                            , description = paste("Between-sample normalization (method:", input$norm_method, ")")
                        )
                        add_log("Between-sample normalization complete")
                        norm_data$normalization_complete <- TRUE

                        # --- Generate Post-Norm QC Plots ---
                        generateMetabQcPlots(
                            theObject = current_s4
                            , experiment_paths = experiment_paths
                            , stage = "post_norm"
                            , grouping_variable = aesthetics$color_var
                            , shape_variable = aesthetics$shape_var
                        )
                        add_log("Post-normalization QC plots generated")

                        # --- Step 4: RUV-III Batch Correction ---
                        shiny::incProgress(1/total_steps, detail = "Running RUV-III batch correction...")
                        if (input$ruv_mode != "skip") {
                            add_log(paste("Running RUV-III (mode:", input$ruv_mode, ")"))

                            ruv_params <- list(
                                percentage_min = input$auto_percentage_min
                                , percentage_max = input$auto_percentage_max
                                , ruv_grouping_variable = input$ruv_grouping_variable
                                , separation_metric = input$separation_metric
                                , k_penalty_weight = input$k_penalty_weight
                                , adaptive_k_penalty = input$adaptive_k_penalty
                                , manual_k = input$ruv_k
                                , manual_percentage = input$ruv_percentage
                            )

                            ruv_results <- runPerAssayRuvOptimization(
                                theObject = current_s4
                                , ruv_mode = input$ruv_mode
                                , params = ruv_params
                                , experiment_paths = experiment_paths
                            )
                            norm_data$ruv_optimization_results <- ruv_results

                            # Extract best k and ctrl per assay
                            best_k_list <- extractBestKPerAssay(ruv_results)
                            ctrl_list <- extractCtrlPerAssay(ruv_results)

                            add_log(paste("RUV optimization complete. Best k per assay:",
                                paste(names(best_k_list), "=", unlist(best_k_list), collapse = ", ")
                            ))

                            # Apply RUV-III-C with per-assay parameters
                            current_s4 <- ruvIII_C_Varying(
                                theObject = current_s4
                                , ruv_grouping_variable = input$ruv_grouping_variable
                                , ruv_number_k = best_k_list
                                , ctrl = ctrl_list
                            )
                            norm_data$ruv_corrected_obj <- current_s4
                            norm_data$ruv_complete <- TRUE

                            workflow_data$state_manager$saveState(
                                state_name = "lipid_ruv_corrected"
                                , s4_data_object = current_s4
                                , config_object = workflow_data$config_list
                                , description = "RUV-III batch correction complete"
                            )
                            add_log("RUV-III correction applied")

                            # --- Generate RUV QC Plots ---
                            shiny::incProgress(1/total_steps, detail = "Generating RUV QC plots...")
                            generateMetabQcPlots(
                                theObject = current_s4
                                , experiment_paths = experiment_paths
                                , stage = "ruv_corrected"
                                , grouping_variable = aesthetics$color_var
                                , shape_variable = aesthetics$shape_var
                            )
                            add_log("RUV QC plots generated")

                        } else {
                            add_log("RUV-III skipped")
                            norm_data$ruv_corrected_obj <- current_s4
                            norm_data$ruv_complete <- TRUE
                        }

                        # --- Generate Composite QC Figure ---
                        add_log("Generating composite QC figure...")
                        tryCatch({
                            qc_dir <- experiment_paths$lipid_qc_dir
                            if (!is.null(qc_dir) && dir.exists(qc_dir) && !is.null(norm_data$assay_names)) {

                                # Determine columns based on RUV mode
                                ncol_composite <- if (input$ruv_mode == "skip") 2 else 3
                                column_labels <- if (input$ruv_mode == "skip") {
                                    c("Pre-Normalisation", "Post-Normalisation")
                                } else {
                                    c("Pre-Normalisation", "Post-Normalisation", "RUV-Corrected")
                                }

                                # Build file list and row labels for all assays combined
                                all_plot_files <- c()
                                all_row_labels <- list()
                                label_counter <- 1

                                for (assay_name in norm_data$assay_names) {
                                    safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))

                                    # Plot types to include
                                    plot_types <- c("pca", "density", "rle", "correlation")

                                    for (plot_type in plot_types) {
                                        if (input$ruv_mode == "skip") {
                                            # 2 columns: Pre-Norm, Post-Norm
                                            files <- c(
                                                file.path(qc_dir, sprintf("%s_pre_norm_%s.png", safe_name, plot_type))
                                                , file.path(qc_dir, sprintf("%s_post_norm_%s.png", safe_name, plot_type))
                                            )
                                            labels <- c(
                                                sprintf("%s)", letters[label_counter])
                                                , sprintf("%s)", letters[label_counter + 1])
                                            )
                                            label_counter <- label_counter + 2
                                        } else {
                                            # 3 columns: Pre-Norm, Post-Norm, RUV-Corrected
                                            files <- c(
                                                file.path(qc_dir, sprintf("%s_pre_norm_%s.png", safe_name, plot_type))
                                                , file.path(qc_dir, sprintf("%s_post_norm_%s.png", safe_name, plot_type))
                                                , file.path(qc_dir, sprintf("%s_ruv_corrected_%s.png", safe_name, plot_type))
                                            )
                                            labels <- c(
                                                sprintf("%s)", letters[label_counter])
                                                , sprintf("%s)", letters[label_counter + 1])
                                                , sprintf("%s)", letters[label_counter + 2])
                                            )
                                            label_counter <- label_counter + 3
                                        }

                                        all_plot_files <- c(all_plot_files, files)
                                        row_key <- sprintf("%s_%s", safe_name, plot_type)
                                        all_row_labels[[row_key]] <- labels
                                    }
                                }

                                # Generate combined composite
                                composite_path <- file.path(qc_dir, "composite_QC_figure.png")
                                generateCompositeFromFiles(
                                    plot_files = all_plot_files
                                    , output_path = composite_path
                                    , ncol = ncol_composite
                                    , row_labels = all_row_labels
                                    , column_labels = column_labels
                                )
                                add_log(sprintf("Composite QC figure saved to: %s", composite_path))

                            }
                        }, error = function(e) {
                            add_log(paste("Warning: Could not generate composite QC figure:", e$message))
                            logger::log_warn(paste("Composite QC generation failed:", e$message))
                        })

                        # --- Trigger plot refresh ---
                        norm_data$plot_refresh_trigger <- norm_data$plot_refresh_trigger + 1

                        add_log("Normalization pipeline complete!")
                        shiny::showNotification("Normalization pipeline complete!", type = "message")

                    }, error = function(e) {
                        add_log(paste("ERROR:", e$message))
                        logger::log_error(paste("Normalization pipeline error:", e$message))
                        shiny::showNotification(paste("Error:", e$message), type = "error")
                    })
                }
            )
        })

        # ================================================================
        # Reset to Pre-Normalization
        # ================================================================
        shiny::observeEvent(input$reset_normalization, {
            shiny::req(workflow_data$state_manager)

            tryCatch({
                # Reset to filtered state
                if (!is.null(norm_data$post_filter_obj)) {
                    workflow_data$state_manager$saveState(
                        state_name = "lipid_reset"
                        , s4_data_object = norm_data$post_filter_obj
                        , config_object = workflow_data$config_list
                        , description = "Reset to pre-normalization state"
                    )
                }

                # Reset local state
                norm_data$normalization_complete <- FALSE
                norm_data$ruv_complete <- FALSE
                norm_data$correlation_filtering_complete <- FALSE
                norm_data$post_norm_obj <- NULL
                norm_data$ruv_corrected_obj <- NULL
                norm_data$correlation_filtered_obj <- NULL
                norm_data$ruv_optimization_results <- list()

                add_log("Reset to pre-normalization state")
                shiny::showNotification("Reset to pre-normalization state", type = "message")

            }, error = function(e) {
                add_log(paste("ERROR during reset:", e$message))
                shiny::showNotification(paste("Error:", e$message), type = "error")
            })
        })

        # ================================================================
        # Render RUV Cancor Plots (per-assay)
        # ================================================================
        shiny::observe({
            shiny::req(norm_data$assay_names)
            shiny::req(length(norm_data$ruv_optimization_results) > 0)

            purrr::walk(norm_data$assay_names, \(assay_name) {
                safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))

                # Render cancor plot
                output[[paste0("cancor_plot_", safe_name)]] <- shiny::renderPlot({
                    result <- norm_data$ruv_optimization_results[[assay_name]]
                    if (!is.null(result) && isTRUE(result$success) && !is.null(result$cancor_plot)) {
                        result$cancor_plot
                    } else {
                        ggplot2::ggplot() +
                            ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data", size = 6) +
                            ggplot2::theme_void()
                    }
                })

                # Render summary
                output[[paste0("ruv_summary_", safe_name)]] <- shiny::renderText({
                    result <- norm_data$ruv_optimization_results[[assay_name]]
                    if (!is.null(result) && isTRUE(result$success)) {
                        paste0(
                            "Best k: ", result$best_k, "\n"
                            , "Best %: ", result$best_percentage, "\n"
                            , "Separation: ", round(result$separation_score, 4), "\n"
                            , "Controls: ", sum(result$control_genes_index, na.rm = TRUE)
                        )
                    } else if (!is.null(result)) {
                        paste0("Failed: ", result$error)
                    } else {
                        "Not yet computed"
                    }
                })

                # Render table
                output[[paste0("ruv_table_", safe_name)]] <- DT::renderDataTable({
                    result <- norm_data$ruv_optimization_results[[assay_name]]
                    if (!is.null(result) && isTRUE(result$success) && !is.null(result$optimization_results)) {
                        DT::datatable(
                            result$optimization_results
                            , options = list(pageLength = 5, dom = "t")
                            , rownames = FALSE
                        )
                    } else {
                        NULL
                    }
                })
            })
        })

        # ================================================================
        # Correlation Filtering
        # ================================================================
        shiny::observeEvent(input$apply_correlation_filter, {
            shiny::req(workflow_data$state_manager)
            shiny::req(norm_data$ruv_complete || norm_data$normalization_complete)

            threshold <- input$min_pearson_correlation_threshold
            add_log(paste("Applying correlation filter (threshold:", threshold, ")"))
            shiny::showNotification("Applying correlation filter...", id = "corr_working", duration = NULL)

            tryCatch({
                current_s4 <- if (!is.null(norm_data$ruv_corrected_obj)) {
                    norm_data$ruv_corrected_obj
                } else {
                    norm_data$post_norm_obj
                }
                shiny::req(current_s4)

                # Calculate correlations per assay
                logger::log_info("Calculating Pearson correlations per sample pair...")
                corr_results <- pearsonCorForSamplePairs(
                    theObject = current_s4
                    , correlation_group = input$ruv_grouping_variable
                )

                norm_data$correlation_results <- corr_results

                # Apply filtering using the new S4 method
                filtered_s4 <- filterSamplesByLipidCorrelationThreshold(
                    theObject = current_s4
                    , pearson_correlation_per_pair = corr_results
                    , min_pearson_correlation_threshold = threshold
                )

                norm_data$correlation_filtered_obj <- filtered_s4
                norm_data$correlation_filtering_complete <- TRUE

                workflow_data$state_manager$saveState(
                    state_name = "lipid_correlation_filtered"
                    , s4_data_object = filtered_s4
                    , config_object = workflow_data$config_list
                    , description = paste("Correlation filtering (threshold:", threshold, ")")
                )

                # Ensure QC is marked complete when normalization completes
                # Must replace entire list to trigger reactivity
                updated_status <- workflow_data$tab_status
                updated_status$quality_control <- "complete"
                updated_status$normalization <- "complete"
                workflow_data$tab_status <- updated_status

                add_log("Correlation filtering complete")
                shiny::removeNotification("corr_working")
                shiny::showNotification("Correlation filtering complete! Ready for DE analysis.", type = "message")

            }, error = function(e) {
                add_log(paste("ERROR in correlation filtering:", e$message))
                logger::log_error(paste("Correlation filtering error:", e$message))
                shiny::removeNotification("corr_working")
                shiny::showNotification(paste("Error:", e$message), type = "error")
            })
        })

        # ================================================================
        # Skip Correlation Filtering
        # ================================================================
        shiny::observeEvent(input$skip_correlation_filter, {
            shiny::req(workflow_data$state_manager)
            shiny::req(norm_data$ruv_complete || norm_data$normalization_complete)

            current_s4 <- if (!is.null(norm_data$ruv_corrected_obj)) {
                norm_data$ruv_corrected_obj
            } else {
                norm_data$post_norm_obj
            }

            if (!is.null(current_s4)) {
                workflow_data$state_manager$saveState(
                    state_name = "lipid_norm_complete"
                    , s4_data_object = current_s4
                    , config_object = workflow_data$config_list
                    , description = "Normalization complete (correlation filtering skipped)"
                )

                # Ensure QC is marked complete when normalization completes
                # Must replace entire list to trigger reactivity
                updated_status <- workflow_data$tab_status
                updated_status$quality_control <- "complete"
                updated_status$normalization <- "complete"
                workflow_data$tab_status <- updated_status

                add_log("Correlation filtering skipped - ready for DE analysis")
                shiny::showNotification("Normalization complete! Proceeding to DE analysis.", type = "message")
            }
        })

        # ================================================================
        # Correlation Filter Summary
        # ================================================================
        output$correlation_filter_summary <- shiny::renderText({
            if (!norm_data$correlation_filtering_complete) {
                return("Apply correlation filter to see results...")
            }

            # Generate summary from correlation_results
            corr_results <- norm_data$correlation_results
            filtered_obj <- norm_data$correlation_filtered_obj
            original_obj <- if (!is.null(norm_data$ruv_corrected_obj)) {
                norm_data$ruv_corrected_obj
            } else {
                norm_data$post_norm_obj
            }
            
            if (is.null(corr_results) || length(corr_results) == 0) {
                return("No correlation results available.")
            }
            
            # Build summary text
            summary_lines <- c("=== Correlation Filtering Summary ===\n")
            
            for (assay_name in names(corr_results)) {
                assay_corr <- corr_results[[assay_name]]
                if (!is.null(assay_corr) && nrow(assay_corr) > 0) {
                    n_pairs <- nrow(assay_corr)
                    mean_corr <- round(mean(assay_corr$pearson_correlation, na.rm = TRUE), 3)
                    min_corr <- round(min(assay_corr$pearson_correlation, na.rm = TRUE), 3)
                    max_corr <- round(max(assay_corr$pearson_correlation, na.rm = TRUE), 3)
                    
                    summary_lines <- c(summary_lines, sprintf(
                        "\n[%s]\n  Sample pairs: %d\n  Correlation: mean=%.3f, min=%.3f, max=%.3f",
                        assay_name, n_pairs, mean_corr, min_corr, max_corr
                    ))
                }
            }
            
            # Add sample count comparison if objects are available
            if (!is.null(original_obj) && !is.null(filtered_obj)) {
                original_samples <- nrow(original_obj@design_matrix)
                filtered_samples <- nrow(filtered_obj@design_matrix)
                removed <- original_samples - filtered_samples
                
                summary_lines <- c(summary_lines, sprintf(
                    "\n\n[Sample Filtering]\n  Original: %d samples\n  After filtering: %d samples\n  Removed: %d samples",
                    original_samples, filtered_samples, removed
                ))
            }
            
            paste(summary_lines, collapse = "")
        })

        # ================================================================
        # Final QC Plot
        # ================================================================
        output$final_qc_plot <- shiny::renderPlot({
            shiny::req(norm_data$correlation_filtering_complete || norm_data$ruv_complete)

            current_s4 <- if (!is.null(norm_data$correlation_filtered_obj)) {
                norm_data$correlation_filtered_obj
            } else if (!is.null(norm_data$ruv_corrected_obj)) {
                norm_data$ruv_corrected_obj
            } else {
                norm_data$post_norm_obj
            }

            if (is.null(current_s4)) {
                return(ggplot2::ggplot() +
                    ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data", size = 6) +
                    ggplot2::theme_void())
            }

            aesthetics <- getPlotAesthetics()

            tryCatch({
                pca_plots <- plotPca(
                    current_s4
                    , grouping_variable = aesthetics$color_var
                    , shape_variable = aesthetics$shape_var
                    , title = "Final QC - PCA"
                )

                if (is.list(pca_plots) && length(pca_plots) > 0) {
                    pca_plots[[1]]
                } else if (inherits(pca_plots, "ggplot")) {
                    pca_plots
                } else {
                    ggplot2::ggplot() + ggplot2::theme_void()
                }
            }, error = function(e) {
                ggplot2::ggplot() +
                    ggplot2::annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message), size = 4) +
                    ggplot2::theme_void()
            })
        })

        # ================================================================
        # Export Session (Comprehensive - matches proteomics pattern)
        # ================================================================
        shiny::observeEvent(input$export_session, {
            logger::log_info("=== EXPORT NORMALIZED SESSION BUTTON CLICKED ===")
            shiny::req(workflow_data$state_manager)

            # Check if normalization is complete
            if (!norm_data$normalization_complete) {
                shiny::showNotification(
                    "Please complete normalization before exporting session data."
                    , type = "warning"
                    , duration = 5
                )
                return()
            }

            tryCatch({
                # Get the source directory for saving (consistent with proteomics)
                source_dir <- experiment_paths$source_dir
                if (is.null(source_dir) || !dir.exists(source_dir)) {
                    # Fallback to export_dir if source_dir not available
                    source_dir <- experiment_paths$export_dir
                    if (is.null(source_dir)) {
                        stop("Could not find a valid directory to save session data.")
                    }
                    if (!dir.exists(source_dir)) {
                        dir.create(source_dir, recursive = TRUE)
                    }
                }

                shiny::withProgress(message = "Exporting normalized session data...", value = 0, {

                    # Step 1: Gather all essential data for DE analysis
                    shiny::incProgress(0.2, detail = "Gathering data...")

                    # Get current state from R6 manager
                    current_state_name <- workflow_data$state_manager$current_state
                    current_s4 <- workflow_data$state_manager$getState(current_state_name)

                    # Calculate feature counts per assay
                    feature_counts_per_assay <- NULL
                    if (!is.null(current_s4) && inherits(current_s4, "LipidomicsAssayData")) {
                        feature_counts_per_assay <- purrr::map(names(current_s4@lipid_data), \(assay_name) {
                            assay_data <- current_s4@lipid_data[[assay_name]]
                            if (!is.null(assay_data)) {
                                n_features <- length(unique(assay_data[[current_s4@lipid_id_column]]))
                                n_samples <- ncol(assay_data) - length(c(current_s4@lipid_id_column, current_s4@annotation_id_column))
                                list(features = n_features, samples = n_samples)
                            } else {
                                list(features = 0, samples = 0)
                            }
                        }) |> purrr::set_names(names(current_s4@lipid_data))
                    }

                    # Build comprehensive session data
                    session_data <- list(
                        # --- R6 State Info ---
                        r6_current_state_name = current_state_name
                        , current_s4_object = current_s4

                        # --- Workflow artifacts ---
                        , contrasts_tbl = workflow_data$contrasts_tbl
                        , design_matrix = workflow_data$design_matrix
                        , config_list = workflow_data$config_list

                        # --- Lipidomics-specific ---
                        , itsd_selections = norm_data$itsd_selections
                        , ruv_optimization_results = norm_data$ruv_optimization_results
                        , correlation_results = norm_data$correlation_results
                        , assay_names = norm_data$assay_names

                        # --- Export metadata ---
                        , export_timestamp = Sys.time()
                        , omic_type = "lipidomics"
                        , experiment_label = experiment_label

                        # --- Normalization parameters ---
                        , normalization_method = input$norm_method
                        , ruv_mode = input$ruv_mode
                        , itsd_applied = input$apply_itsd
                        , itsd_aggregation = if (isTRUE(input$apply_itsd)) input$itsd_aggregation else NA
                        , log_offset = input$log_offset
                        , correlation_threshold = input$min_pearson_correlation_threshold
                        , ruv_grouping_variable = input$ruv_grouping_variable

                        # --- Feature counts per assay ---
                        , feature_counts = feature_counts_per_assay
                        , lipid_counts = workflow_data$lipid_counts

                        # --- QC parameters ---
                        , qc_params = workflow_data$qc_params

                        # --- Processing flags ---
                        , normalization_complete = norm_data$normalization_complete
                        , ruv_complete = norm_data$ruv_complete
                        , correlation_filtering_complete = norm_data$correlation_filtering_complete
                    )

                    logger::log_info("*** EXPORT: Gathered session data successfully ***")
                    logger::log_info(sprintf("*** EXPORT: Assays: %s ***", paste(norm_data$assay_names, collapse = ", ")))
                    logger::log_info(sprintf("*** EXPORT: Contrasts available: %d ***", 
                                            ifelse(is.null(session_data$contrasts_tbl), 0, nrow(session_data$contrasts_tbl))))

                    # Step 2: Save to RDS file with timestamp
                    shiny::incProgress(0.3, detail = "Saving to file...")

                    timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
                    session_filename <- sprintf("lipid_filtered_session_data_%s.rds", timestamp_str)
                    session_filepath <- file.path(source_dir, session_filename)

                    saveRDS(session_data, session_filepath)
                    logger::log_info(sprintf("*** EXPORT: Session data saved to: %s ***", session_filepath))

                    # Step 3: Save "latest" version for easy access
                    shiny::incProgress(0.1, detail = "Creating latest version...")

                    latest_filename <- "lipid_filtered_session_data_latest.rds"
                    latest_filepath <- file.path(source_dir, latest_filename)

                    saveRDS(session_data, latest_filepath)
                    logger::log_info(sprintf("*** EXPORT: Latest version saved to: %s ***", latest_filepath))

                    # Step 4: Save individual metadata files for redundancy
                    shiny::incProgress(0.1, detail = "Saving metadata files...")

                    tryCatch({
                        # Save RUV optimization results (per-assay)
                        if (!is.null(session_data$ruv_optimization_results) && length(session_data$ruv_optimization_results) > 0) {
                            ruv_file <- file.path(source_dir, "lipid_ruv_optimization_results.RDS")
                            saveRDS(session_data$ruv_optimization_results, ruv_file)
                            logger::log_info("*** EXPORT: Saved lipid_ruv_optimization_results.RDS ***")
                        }

                        # Save ITSD selections
                        if (!is.null(session_data$itsd_selections) && length(session_data$itsd_selections) > 0) {
                            itsd_file <- file.path(source_dir, "lipid_itsd_selections.RDS")
                            saveRDS(session_data$itsd_selections, itsd_file)
                            logger::log_info("*** EXPORT: Saved lipid_itsd_selections.RDS ***")
                        }

                        # Save QC parameters
                        if (!is.null(session_data$qc_params)) {
                            qc_params_file <- file.path(source_dir, "lipid_qc_params.RDS")
                            saveRDS(session_data$qc_params, qc_params_file)
                            logger::log_info("*** EXPORT: Saved lipid_qc_params.RDS ***")
                        }

                    }, error = function(e) {
                        logger::log_warn(sprintf("*** WARNING: Some metadata files could not be saved: %s ***", e$message))
                    })

                    # Step 5: Create human-readable summary file
                    shiny::incProgress(0.2, detail = "Creating summary...")

                    # Build RUV summary per assay
                    ruv_summary_lines <- ""
                    if (!is.null(session_data$ruv_optimization_results)) {
                        for (assay_name in names(session_data$ruv_optimization_results)) {
                            result <- session_data$ruv_optimization_results[[assay_name]]
                            if (!is.null(result) && isTRUE(result$success)) {
                                ruv_summary_lines <- paste0(ruv_summary_lines, sprintf(
                                    "\n  %s: k=%d, %%=%.1f, controls=%d"
                                    , assay_name
                                    , result$best_k
                                    , result$best_percentage
                                    , sum(result$control_genes_index, na.rm = TRUE)
                                ))
                            }
                        }
                    }
                    if (ruv_summary_lines == "") ruv_summary_lines <- "\n  (RUV skipped or not applied)"

                    # Build feature counts summary
                    feature_summary_lines <- ""
                    if (!is.null(feature_counts_per_assay)) {
                        for (assay_name in names(feature_counts_per_assay)) {
                            counts <- feature_counts_per_assay[[assay_name]]
                            feature_summary_lines <- paste0(feature_summary_lines, sprintf(
                                "\n  %s: %d features, %d samples"
                                , assay_name
                                , counts$features
                                , counts$samples
                            ))
                        }
                    }

                    summary_content <- sprintf(
                        "Lipidomics Normalized Session Data Export Summary\n===================================================\n\nExport Timestamp: %s\nSession File: %s\n\nData Summary:%s\n\nNormalization Parameters:\n- Method: %s\n- ITSD applied: %s\n- ITSD aggregation: %s\n- Log2 offset: %s\n- RUV mode: %s\n- RUV grouping variable: %s\n- Correlation threshold: %s\n\nRUV Optimization Results (per-assay):%s\n\nContrasts:\n%s\n\nThis data is ready for differential expression analysis.\nUse 'Load Filtered Session' in the DE tab to import.\n"
                        , format(Sys.time(), "%Y-%m-%d %H:%M:%S")
                        , session_filename
                        , feature_summary_lines
                        , session_data$normalization_method
                        , ifelse(isTRUE(session_data$itsd_applied), "Yes", "No")
                        , ifelse(is.na(session_data$itsd_aggregation), "N/A", session_data$itsd_aggregation)
                        , session_data$log_offset
                        , session_data$ruv_mode
                        , session_data$ruv_grouping_variable
                        , ifelse(is.null(session_data$correlation_threshold), "N/A", session_data$correlation_threshold)
                        , ruv_summary_lines
                        , if (!is.null(session_data$contrasts_tbl)) paste(session_data$contrasts_tbl$friendly_names, collapse = "\n") else "None defined"
                    )

                    summary_filepath <- file.path(source_dir, "lipid_filtered_session_summary.txt")
                    writeLines(summary_content, summary_filepath)
                    logger::log_info(sprintf("*** EXPORT: Summary saved to: %s ***", summary_filepath))

                })

                add_log(paste("Exported comprehensive session data to:", session_filepath))
                shiny::showNotification(
                    sprintf("Session data exported successfully!\nSaved as: %s\nSee summary file for details.", session_filename)
                    , type = "message"
                    , duration = 10
                )

                logger::log_info("=== EXPORT NORMALIZED SESSION COMPLETED SUCCESSFULLY ===")

            }, error = function(e) {
                logger::log_error(paste("*** ERROR in session export:", e$message, "***"))
                add_log(paste("Export error:", e$message))
                shiny::showNotification(paste("Export error:", e$message), type = "error", duration = 10)
            })
        })
    })
}
