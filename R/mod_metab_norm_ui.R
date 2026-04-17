#' @rdname mod_metab_norm
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 h4 h5 p hr br radioButtons selectInput textInput actionButton uiOutput verbatimTextOutput icon tags checkboxInput numericInput plotOutput conditionalPanel tabsetPanel tabPanel sliderInput helpText
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
                                , "% Metabolites as Negative Controls:"
                                , min = 1
                                , max = 20
                                , value = 5
                                , step = 1
                                , width = "100%"
                            )
                            , shiny::helpText("Percentage of most stable metabolites used as negative controls")

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
                        # Tab 3: PCA (STATIC 2-row x 3-column grid)
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
                        # Tab 4: Density (STATIC 2-row x 3-column grid)
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
                        # Tab 5: RLE (STATIC 2-row x 3-column grid)
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
                        # Tab 6: Correlation Heatmap (STATIC 2-row x 3-column grid)
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

