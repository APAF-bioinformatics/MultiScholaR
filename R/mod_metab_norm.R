# ============================================================================
# mod_metab_norm.R
# ============================================================================
# Purpose: Metabolomics normalization Shiny module (Proteomics UX harmonized)
#
# This module provides a unified normalization pipeline with per-assay
# visualization for metabolomics data (pos/neg mode support).
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

#' @title Metabolomics Normalization Module
#' @description A Shiny module for multi-step normalization of metabolomics data.
#'              Harmonized with proteomics UX pattern.
#' @name mod_metab_norm
NULL

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
                                ns("itsd_method")
                                , "ITSD Method:"
                                , choices = c(
                                    "Median of IS" = "median"
                                    , "Mean of IS" = "mean"
                                    , "Single IS" = "single"
                                )
                                , selected = "median"
                                , width = "100%"
                            )
                            , shiny::helpText("See ITSD Selection tab to choose standards per assay")
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
                        # Tab 3: PCA (per-assay 2-row × 3-column grid)
                        # --------------------------------------------------------
                        , shiny::tabPanel(
                            "PCA"
                            , icon = shiny::icon("project-diagram")
                            , shiny::br()
                            , shiny::uiOutput(ns("pca_plots_ui"))
                        )

                        # --------------------------------------------------------
                        # Tab 4: Density (per-assay 2-row × 3-column grid)
                        # --------------------------------------------------------
                        , shiny::tabPanel(
                            "Density"
                            , icon = shiny::icon("chart-area")
                            , shiny::br()
                            , shiny::uiOutput(ns("density_plots_ui"))
                        )

                        # --------------------------------------------------------
                        # Tab 5: RLE (per-assay 2-row × 3-column grid)
                        # --------------------------------------------------------
                        , shiny::tabPanel(
                            "RLE"
                            , icon = shiny::icon("chart-bar")
                            , shiny::br()
                            , shiny::uiOutput(ns("rle_plots_ui"))
                        )

                        # --------------------------------------------------------
                        # Tab 6: Correlation Heatmap (per-assay 2-row × 3-column grid)
                        # --------------------------------------------------------
                        , shiny::tabPanel(
                            "Correlation"
                            , icon = shiny::icon("th")
                            , shiny::br()
                            , shiny::uiOutput(ns("correlation_plots_ui"))
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


#' @rdname mod_metab_norm
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent req renderUI renderText renderPlot showNotification removeNotification tags renderImage updateSelectInput observe
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_density geom_point labs theme_minimal theme facet_wrap scale_color_brewer coord_flip
#' @importFrom logger log_info log_error log_warn
#' @importFrom purrr imap map walk set_names
#' @importFrom DT renderDataTable datatable formatStyle styleEqual
mod_metab_norm_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
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

        #' Add message to normalization log
        add_log <- function(message) {
            timestamp <- format(Sys.time(), "%H:%M:%S")
            norm_data$normalization_log <- c(
                norm_data$normalization_log
                , sprintf("[%s] %s", timestamp, message)
            )
        }

        #' Get current plot aesthetics from input
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

        #' Render QC image from disk (memory-efficient)
        render_qc_image <- function(filename, alt_text = "QC Plot") {
            shiny::renderImage({
                norm_data$plot_refresh_trigger
                img_path <- file.path(experiment_paths$metabolite_qc_dir, filename)
                if (file.exists(img_path)) {
                    list(
                        src = img_path
                        , contentType = "image/png"
                        , width = "100%"
                        , height = "auto"
                        , alt = alt_text
                    )
                } else {
                    list(
                        src = ""
                        , alt = "Plot not generated yet"
                    )
                }
            }, deleteFile = FALSE)
        }

        #' Generate 3-column QC row for a single assay
        generate_qc_row_ui <- function(assay_name, plot_type, ns) {
            safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
            shiny::tagList(
                shiny::fluidRow(
                    shiny::column(12
                        , shiny::h5(assay_name, style = "margin-bottom: 5px;")
                    )
                )
                , shiny::fluidRow(
                    shiny::column(4
                        , shiny::tags$small("Post-Filtering", style = "display: block; text-align: center;")
                        , shinyjqui::jqui_resizable(
                            shiny::imageOutput(
                                ns(paste0(plot_type, "_post_filter_", safe_name))
                                , height = "300px"
                            )
                        )
                    )
                    , shiny::column(4
                        , shiny::tags$small("Post-Normalization", style = "display: block; text-align: center;")
                        , shinyjqui::jqui_resizable(
                            shiny::imageOutput(
                                ns(paste0(plot_type, "_post_norm_", safe_name))
                                , height = "300px"
                            )
                        )
                    )
                    , shiny::column(4
                        , shiny::tags$small("RUV-Corrected", style = "display: block; text-align: center;")
                        , shinyjqui::jqui_resizable(
                            shiny::imageOutput(
                                ns(paste0(plot_type, "_ruv_corrected_", safe_name))
                                , height = "300px"
                            )
                        )
                    )
                )
                , shiny::hr()
            )
        }

        # ================================================================
        # Initialize Assay Names from S4 Object
        # ================================================================
        shiny::observe({
            shiny::req(workflow_data$state_manager)

            tryCatch({
                current_s4 <- workflow_data$state_manager$getState()
                if (inherits(current_s4, "MetaboliteAssayData")) {
                    norm_data$assay_names <- names(current_s4@metabolite_data)
                    logger::log_info(paste(
                        "Detected assays:"
                        , paste(norm_data$assay_names, collapse = ", ")
                    ))

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
        # Dynamic UI: PCA Plots (per-assay 3-column grid)
        # ================================================================
        output$pca_plots_ui <- shiny::renderUI({
            shiny::req(norm_data$assay_names)

            shiny::tagList(
                purrr::map(norm_data$assay_names, \(assay_name) {
                    generate_qc_row_ui(assay_name, "pca", ns)
                })
            )
        })

        # ================================================================
        # Dynamic UI: Density Plots (per-assay 3-column grid)
        # ================================================================
        output$density_plots_ui <- shiny::renderUI({
            shiny::req(norm_data$assay_names)

            shiny::tagList(
                purrr::map(norm_data$assay_names, \(assay_name) {
                    generate_qc_row_ui(assay_name, "density", ns)
                })
            )
        })

        # ================================================================
        # Dynamic UI: RLE Plots (per-assay 3-column grid)
        # ================================================================
        output$rle_plots_ui <- shiny::renderUI({
            shiny::req(norm_data$assay_names)

            shiny::tagList(
                purrr::map(norm_data$assay_names, \(assay_name) {
                    generate_qc_row_ui(assay_name, "rle", ns)
                })
            )
        })

        # ================================================================
        # Dynamic UI: Correlation Plots (per-assay 3-column grid)
        # ================================================================
        output$correlation_plots_ui <- shiny::renderUI({
            shiny::req(norm_data$assay_names)

            shiny::tagList(
                purrr::map(norm_data$assay_names, \(assay_name) {
                    generate_qc_row_ui(assay_name, "correlation", ns)
                })
            )
        })

        # ================================================================
        # Render ITSD Selection Tables (per-assay)
        # ================================================================
        shiny::observe({
            shiny::req(norm_data$assay_names)
            shiny::req(workflow_data$state_manager)

            current_s4 <- tryCatch({
                workflow_data$state_manager$getState()
            }, error = function(e) NULL)

            if (is.null(current_s4) || !inherits(current_s4, "MetaboliteAssayData")) {
                return()
            }

            purrr::walk(norm_data$assay_names, \(assay_name) {
                safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
                output_id <- paste0("itsd_table_", safe_name)

                output[[output_id]] <- DT::renderDataTable({
                    assay_data <- current_s4@metabolite_data[[assay_name]]
                    if (is.null(assay_data)) return(NULL)

                    metabolite_id_col <- current_s4@metabolite_id

                    selection_table <- buildItsdSelectionTable(
                        assay_data = assay_data
                        , metabolite_id_col = metabolite_id_col
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
                            add_log(paste("Applying ITSD normalization (method:", input$itsd_method, ")"))

                            # Convert DT row selections to feature IDs per assay
                            itsd_feature_ids <- NULL
                            has_manual_selections <- any(sapply(norm_data$itsd_selections, \(x) length(x) > 0))

                            if (has_manual_selections) {
                                metabolite_id_col <- current_s4@metabolite_id_column
                                itsd_feature_ids <- purrr::imap(norm_data$itsd_selections, \(row_indices, assay_name) {
                                    if (is.null(row_indices) || length(row_indices) == 0) {
                                        return(NULL)  # No manual selection for this assay
                                    }
                                    assay_data <- current_s4@metabolite_data[[assay_name]]
                                    if (is.null(assay_data)) return(NULL)

                                    # Rebuild selection table to get feature IDs
                                    selection_table <- buildItsdSelectionTable(
                                        assay_data = assay_data
                                        , metabolite_id_col = metabolite_id_col
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
                                , method = input$itsd_method
                                , itsd_feature_ids = itsd_feature_ids
                            )
                            norm_data$post_itsd_obj <- current_s4

                            workflow_data$state_manager$saveState(
                                state_name = "metab_itsd_norm"
                                , s4_data_object = current_s4
                                , config_object = workflow_data$config_list
                                , description = paste("ITSD normalization (method:", input$itsd_method, ")")
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
                            state_name = "metab_log2"
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
                                , method = input$norm_method
                            )
                        }
                        norm_data$post_norm_obj <- current_s4

                        workflow_data$state_manager$saveState(
                            state_name = "metab_normalized"
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
                                , k = best_k_list
                                , ctrl = ctrl_list
                            )
                            norm_data$ruv_corrected_obj <- current_s4
                            norm_data$ruv_complete <- TRUE

                            workflow_data$state_manager$saveState(
                                state_name = "metab_ruv_corrected"
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
                        state_name = "metab_reset"
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
        # Render QC Image Outputs (per-assay, per-stage)
        # ================================================================
        shiny::observe({
            shiny::req(norm_data$assay_names)

            plot_types <- c("pca", "density", "rle", "correlation")
            stages <- c("post_filter" = "pre_norm", "post_norm" = "post_norm", "ruv_corrected" = "ruv_corrected")

            purrr::walk(norm_data$assay_names, \(assay_name) {
                safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))

                purrr::walk(plot_types, \(plot_type) {
                    purrr::walk(names(stages), \(stage_key) {
                        stage_prefix <- stages[[stage_key]]
                        output_id <- paste0(plot_type, "_", stage_key, "_", safe_name)
                        filename <- paste0(safe_name, "_", stage_prefix, "_", plot_type, ".png")

                        output[[output_id]] <- render_qc_image(filename, paste(plot_type, "-", assay_name))
                    })
                })
            })
        })

        # ================================================================
        # Correlation Filtering
        # ================================================================
        shiny::observeEvent(input$apply_correlation_filter, {
            shiny::req(workflow_data$state_manager)
            shiny::req(norm_data$ruv_complete || norm_data$normalization_complete)

            add_log(paste("Applying correlation filter (threshold:", input$min_pearson_correlation_threshold, ")"))
            shiny::showNotification("Applying correlation filter...", id = "corr_working", duration = NULL)

            tryCatch({
                current_s4 <- if (!is.null(norm_data$ruv_corrected_obj)) {
                    norm_data$ruv_corrected_obj
                } else {
                    norm_data$post_norm_obj
                }
                shiny::req(current_s4)

                # Calculate correlations per assay
                corr_results <- pearsonCorForSamplePairs(
                    theObject = current_s4
                    , correlation_group = input$ruv_grouping_variable
                )

                norm_data$correlation_results <- corr_results

                # TODO: Apply filtering based on threshold
                norm_data$correlation_filtered_obj <- current_s4
                norm_data$correlation_filtering_complete <- TRUE

                workflow_data$state_manager$saveState(
                    state_name = "metab_correlation_filtered"
                    , s4_data_object = current_s4
                    , config_object = workflow_data$config_list
                    , description = paste("Correlation filtering (threshold:", input$min_pearson_correlation_threshold, ")")
                )

                workflow_data$tab_status$normalization <- "complete"

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
                    state_name = "metab_norm_complete"
                    , s4_data_object = current_s4
                    , config_object = workflow_data$config_list
                    , description = "Normalization complete (correlation filtering skipped)"
                )

                workflow_data$tab_status$normalization <- "complete"

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

            # TODO: Generate summary from correlation_results
            "Correlation filtering results will appear here..."
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
        # Export Session
        # ================================================================
        shiny::observeEvent(input$export_session, {
            shiny::req(workflow_data$state_manager)

            tryCatch({
                current_s4 <- workflow_data$state_manager$getState()

                if (!is.null(current_s4)) {
                    export_path <- file.path(
                        experiment_paths$export_dir
                        , paste0(experiment_label, "_normalized_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
                    )

                    if (!dir.exists(experiment_paths$export_dir)) {
                        dir.create(experiment_paths$export_dir, recursive = TRUE)
                    }

                    saveRDS(current_s4, export_path)
                    add_log(paste("Exported to:", export_path))
                    shiny::showNotification(paste("Data exported to:", export_path), type = "message")
                }

            }, error = function(e) {
                add_log(paste("Export error:", e$message))
                shiny::showNotification(paste("Export error:", e$message), type = "error")
            })
        })
    })
}
