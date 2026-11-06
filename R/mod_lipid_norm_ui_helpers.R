buildLipidNormOptionsControlPanel <- function(ns) {
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
            , testid(
                shiny::selectInput(
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
                , "lipid-norm-method"
            )
            , shiny::helpText("Method for between-sample normalization")

            # --- RUV-III Options ---
            , shiny::hr()
            , shiny::h4("RUV-III Batch Correction")
            , testid(
                shiny::radioButtons(
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
                , "lipid-norm-ruv-mode"
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
                        , "Weighted Difference (deprecated)" = "weighted_difference"
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
            , testid(
                shiny::actionButton(
                    ns("run_normalization")
                    , "Run Normalization Pipeline"
                    , class = "btn-primary"
                    , width = "100%"
                    , icon = shiny::icon("play")
                )
                , "lipid-norm-run"
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

            , testid(
                shiny::actionButton(
                    ns("export_session")
                    , "Export Normalized Data"
                    , class = "btn-success"
                    , width = "100%"
                    , icon = shiny::icon("download")
                )
                , "lipid-norm-export-session"
            )
            , shiny::helpText("Save normalized data for later analysis")
        )
    )
}

#' One figure tab's rows: per assay its label and the three stages' images (`<kind>_<stage>_assay<n>`).
#' @param labels The three stages' captions.
#' @param assays The assay rows (1 and 2 for lipidomics; KIT-097: 1 under the transcriptomics profile).
#' @noRd
buildLipidNormStageRows <- function(ns, kind, labels = c("Post-Filtering", "Post-Normalization", "RUV-Corrected"), assays = 1:2) {
    stages <- c("post_filter", "post_norm", "ruv_corrected")
    rows <- lapply(assays, function(n) shiny::tagList(
        shiny::fluidRow(shiny::column(12, shiny::h5(shiny::textOutput(ns(sprintf("assay%d_label_%s", n, kind))), style = "margin-bottom: 5px;"))),
        shiny::fluidRow(lapply(seq_along(stages), function(i) shiny::column(4,
            shiny::tags$small(labels[[i]], style = "display: block; text-align: center;"),
            shinyjqui::jqui_resizable(shiny::imageOutput(ns(sprintf("%s_%s_assay%d", kind, stages[[i]], n)), height = "300px")))))))
    if (length(rows) > 1L) rows <- c(rows[1L], unlist(lapply(rows[-1L], function(row) list(shiny::hr(), row)), recursive = FALSE))
    shiny::tagList(rows)
}

#' The normalisation tab's QC panels. KIT-097: under an omic `profile` only `profile$norm_tabs` are shown (by name:
#' `itsd`, `ruv`, `pca`, `density`, `rle`, `correlation`, `correlation_filter`, `log`), with one assay row and the profile's
#' stage captions; `NULL` (lipidomics) shows every panel as before.
#' @noRd
buildLipidNormQcTabsetPanel <- function(ns, profile = NULL) {
    assays <- if (is.null(profile)) 1:2 else 1L
    labels <- profile$norm_stage_labels %||% c("Post-Filtering", "Post-Normalization", "RUV-Corrected")
    figure <- function(title, icon, kind) shiny::tabPanel(title, icon = shiny::icon(icon), shiny::br(), buildLipidNormStageRows(ns, kind, labels, assays))
    centred <- function(title, icon, heading, help, ...) shiny::tabPanel(title, icon = shiny::icon(icon), shiny::br(),
        shiny::fluidRow(shiny::column(12, shiny::h5(heading, style = "text-align: center;"), shiny::helpText(help), shiny::hr(), ...)))
    tabs <- list(
        itsd = centred("ITSD Selection", "flask", "Internal Standard Selection per Assay",
            "Select internal standards for each assay. Pre-selected candidates are auto-detected from naming patterns (ITSD, IS, labeled compounds).",
            shiny::uiOutput(ns("itsd_selection_ui"))),
        ruv = centred("RUV QC", "chart-line", profile$norm_ruv_heading %||% "RUV Parameter Optimization Results",
            profile$norm_ruv_help %||% "Canonical correlation plots and optimization summary per assay. Each assay is optimized independently.",
            shiny::uiOutput(ns("ruv_qc_ui"))),
        pca = figure("PCA", "project-diagram", "pca"),
        density = figure("Density", "chart-area", "density"),
        rle = figure("RLE", "chart-bar", "rle"),
        correlation = figure("Correlation", "th", "correlation"),
        correlation_filter = buildLipidNormCorrelationFilterTab(ns),
        log = shiny::tabPanel("Log", icon = shiny::icon("list-alt"), shiny::br(), shiny::fluidRow(shiny::column(12,
            shiny::h5("Normalization Pipeline Log", style = "text-align: center;"), shiny::verbatimTextOutput(ns("norm_log"))))))
    keep <- profile$norm_tabs %||% names(tabs)
    shiny::column(9, do.call(shiny::tabsetPanel, c(list(id = ns("norm_qc_tabs")), unname(tabs[keep]))))
}

#' The correlation filtering panel (the lipid normalisation's final step).
#' @noRd
buildLipidNormCorrelationFilterTab <- function(ns) {
    shiny::tabPanel(
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
}
