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

# KIT-105 W4: the normalisation tab's panel library (kit-101 design 4.7), from the lipid tab's
# helper files at f8d4c5fc (mod_lipid_norm_ui_helpers.R, mod_lipid_norm_support_helpers.R,
# mod_lipid_norm_runtime_helpers.R and mod_lipid_norm_observer_helpers.R; deleted by W4): the
# option sections of `config$norm$panels`, the QC tabset with its stage grid, the correlation-filter
# tab, and the renderers and registrars of the tab's outputs. The shared frame (mod_omic_norm.R)
# assembles them. Ids, texts and markup are the lipid tab's, byte for byte; KIT-110 added trx's
# TMM/logCPM and RUV grid sections and the reset's optional test id, trx's markup byte for byte.
# Recorded divergence from conventions/R-shiny.md: the wellPanel layout is kept (existing-code
# scope rule).

# ------------------------------------------------------------------------------------------------
# The option sections (column 3).
# ------------------------------------------------------------------------------------------------

#' The design-driven select of the aesthetics and the RUV grouping (`ns(id)`, its four default choices).
#' @noRd
omicNormDesignSelect <- function(ns, id, label) {
    shiny::selectInput(ns(id), label, choices = c("group", "factor1", "factor2", "batch"), selected = "group", width = "100%")
}

#' The option sections `norm$panels` names, in its order: each a list of tags in the well panel
#' (mod_lipid_norm_ui_helpers.R:6-261 at f8d4c5fc; KIT-110: trx's `tmm_logcpm` and `ruv_grid`, mod_trx_norm.R:26-47 at
#' 38d9d25e). Only the named sections are built. `norm` is the omic's `norm` group.
#' @noRd
omicNormOptionSections <- function(ns, norm) {
    panels <- norm$panels
    section <- function(name) switch(name,
        aesthetics = list(shiny::hr(), shiny::h5("Plot Aesthetics"), omicNormDesignSelect(ns, "color_variable", "Color by:"),
            shiny::helpText("Variable to use for plot coloring"), omicNormDesignSelect(ns, "shape_variable", "Shape by:"),
            shiny::helpText("Variable to use for plot shapes")),
        itsd = list(shiny::hr(), shiny::h4("ITSD Normalization"),
            shiny::checkboxInput(ns("apply_itsd"), "Apply ITSD normalization", value = TRUE, width = "100%"),
            shiny::helpText("Normalize using internal standards (per assay)"),
            shiny::conditionalPanel(condition = "input.apply_itsd == true", ns = ns,
                shiny::selectInput(ns("itsd_aggregation"), "ITSD Aggregation:", choices = panels$itsd$aggregation,
                    selected = unname(panels$itsd$aggregation[[1L]]), width = "100%"),
                shiny::helpText("Method to combine multiple internal standards per sample. Select standards in ITSD Selection tab."))),
        log2 = list(shiny::hr(), shiny::h4("Log2 Transform"),
            shiny::numericInput(ns("log_offset"), "Offset for zeros:", value = panels$log2$offset, min = 0, step = 0.1, width = "100%"),
            shiny::helpText("Applied after ITSD normalization")),
        between_sample = list(shiny::hr(), shiny::h4("Between-Sample Normalization"),
            testid(shiny::selectInput(ns("norm_method"), "Normalization Method:", choices = panels$between_sample$methods,
                selected = panels$between_sample$default, width = "100%"), norm$testids[["method"]]),
            shiny::helpText("Method for between-sample normalization")),
        tmm_logcpm = list(shiny::hr(), shiny::h4("TMM and logCPM"),
            shiny::helpText("edgeR TMM normalisation factors of the filtered counts; log2 counts per million with a prior count of 2.")),
        ruv_iii = omicNormRuvSection(ns, panels$ruv_iii, norm$testids[["ruv_mode"]]),
        ruv_grid = omicNormRuvGridSection(ns, panels$ruv_grid, norm$testids[["ruv_mode"]]),
        actions = list(shiny::hr(),
            testid(shiny::actionButton(ns("run_normalization"), "Run Normalization Pipeline", class = "btn-primary", width = "100%",
                icon = shiny::icon("play")), norm$testids[["run"]]),
            shiny::br(), shiny::br(),
            omicNormResetButton(ns, if ("reset" %in% names(norm$testids)) norm$testids[["reset"]])),
        export = list(shiny::br(), shiny::br(),
            testid(shiny::actionButton(ns(panels$export$button), "Export Normalized Data", class = "btn-success", width = "100%",
                icon = shiny::icon("download")), norm$testids[["export"]]),
            shiny::helpText("Save normalized data for later analysis")))
    lapply(stats::setNames(names(panels), names(panels)), section)
}

#' The reset button, with its test id when the omic has one (KIT-110: trx's; lipid and metab have none).
#' @noRd
omicNormResetButton <- function(ns, testId = NULL) {
    button <- shiny::actionButton(ns("reset_normalization"), "Reset to Pre-Normalization", class = "btn-warning", width = "100%",
        icon = shiny::icon("undo"))
    if (is.null(testId)) button else testid(button, testId)
}

#' The RUV grid section (KIT-110, trx; mod_trx_norm.R:26-47 at 38d9d25e): RUVs' modes, and the automatic grid's and the
#' manual run's inputs with the grid's values. The input bounds and help texts are trx's, literal.
#' @noRd
omicNormRuvGridSection <- function(ns, grid, testId) {
    automatic <- shiny::conditionalPanel(condition = "input.ruv_mode == 'automatic'", ns = ns,
        shiny::h5("Automatic Optimization Parameters"),
        shiny::fluidRow(shiny::column(4, shiny::numericInput(ns("auto_percentage_min"), "Min %:", value = grid$percentage_min, min = 1, max = 50, width = "100%")),
            shiny::column(4, shiny::numericInput(ns("auto_percentage_max"), "Max %:", value = grid$percentage_max, min = 1, max = 50, width = "100%")),
            shiny::column(4, shiny::numericInput(ns("auto_percentage_step"), "Step %:", value = grid$percentage_step, min = 1, max = 25, width = "100%"))),
        shiny::numericInput(ns("k_penalty_weight"), "K Value Penalty Weight:", value = grid$k_penalty_weight, min = 0, max = 0.9, step = 0.01,
            width = "100%"),
        shiny::helpText("The notebook's grid: 5 to 30 % by 5, k capped at 3 below 15 samples and 4 otherwise, penalty 0.05; the",
            "percentage and k with the maximum composite score are chosen."))
    manual <- shiny::conditionalPanel(condition = "input.ruv_mode == 'manual'", ns = ns,
        shiny::sliderInput(ns("ruv_percentage"), "% Genes as Inferred Negative Controls:", min = 1, max = 50, value = grid$percentage, step = 1,
            width = "100%"),
        shiny::numericInput(ns("ruv_k"), "Number of Factors (k):", value = grid$k, min = 1, max = 4, width = "100%"),
        shiny::helpText("k must not exceed the cap (3 below 15 samples, 4 otherwise); a larger k is refused."))
    list(shiny::hr(), shiny::h4("RUVs (RUVSeq)"),
        testid(shiny::radioButtons(ns("ruv_mode"), "RUV Parameter Tuning:", choices = grid$modes, selected = grid$default,
            width = "100%"), testId),
        shiny::helpText("RUVs uses inferred negative controls (the least significant genes in a preliminary group GLM). Its factors",
            "W_1..W_k enter the DA model as covariates; the RUVs-corrected logCPM is for display only."),
        automatic, manual)
}

#' The RUV-III section: the modes, the grouping, and the automatic, manual and skip panels
#' (mod_lipid_norm_ui_helpers.R:88-222 at f8d4c5fc), with the omic's control noun.
#' @noRd
omicNormRuvSection <- function(ns, ruv, testId) {
    list(shiny::hr(), shiny::h4("RUV-III Batch Correction"),
        testid(shiny::radioButtons(ns("ruv_mode"), "RUV Parameter Tuning:", choices = ruv$modes, selected = ruv$default,
            width = "100%"), testId),
        omicNormDesignSelect(ns, "ruv_grouping_variable", "RUV Grouping Variable:"),
        shiny::helpText("Experimental factor to preserve during batch correction (shared across assays)"),
        shiny::conditionalPanel(condition = "input.ruv_mode == 'automatic'", ns = ns,
            shiny::h5("Automatic Optimization Parameters"),
            shiny::helpText("Auto-optimization runs independently for each assay (pos/neg)"),
            shiny::fluidRow(
                shiny::column(6, shiny::numericInput(ns("auto_percentage_min"), "Min %:", value = 1, min = 1, max = 50, width = "100%")),
                shiny::column(6, shiny::numericInput(ns("auto_percentage_max"), "Max %:", value = 20, min = 1, max = 50, width = "100%"))),
            shiny::helpText("Range of percentages to test (1-50%)"),
            shiny::selectInput(ns("separation_metric"), "Separation Quality Metric:", choices = list(
                "Maximum Difference (recommended)" = "max_difference", "Mean Difference" = "mean_difference",
                "Area Under Curve" = "auc", "Weighted Difference (deprecated)" = "weighted_difference"),
                selected = "max_difference", width = "100%"),
            shiny::sliderInput(ns("k_penalty_weight"), "K Value Penalty Weight:", min = 0.1, max = 0.9, value = 0.5, step = 0.1,
                width = "100%"),
            shiny::helpText("Controls penalty for high k values (0.1 = lenient, 0.9 = strict)"),
            shiny::checkboxInput(ns("adaptive_k_penalty"), "Adaptive K Penalty (sample size aware)", value = TRUE, width = "100%"),
            shiny::br()),
        shiny::conditionalPanel(condition = "input.ruv_mode == 'manual'", ns = ns,
            shiny::sliderInput(ns("ruv_percentage"), sprintf("%% %s as Negative Controls:", ruv$noun), min = 1, max = 20,
                value = ruv$percentage, step = 1, width = "100%"),
            shiny::helpText(sprintf("Percentage of most stable %s used as negative controls", tolower(ruv$noun))),
            shiny::numericInput(ns("ruv_k"), "Number of Factors (k):", value = ruv$k, min = 1, max = 10, width = "100%"),
            shiny::helpText("Number of unwanted variation factors to remove"),
            shiny::br()),
        shiny::conditionalPanel(condition = "input.ruv_mode == 'skip'", ns = ns,
            shiny::tags$div(style = "background-color: #fff3cd; border: 1px solid #ffc107; padding: 10px; border-radius: 4px;",
                shiny::h5("RUV-III Will Be Skipped", style = "color: #856404; margin-top: 0;"),
                shiny::helpText("Batch correction will not be applied. Consider including batch as a covariate in DE model if batch effects are present.",
                    style = "margin-bottom: 0;")),
            shiny::br()))
}

#' The options column: `norm$panels`' sections in their order.
#' @noRd
omicNormOptionsColumn <- function(ns, norm) {
    sections <- omicNormOptionSections(ns, norm)
    shiny::column(3, do.call(shiny::wellPanel, c(list(shiny::h4("Normalization Options")), do.call(c, unname(sections)))))
}

# ------------------------------------------------------------------------------------------------
# The QC tabset (column 9).
# ------------------------------------------------------------------------------------------------

#' One figure tab's rows: per assay its label and the three stages' images (`<kind>_<stage>_assay<n>`;
#' mod_lipid_norm_ui_helpers.R:270-279 at f8d4c5fc).
#' @param labels The three stages' captions.
#' @param assays The assay rows.
#' @noRd
omicNormStageRows <- function(ns, kind, labels, assays) {
    stages <- c("post_filter", "post_norm", "ruv_corrected")
    rows <- lapply(assays, function(n) shiny::tagList(
        shiny::fluidRow(shiny::column(12, shiny::h5(shiny::textOutput(ns(sprintf("assay%d_label_%s", n, kind))), style = "margin-bottom: 5px;"))),
        shiny::fluidRow(lapply(seq_along(stages), function(i) shiny::column(4,
            shiny::tags$small(labels[[i]], style = "display: block; text-align: center;"),
            shinyjqui::jqui_resizable(shiny::imageOutput(ns(sprintf("%s_%s_assay%d", kind, stages[[i]], n)), height = "300px")))))))
    if (length(rows) > 1L) rows <- c(rows[1L], unlist(lapply(rows[-1L], function(row) list(shiny::hr(), row)), recursive = FALSE))
    shiny::tagList(rows)
}

#' The QC tabset `norm_qc_tabs`: the tabs of `view$plot_tabs` in their order, the stage grid of `view$assay_slots` rows
#' with `view$stage_labels`, and the RUV tab's `view$ruv_heading` and `view$ruv_help` (mod_lipid_norm_ui_helpers.R:285-307
#' at f8d4c5fc). `view` is the omic's `norm` group.
#' @noRd
omicNormQcTabset <- function(ns, view) {
    assays <- seq_len(view$assay_slots)
    figure <- function(title, icon, kind) shiny::tabPanel(title, icon = shiny::icon(icon), shiny::br(),
        omicNormStageRows(ns, kind, view$stage_labels, assays))
    centred <- function(title, icon, heading, help, ...) shiny::tabPanel(title, icon = shiny::icon(icon), shiny::br(),
        shiny::fluidRow(shiny::column(12, shiny::h5(heading, style = "text-align: center;"), shiny::helpText(help), shiny::hr(), ...)))
    tabs <- list(
        itsd = centred("ITSD Selection", "flask", "Internal Standard Selection per Assay",
            "Select internal standards for each assay. Pre-selected candidates are auto-detected from naming patterns (ITSD, IS, labeled compounds).",
            shiny::uiOutput(ns("itsd_selection_ui"))),
        ruv = centred("RUV QC", "chart-line", view$ruv_heading, view$ruv_help, shiny::uiOutput(ns("ruv_qc_ui"))),
        pca = figure("PCA", "project-diagram", "pca"),
        density = figure("Density", "chart-area", "density"),
        rle = figure("RLE", "chart-bar", "rle"),
        correlation = figure("Correlation", "th", "correlation"),
        correlation_filter = omicNormCorrelationFilterTab(ns),
        log = shiny::tabPanel("Log", icon = shiny::icon("list-alt"), shiny::br(), shiny::fluidRow(shiny::column(12,
            shiny::h5("Normalization Pipeline Log", style = "text-align: center;"), shiny::verbatimTextOutput(ns("norm_log"))))))
    shiny::column(9, do.call(shiny::tabsetPanel, c(list(id = ns("norm_qc_tabs")), unname(tabs[view$plot_tabs]))))
}

#' The correlation filtering tab, the normalisation's final step (mod_lipid_norm_ui_helpers.R:311-371 at f8d4c5fc).
#' @noRd
omicNormCorrelationFilterTab <- function(ns) {
    shiny::tabPanel("Correlation Filtering", icon = shiny::icon("filter"), shiny::br(), shiny::fluidRow(
        shiny::column(4, shiny::wellPanel(
            shiny::h5("Sample Correlation Filtering"),
            shiny::helpText("FINAL STEP: Filter samples based on their correlation with other samples in the same group. Low correlation indicates potential technical issues."),
            shiny::sliderInput(ns("min_pearson_correlation_threshold"), "Min. Pearson Correlation Threshold:", min = 0.0, max = 1.0,
                value = 0.5, step = 0.05, width = "100%"),
            shiny::helpText("Samples with correlation below this threshold will be flagged"),
            shiny::br(),
            shiny::actionButton(ns("apply_correlation_filter"), "Apply Correlation Filter", class = "btn-primary", width = "100%",
                icon = shiny::icon("play")),
            shiny::br(), shiny::br(),
            shiny::actionButton(ns("skip_correlation_filter"), "Skip Filtering & Proceed to DE", class = "btn-info", width = "100%",
                icon = shiny::icon("forward")),
            shiny::helpText("Bypass sample filtering and proceed directly to DE"),
            shiny::br(), shiny::br(),
            shiny::h5("Filter Results"),
            shiny::verbatimTextOutput(ns("correlation_filter_summary")))),
        shiny::column(8,
            shiny::h5("Final Filtered Data QC", style = "text-align: center;"),
            shiny::helpText("QC plots for the final dataset after correlation filtering"),
            shinyjqui::jqui_resizable(shiny::plotOutput(ns("final_qc_plot"), height = "500px")))))
}

# ------------------------------------------------------------------------------------------------
# The tab's state, log and aesthetics.
# ------------------------------------------------------------------------------------------------

#' The tab's reactive values (mod_lipid_norm_runtime_helpers.R:1-23 at f8d4c5fc). KIT-075: the `*_obj` slots hold state
#' references, not objects.
#' @noRd
omicNormReactiveState <- function() {
    shiny::reactiveValues(assay_names = NULL, itsd_selections = list(), pre_norm_qc_generated = FALSE,
        normalization_complete = FALSE, ruv_complete = FALSE, correlation_filtering_complete = FALSE, qc_plot_paths = list(),
        post_filter_obj = NULL, post_itsd_obj = NULL, post_log2_obj = NULL, post_norm_obj = NULL, ruv_corrected_obj = NULL,
        correlation_filtered_obj = NULL, ruv_optimization_results = list(), correlation_results = list(),
        plot_refresh_trigger = 0, normalization_log = character(0))
}

#' The tab's log: `addLog(message)` appends `[HH:MM:SS] message` to `normData$normalization_log`.
#' @noRd
omicNormAddLog <- function(normData) {
    function(message) {
        normData$normalization_log <- c(normData$normalization_log, sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), message))
        invisible(normData$normalization_log)
    }
}

#' The plot aesthetics: the colour and shape variables, `group` when unset.
#' @noRd
omicNormPlotAestheticsGetter <- function(input) {
    function() {
        list(color_var = if (is.null(input$color_variable) || identical(input$color_variable, "")) "group" else input$color_variable,
            shape_var = if (is.null(input$shape_variable) || identical(input$shape_variable, "")) "group" else input$shape_variable)
    }
}

# ------------------------------------------------------------------------------------------------
# The stage grid's renderers and registrars (with an omic's file stem and stage prefixes: trx's, KIT-110).
# ------------------------------------------------------------------------------------------------

#' An image slot's renderer factory: `function(assaySlot, plotType, stagePrefix)` gives the renderImage of
#' `<sanitised assay>_<stagePrefix>_<plotType>.png` in `experimentPaths[[qcDirKey]]`.
#' @param sanitizeAssayNameFn The file stem of an assay name (trx: `trx`, its `norm$image_stem`).
#' @param qcDirKey The QC directory's key in `experimentPaths`.
#' @noRd
omicNormQcImageRenderer <- function(normData, experimentPaths,
    sanitizeAssayNameFn = function(assayName) gsub("[^A-Za-z0-9]", "_", tolower(assayName)), qcDirKey) {
    function(assaySlot, plotType, stagePrefix) {
        # Bind this output's own spec now: the registrars call this in a loop (F105-01).
        force(assaySlot)
        force(plotType)
        force(stagePrefix)
        shiny::renderImage({
            normData$plot_refresh_trigger
            assay_names <- normData$assay_names
            if (is.null(assay_names) || length(assay_names) < assaySlot) return(list(src = "", alt = "Assay not detected yet"))
            assay_name <- assay_names[[assaySlot]]
            filename <- paste0(sanitizeAssayNameFn(assay_name), "_", stagePrefix, "_", plotType, ".png")
            qc_dir <- experimentPaths[[qcDirKey]]
            if (is.null(qc_dir)) return(list(src = "", alt = "QC directory not configured"))
            img_path <- file.path(qc_dir, filename)
            if (!file.exists(img_path)) return(list(src = "", alt = paste("Plot not generated yet:", filename)))
            list(src = img_path, contentType = "image/png", width = "100%", height = "auto", alt = paste(plotType, "-", assay_name))
        }, deleteFile = FALSE)
    }
}

#' The stage grid's image slots `<kind>_<post_filter|post_norm|ruv_corrected>_assay<n>`, for `assaySlots` assays, each the
#' renderer of its assay, kind and stage file prefix.
#' @param stagePrefixes The file prefix of each stage (KIT-097: trx's `raw`, `tmm`, `ruvs`, its `norm$stage_prefixes`).
#' @noRd
registerOmicNormQcImageOutputs <- function(output, renderQcImageForAssay,
    stagePrefixes = c(pre_norm = "pre_norm", post_norm = "post_norm", ruv_corrected = "ruv_corrected"), assaySlots = 2L) {
    stages <- c(post_filter = "pre_norm", post_norm = "post_norm", ruv_corrected = "ruv_corrected")
    for (kind in c("pca", "density", "rle", "correlation")) {
        for (n in seq_len(assaySlots)) {
            for (stage in names(stages)) {
                output[[sprintf("%s_%s_assay%d", kind, stage, n)]] <- renderQcImageForAssay(n, kind, stagePrefixes[[stages[[stage]]]])
            }
        }
    }
    invisible(output)
}

#' An assay label's renderer factory: `function(assaySlot)` gives `Assay: <name>`, or `Assay <n>: (detecting...)`.
#' @noRd
omicNormAssayLabelRenderer <- function(normData) {
    function(assaySlot) {
        # Bind this output's own slot now: the registrars call this in a loop (F105-01).
        force(assaySlot)
        shiny::renderText({
            if (!is.null(normData$assay_names) && length(normData$assay_names) >= assaySlot) {
                paste("Assay:", normData$assay_names[[assaySlot]])
            } else {
                paste0("Assay ", assaySlot, ": (detecting...)")
            }
        })
    }
}

#' The stage grid's labels `assay<n>_label_<kind>`.
#' @noRd
registerOmicNormAssayLabelOutputs <- function(output, renderAssayLabel, assaySlots = 2L) {
    for (kind in c("pca", "density", "rle", "correlation")) {
        for (n in seq_len(assaySlots)) output[[sprintf("assay%d_label_%s", n, kind)]] <- renderAssayLabel(n)
    }
    invisible(output)
}

#' The log tab's renderer factory.
#' @noRd
omicNormLogRenderer <- function(normData) {
    function() {
        shiny::renderText({
            if (length(normData$normalization_log) == 0) return("Normalization log will appear here as you apply steps...")
            paste(normData$normalization_log, collapse = "\n")
        })
    }
}

#' @noRd
registerOmicNormLogOutput <- function(output, renderNormLog) {
    output$norm_log <- renderNormLog()
    invisible(output)
}

#' The RUV QC tab's renderer factory: per assay its canonical-correlation plot, summary and results table.
#' @noRd
omicNormRuvQcUiRenderer <- function(normData, ns) {
    function() {
        shiny::renderUI({
            shiny::req(normData$assay_names)
            shiny::tagList(purrr::map(normData$assay_names, function(assay_name) {
                safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
                shiny::tagList(
                    shiny::fluidRow(shiny::column(12, shiny::h5(paste("Assay:", assay_name),
                        style = "border-bottom: 1px solid #ddd; padding-bottom: 5px;"))),
                    shiny::fluidRow(
                        shiny::column(8, shinyjqui::jqui_resizable(shiny::plotOutput(ns(paste0("cancor_plot_", safe_name)), height = "400px"))),
                        shiny::column(4, shiny::wellPanel(shiny::h6("Optimization Summary"),
                            shiny::verbatimTextOutput(ns(paste0("ruv_summary_", safe_name))), shiny::br(),
                            shiny::h6("Results Table"), DT::dataTableOutput(ns(paste0("ruv_table_", safe_name)))))),
                    shiny::hr())
            }))
        })
    }
}

#' @noRd
registerOmicNormRuvQcOutput <- function(output, renderRuvQcUi) {
    output$ruv_qc_ui <- renderRuvQcUi()
    invisible(output)
}

#' Each assay's RUV outputs once the optimisation has results: the canonical-correlation plot, the summary and the table
#' (mod_lipid_norm_observer_helpers.R:87-148 at f8d4c5fc).
#' @noRd
registerOmicNormRuvCancorOutputs <- function(output, normData) {
    shiny::observe({
        shiny::req(normData$assay_names)
        shiny::req(length(normData$ruv_optimization_results) > 0)
        purrr::walk(normData$assay_names, function(assay_name) {
            safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
            output[[paste0("cancor_plot_", safe_name)]] <- shiny::renderPlot({
                result <- normData$ruv_optimization_results[[assay_name]]
                if (!is.null(result) && isTRUE(result$success) && !is.null(result$cancor_plot)) return(result$cancor_plot)
                omicNormTextPlot("No data", 6)
            })
            output[[paste0("ruv_summary_", safe_name)]] <- shiny::renderText(omicNormRuvSummary(normData$ruv_optimization_results[[assay_name]]))
            output[[paste0("ruv_table_", safe_name)]] <- DT::renderDataTable({
                result <- normData$ruv_optimization_results[[assay_name]]
                if (!is.null(result) && isTRUE(result$success) && !is.null(result$optimization_results)) {
                    DT::datatable(result$optimization_results, options = list(pageLength = 5, dom = "t"), rownames = FALSE)
                } else {
                    NULL
                }
            })
        })
    })
    invisible(output)
}

#' One assay's RUV optimisation summary text.
#' @noRd
omicNormRuvSummary <- function(result) {
    if (is.null(result)) return("Not yet computed")
    if (!isTRUE(result$success)) return(paste0("Failed: ", result$error))
    paste0("Best k: ", result$best_k, "\n", "Best %: ", result$best_percentage, "\n", "Separation: ",
        round(result$separation_score, 4), "\n", "Controls: ", sum(result$control_genes_index, na.rm = TRUE))
}

#' A text placeholder plot.
#' @noRd
omicNormTextPlot <- function(label, size) {
    ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = label, size = size) + ggplot2::theme_void()
}

# ------------------------------------------------------------------------------------------------
# The ITSD selection, the correlation summary and the final QC plot.
# ------------------------------------------------------------------------------------------------

#' The ITSD selection tab's UI: one table per assay.
#' @noRd
omicNormItsdSelectionUi <- function(normData, ns) {
    shiny::renderUI({
        shiny::req(normData$assay_names)
        shiny::tagList(purrr::map(normData$assay_names, function(assay_name) {
            safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
            shiny::wellPanel(shiny::h5(paste("Assay:", assay_name)), DT::dataTableOutput(ns(paste0("itsd_table_", safe_name))), shiny::br())
        }))
    })
}

#' The ITSD selection tables, built from the head's state a page at a time (KIT-082; nothing whole kept), and the shown
#' rows' feature ids (`itsd_candidate_ids`, KIT-078: a managed run resolves a selection without the matrix). Synchronous,
#' so that Run cannot race the selection (mod_lipid_norm_observer_helpers.R:150-230 at f8d4c5fc).
#' @param config The omic's GUI configuration (its binding and `norm$itsd_table`).
#' @noRd
omicNormItsdTables <- function(output, workflowData, normData, config) {
    shiny::req(normData$assay_names)
    shiny::req(workflowData$state_manager)
    binding <- config$binding$binding
    shell <- tryCatch(readWorkflowOmicsShell(binding, workflowData$state_manager), error = function(e) NULL)
    if (is.null(shell) || !inherits(shell, binding$class)) return(invisible(NULL))
    id_col <- methods::slot(shell, binding$id_slot)
    annotation_cols <- shell@annotation_id_column
    rm(shell)
    build <- config$norm$itsd_table
    stored_tables <- tryCatch(readWorkflowOmicsItsdTables(binding, workflowData$state_manager, normData$assay_names,
        function(page) build(page, id_col, annotation_cols)), error = function(e) {
        logger::log_warn(paste("Could not read the ITSD candidate tables:", conditionMessage(e)))
        NULL
    })
    if (is.null(stored_tables)) return(invisible(NULL))
    selection_tables <- lapply(stats::setNames(normData$assay_names, normData$assay_names), function(assay_name) stored_tables[[assay_name]])
    normData$itsd_candidate_ids <- lapply(selection_tables, function(table) table$feature_id)
    purrr::walk(normData$assay_names, function(assay_name) {
        selection_table <- selection_tables[[assay_name]]
        output[[paste0("itsd_table_", gsub("[^A-Za-z0-9]", "_", tolower(assay_name)))]] <- DT::renderDataTable({
            if (is.null(selection_table)) return(NULL)
            DT::datatable(selection_table, selection = list(mode = "multiple", selected = which(selection_table$is_candidate)),
                filter = "top", options = list(pageLength = 10, scrollX = TRUE, order = list(list(4, "desc"), list(3, "asc"))),
                rownames = FALSE) |>
                DT::formatStyle("is_candidate", backgroundColor = DT::styleEqual(TRUE, "#d4edda")) |>
                DT::formatRound(columns = c("mean_intensity", "cv_percent"), digits = 2)
        })
    })
    invisible(NULL)
}

#' The ITSD selection's tracking: each assay table's selected rows into `normData$itsd_selections`.
#' @noRd
registerOmicNormItsdSelectionTracking <- function(input, normData) {
    shiny::observe({
        shiny::req(normData$assay_names)
        purrr::walk(normData$assay_names, function(assay_name) {
            input_id <- paste0("itsd_table_", gsub("[^A-Za-z0-9]", "_", tolower(assay_name)), "_rows_selected")
            shiny::observeEvent(input[[input_id]], {
                normData$itsd_selections[[assay_name]] <- input[[input_id]]
                logger::log_info(paste("ITSD selection updated for", assay_name, ":", length(input[[input_id]]), "features selected"))
            }, ignoreNULL = FALSE)
        })
    })
    invisible(normData)
}

#' The correlation filter's summary: per assay the pairs and their correlations, and the samples kept (the designs from
#' the states' metadata, KIT-075).
#' @param readMetadata The omic's metadata reader (`config$binding$read_metadata`).
#' @noRd
omicNormCorrelationSummary <- function(normData, workflowData, readMetadata) {
    shiny::renderText({
        if (!normData$correlation_filtering_complete) return("Apply correlation filter to see results...")
        corr_results <- normData$correlation_results
        filtered_obj <- normData$correlation_filtered_obj
        original_obj <- if (!is.null(normData$ruv_corrected_obj)) normData$ruv_corrected_obj else normData$post_norm_obj
        if (is.null(corr_results) || length(corr_results) == 0) return("No correlation results available.")
        summary_lines <- c("=== Correlation Filtering Summary ===\n")
        for (assay_name in names(corr_results)) {
            assay_corr <- corr_results[[assay_name]]
            if (!is.null(assay_corr) && nrow(assay_corr) > 0) {
                values <- assay_corr$pearson_correlation
                summary_lines <- c(summary_lines, sprintf("\n[%s]\n  Sample pairs: %d\n  Correlation: mean=%.3f, min=%.3f, max=%.3f",
                    assay_name, nrow(assay_corr), round(mean(values, na.rm = TRUE), 3), round(min(values, na.rm = TRUE), 3),
                    round(max(values, na.rm = TRUE), 3)))
            }
        }
        if (!is.null(original_obj) && !is.null(filtered_obj)) {
            original_samples <- nrow(readMetadata(workflowData$state_manager, original_obj)$design_matrix)
            filtered_samples <- nrow(readMetadata(workflowData$state_manager, filtered_obj)$design_matrix)
            summary_lines <- c(summary_lines, sprintf(
                "\n\n[Sample Filtering]\n  Original: %d samples\n  After filtering: %d samples\n  Removed: %d samples",
                original_samples, filtered_samples, original_samples - filtered_samples))
        }
        paste(summary_lines, collapse = "")
    })
}

#' The final QC plot: the PCA of the correlation-filtered state (else the RUV-corrected, else the normalised). KIT-082: a
#' managed session's worker draws it (`workflowOmicsFinalPca`), with the omic's placeholder sizes (D7); in process the
#' state is read for this render only (KIT-075) and the binding's `final_pca` draws it.
#' @param getPlotAesthetics The tab's aesthetics getter.
#' @param config The omic's GUI configuration.
#' @noRd
omicNormFinalQcPlot <- function(normData, getPlotAesthetics, workflowData, config) {
    binding <- config$binding$binding
    sizes <- config$norm$final_pca
    shiny::renderPlot({
        shiny::req(normData$correlation_filtering_complete || normData$ruv_complete)
        current_reference <- normData$correlation_filtered_obj %||% normData$ruv_corrected_obj %||% normData$post_norm_obj
        if (is.null(current_reference)) return(omicNormTextPlot("No data", 6))
        aesthetics <- getPlotAesthetics()
        runner <- workflowOmicsManagedActions(binding$omic_type)
        if (!is.null(runner)) {
            entry <- workflowOmicsFinalPca(binding, runner, workflowData$state_manager, current_reference, aesthetics$color_var,
                aesthetics$shape_var)
            if (is.null(entry)) return(omicNormTextPlot("Preparing the final QC plot...", sizes$preparing_size))
            if (!is.null(entry$message)) return(omicNormTextPlot(paste("Error:", entry$message), sizes$error_size))
            return(grid::grid.draw(entry$plot))
        }
        current_s4 <- config$binding$read_state(workflowData$state_manager, current_reference)
        tryCatch(binding$final_pca(current_s4, aesthetics$color_var, aesthetics$shape_var),
            error = function(e) omicNormTextPlot(paste("Error:", e$message), 4))
    })
}
