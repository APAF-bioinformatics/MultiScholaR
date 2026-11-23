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

# KIT-104 W5: the four shared QC sub-tabs (kit-101 design 4.6): the intensity filter, the duplicate
# resolution, the internal standards and the finalise. Each renders the lipid sub-tab's markup
# (mod_lipid_qc_{intensity,duplicates,itsd,s4}.R and their helper files at d8f6dde4) with the
# omic's texts, notification ids and settings from `config$qc`, and keeps its inputs, outputs,
# observers and notifications in the lipid order. The computation is the core's
# (func_omics_qc_shared.R, `runWorkflowOmicsQcAction`); a sub-tab holds the numbers it shows, never
# a state object (KIT-075, KIT-082). Every action ignores a click made before its server started
# (`qcClickAfterStart`, KIT-101 F-17 and F-24). Each server's module frame holds only its input,
# output and session, so a test's variables are not masked under testServer. The file:line notes
# below cite the lipid files at d8f6dde4; KIT-104 W5 deleted the two helper files
# (mod_lipid_qc_intensity_server_helpers.R, mod_lipid_qc_duplicates_server_helpers.R).

# ------------------------------------------------------------------------------------------------
# Shared pieces.
# ------------------------------------------------------------------------------------------------

#' The title of the qc tab whose component is `component` (the tabPanel's title).
#' @noRd
omicQcTabTitle <- function(config, component) {
    for (group in config$qc$groups) {
        for (tab in group$tabs) if (identical(tab$component, component)) return(tab$title)
    }
    cli::cli_abort("the qc configuration has no {.val {component}} tab.", class = "omic_gui_config_invalid")
}

#' The QC progress figure's output: a grob or gtable drawn, a ggplot printed (the lipid sub-tabs' `filter_plot`).
#' @noRd
registerOmicQcPlotOutput <- function(output, plot) {
    output$filter_plot <- shiny::renderPlot({
        drawn <- plot()
        shiny::req(drawn)
        if (inherits(drawn, "grob") || inherits(drawn, "gtable")) {
            grid::grid.draw(drawn)
        } else if (inherits(drawn, "ggplot")) {
            print(drawn)
        }
    })
    invisible(output)
}

#' A figure drawn for the tab, or NULL with the reason logged, as the sub-tabs did.
#' @noRd
omicQcDrawn <- function(figure) {
    tryCatch(figure(), error = function(e) {
        logger::log_warn(paste("Could not generate QC plot:", e$message))
        NULL
    })
}

#' The filtering-progress figure's drawing for a QC step: `config$qc$filtering_figure` with the step's label.
#' @noRd
omicQcDraw <- function(config, omic_type, step_name) {
    function(arguments) {
        do.call(config$qc$filtering_figure, c(arguments, list(step_name = step_name, omics_type = omic_type,
            return_grid = TRUE, overwrite = TRUE)))
    }
}

#' The head state's metadata, refused unless it is of the omic's class.
#' @noRd
omicQcHeadMetadata <- function(workflow_data, config) {
    current <- config$binding$read_metadata(workflow_data$state_manager)
    shiny::req(current)
    class <- config$binding$binding$class
    if (!identical(current$class, class)) stop(sprintf("Current state is not a %s object", class))
    current
}

# ------------------------------------------------------------------------------------------------
# The intensity filter.
# ------------------------------------------------------------------------------------------------

#' The intensity filter sub-tab's UI (mod_lipid_qc_intensity.R:37-102).
#' @param id Module id (`intensity`).
#' @param config An omic GUI configuration.
#' @noRd
mod_omic_qc_intensity_ui <- function(id, config) {
    ns <- shiny::NS(id)
    texts <- config$qc$texts
    shiny::tabPanel(omicQcTabTitle(config, "omic_qc_intensity"), shiny::br(), shiny::fluidRow(
        shiny::column(4, shiny::wellPanel(
            shiny::h4(texts[["intensity_heading"]]),
            shiny::p(texts[["intensity_help"]]),
            shiny::hr(),
            shiny::sliderInput(ns("intensity_cutoff_percentile"), "Intensity Cutoff Percentile (%)", min = 1, max = 50,
                value = 10, step = 1),
            shiny::helpText("Intensities below this percentile are considered 'low' (default: 10%)"),
            shiny::sliderInput(ns("proportion_below_cutoff"), "Max Proportion of Samples Below Threshold", min = 0.1,
                max = 1.0, value = 0.5, step = 0.05),
            shiny::helpText(texts[["intensity_proportion_help"]]),
            shiny::hr(),
            shiny::div(
                shiny::actionButton(ns("apply_filter"), "Apply Filter", class = "btn-primary", width = "48%"),
                shiny::actionButton(ns("revert_filter"), "Revert", class = "btn-warning", width = "48%",
                    style = "margin-left: 4%")))),
        shiny::column(8,
            shiny::verbatimTextOutput(ns("filter_results")),
            shiny::br(),
            shiny::uiOutput(ns("assay_results_tabs")),
            shiny::br(),
            shinyjqui::jqui_resizable(shiny::plotOutput(ns("filter_plot"), height = "600px", width = "100%")))))
}

#' The intensity filter sub-tab's server (mod_lipid_qc_intensity.R:113-148).
#' @param id Module id.
#' @param workflow_data The workflow bus.
#' @param omic_type,experiment_label The QC tab's omic and label.
#' @param config An omic GUI configuration.
#' @noRd
mod_omic_qc_intensity_server <- function(id, workflow_data, omic_type, experiment_label, config) {
    shiny::moduleServer(id, function(input, output, session) {
        registerOmicQcIntensity(input, output, session, workflow_data, omic_type, config)
    })
}

#' The intensity filter's per-assay result tabs (mod_lipid_qc_intensity_server_helpers.R:1-46).
#' @noRd
omicQcIntensityAssayTabsUi <- function(stats, ns, texts) {
    if (is.null(stats) || nrow(stats) == 0) return(NULL)
    row <- function(label, value) shiny::tags$tr(shiny::tags$td(shiny::strong(label)), shiny::tags$td(value))
    tabs <- lapply(seq_len(nrow(stats)), function(i) {
        shiny::tabPanel(stats$Assay[i], shiny::br(), shiny::tags$table(class = "table table-striped table-condensed",
            style = "width: auto;", shiny::tags$tbody(
                row(texts[["intensity_original_label"]], stats$Original[i]),
                row("After filtering:", stats$Filtered[i]),
                row("Removed:", stats$Removed[i]),
                row("Percent retained:", paste0(stats$Percent_Retained[i], "%")))))
    })
    do.call(shiny::tabsetPanel, c(list(id = ns("assay_stats_tabs")), tabs))
}

#' The intensity filter's results text (mod_lipid_qc_intensity_server_helpers.R:197-223).
#' @noRd
omicQcIntensityText <- function(stats, percentile, proportion, texts, prefix) {
    paste(texts[["intensity_result_heading"]], "============================================",
        sprintf("Intensity cutoff percentile: %d%%", percentile),
        sprintf("Max proportion below cutoff: %.2f", proportion),
        "", "Per-Assay Results:",
        paste(vapply(seq_len(nrow(stats)), function(i) {
            sprintf("  %s: %d -> %d (removed %d, %.1f%% retained)", stats$Assay[i], stats$Original[i], stats$Filtered[i],
                stats$Removed[i], stats$Percent_Retained[i])
        }, character(1)), collapse = "\n"),
        "", sprintf(texts[["intensity_total"]], sum(stats$Original), sum(stats$Filtered), sum(stats$Removed)),
        "", sprintf("State saved as: '%s_intensity_filtered'", prefix), sep = "\n")
}

#' The intensity filter sub-tab inside its module session: Apply (the declared filter through the core's QC action,
#' its counts, text and figure) and Revert, then the outputs (mod_lipid_qc_intensity_server_helpers.R:126-271).
#' @noRd
registerOmicQcIntensity <- function(input, output, session, workflow_data, omic_type, config) {
    binding <- config$binding$binding
    texts <- config$qc$texts
    working <- config$qc$notification_ids[["intensity"]]
    stats <- shiny::reactiveVal(NULL)
    plot <- shiny::reactiveVal(NULL)

    apply_after_start <- qcClickAfterStart(input, "apply_filter")
    shiny::observeEvent(input$apply_filter, {
        if (!apply_after_start()) return(invisible(NULL))
        shiny::req(workflow_data$state_manager)
        shiny::showNotification(texts[["intensity_working"]], id = working, duration = NULL)
        tryCatch({
            # KIT-075: the counts come from the declaration's result and the figure from the committed state read for
            # this action (or the worker's metrics); no object outlives the observer.
            omicQcHeadMetadata(workflow_data, config)
            percentile <- input$intensity_cutoff_percentile
            proportion <- input$proportion_below_cutoff
            logger::log_info(paste("Applying", binding$feature_noun, "intensity filter: percentile =", percentile,
                ", proportion =", proportion))
            deliver <- function(result, figure) {
                table <- omicsQcCountsTable(result$value$counts)
                stats(table)
                plot(omicQcDrawn(figure))
                output$filter_results <- shiny::renderText(omicQcIntensityText(table, percentile, proportion, texts,
                    binding$prefix))
                noun <- binding$feature_noun
                logger::log_info(paste(paste0(toupper(substr(noun, 1L, 1L)), substring(noun, 2L)),
                    "intensity filter applied: removed", sum(table$Removed), paste0(noun, "s")))
                shiny::removeNotification(working)
                shiny::showNotification(sprintf(texts[["intensity_success"]], sum(table$Filtered)), type = "message")
            }
            runWorkflowOmicsQcAction(binding, workflow_data, omicsIntensityFilterStep(binding, percentile, proportion),
                deliver, omicQcDraw(config, omic_type, "2_Intensity_Filtered"),
                failed = function(error) shiny::removeNotification(working))
        }, error = function(e) {
            message <- paste(texts[["intensity_error"]], e$message)
            logger::log_error(message)
            shiny::showNotification(message, type = "error", duration = 15)
            shiny::removeNotification(working)
        })
    })

    registerOmicQcRevertObserver(input, output, "revert_filter", workflow_data, "filter_results", function() {
        stats(NULL)
        plot(NULL)
    }, paste("Reverted", binding$feature_noun, "intensity filter to"))

    output$assay_results_tabs <- shiny::renderUI(omicQcIntensityAssayTabsUi(stats(), session$ns, texts))
    registerOmicQcPlotOutput(output, plot)
    invisible(NULL)
}

# ------------------------------------------------------------------------------------------------
# The duplicate resolution.
# ------------------------------------------------------------------------------------------------

#' The duplicate resolution sub-tab's UI (mod_lipid_qc_duplicates.R:38-95).
#' @param id Module id (`duplicates`).
#' @param config An omic GUI configuration.
#' @noRd
mod_omic_qc_duplicates_ui <- function(id, config) {
    ns <- shiny::NS(id)
    shiny::tabPanel(omicQcTabTitle(config, "omic_qc_duplicates"), shiny::br(), shiny::fluidRow(
        shiny::column(4, shiny::wellPanel(
            shiny::h4("Duplicate Feature Resolution"),
            shiny::p(config$qc$texts[["duplicates_help"]]),
            shiny::hr(),
            shiny::h5("Detection Summary"),
            shiny::uiOutput(ns("duplicate_summary")),
            shiny::hr(),
            shiny::div(shiny::actionButton(ns("detect_duplicates"), "Detect Duplicates", class = "btn-info",
                width = "100%")),
            shiny::br(),
            shiny::div(
                shiny::actionButton(ns("resolve_duplicates"), "Resolve Duplicates", class = "btn-primary", width = "48%"),
                shiny::actionButton(ns("revert_duplicates"), "Revert", class = "btn-warning", width = "48%",
                    style = "margin-left: 4%")))),
        shiny::column(8,
            shiny::verbatimTextOutput(ns("resolution_results")),
            shiny::br(),
            shiny::uiOutput(ns("duplicate_tables")),
            shiny::br(),
            shinyjqui::jqui_resizable(shiny::plotOutput(ns("filter_plot"), height = "600px", width = "100%")))))
}

#' The duplicate resolution sub-tab's server (mod_lipid_qc_duplicates.R:124-134).
#' @inheritParams mod_omic_qc_intensity_server
#' @noRd
mod_omic_qc_duplicates_server <- function(id, workflow_data, omic_type, experiment_label, config) {
    shiny::moduleServer(id, function(input, output, session) {
        registerOmicQcDuplicates(input, output, session, workflow_data, omic_type, config)
    })
}

#' The output id of an assay's duplicate table.
#' @noRd
omicQcDuplicateTableId <- function(assay) paste0("dup_table_", gsub("[^a-zA-Z0-9]", "_", assay))

#' The detection summary (mod_lipid_qc_duplicates_server_helpers.R:110-133).
#' @noRd
omicQcDuplicateSummaryUi <- function(duplicates) {
    if (is.null(duplicates)) {
        return(shiny::p(shiny::icon("info-circle"), " Click 'Detect Duplicates' to scan for duplicate features.",
            style = "color: #666;"))
    }
    items <- lapply(names(duplicates), function(assay) {
        count <- if (is.null(duplicates[[assay]])) 0 else nrow(duplicates[[assay]])
        shiny::tags$li(shiny::icon(if (count == 0) "check-circle" else "exclamation-triangle",
            style = paste("color:", if (count == 0) "green" else "orange")), sprintf(" %s: %d duplicate IDs", assay, count))
    })
    shiny::tags$ul(items, style = "list-style: none; padding-left: 0;")
}

#' The per-assay duplicate tables' tabs (mod_lipid_qc_duplicates_server_helpers.R:11-39).
#' @noRd
omicQcDuplicateTablesUi <- function(duplicates, ns, texts) {
    if (is.null(duplicates)) return(NULL)
    found <- names(duplicates)[vapply(duplicates, function(table) !is.null(table) && nrow(table) > 0, logical(1))]
    if (length(found) == 0) {
        return(shiny::wellPanel(shiny::icon("check-circle", style = "color: green; font-size: 24px;"),
            shiny::h5("No duplicates found in any assay!"), shiny::p(texts[["duplicates_none"]])))
    }
    tabs <- lapply(found, function(assay) {
        shiny::tabPanel(assay, shiny::br(), DT::DTOutput(ns(omicQcDuplicateTableId(assay))))
    })
    do.call(shiny::tabsetPanel, c(list(id = ns("dup_tables_tabs")), tabs))
}

#' The duplicate resolution's results text (mod_lipid_qc_duplicates_server_helpers.R:278-300).
#' @noRd
omicQcDuplicatesText <- function(resolved, prefix) {
    removed <- sum(vapply(resolved, function(stat) stat$removed, numeric(1)))
    paste("Duplicate Resolution Complete", "==============================",
        "Strategy: Keep feature with highest average intensity", "", "Per-Assay Results:",
        paste(vapply(names(resolved), function(assay) {
            stat <- resolved[[assay]]
            sprintf("  %s: %d -> %d rows (removed %d duplicates)", assay, stat$original, stat$resolved, stat$removed)
        }, character(1)), collapse = "\n"),
        "", sprintf("Total duplicate rows removed: %d", removed),
        "", sprintf("State saved as: '%s_duplicates_resolved'", prefix), sep = "\n")
}

#' The duplicate resolution sub-tab inside its module session: Detect (the core's detection on the id column), the
#' summary and tables, Resolve (the declared resolution through the core's QC action), Revert and the figure
#' (mod_lipid_qc_duplicates_server_helpers.R:552-642).
#' @noRd
registerOmicQcDuplicates <- function(input, output, session, workflow_data, omic_type, config) {
    binding <- config$binding$binding
    texts <- config$qc$texts
    working <- config$qc$notification_ids[["duplicates"]]
    duplicates <- shiny::reactiveVal(NULL)
    resolved <- shiny::reactiveVal(NULL)
    plot <- shiny::reactiveVal(NULL)

    # KIT-101 F-24: Detect is display-only, and a click made before the start is not run either.
    detect_after_start <- qcClickAfterStart(input, "detect_duplicates")
    shiny::observeEvent(input$detect_duplicates, {
        if (!detect_after_start()) return(invisible(NULL))
        shiny::req(workflow_data$state_manager)
        tryCatch({
            # KIT-082: the finder reads the id column and nothing else.
            found <- detectWorkflowOmicsDuplicates(binding, workflow_data$state_manager)
            logger::log_info(paste("Detected", found$total, "duplicate feature IDs across assays"))
            duplicates(found$duplicates)
            shiny::showNotification(sprintf("Detection complete: %d duplicate IDs found", found$total),
                type = if (found$total > 0) "warning" else "message")
        }, error = function(e) {
            message <- paste("Error detecting duplicates:", e$message)
            logger::log_error(message)
            shiny::showNotification(message, type = "error")
        })
    })

    output$duplicate_summary <- shiny::renderUI(omicQcDuplicateSummaryUi(duplicates()))
    output$duplicate_tables <- shiny::renderUI(omicQcDuplicateTablesUi(duplicates(), session$ns, texts))
    shiny::observe({
        found <- duplicates()
        shiny::req(found)
        lapply(names(found), function(assay) {
            table <- found[[assay]]
            if (is.null(table) || nrow(table) == 0) return(NULL)
            output[[omicQcDuplicateTableId(assay)]] <- DT::renderDT(DT::datatable(table,
                options = list(pageLength = 10, scrollX = TRUE, dom = "frtip"), rownames = FALSE,
                class = "compact stripe"))
            NULL
        })
    })

    resolve_after_start <- qcClickAfterStart(input, "resolve_duplicates") # KIT-101 F-24
    shiny::observeEvent(input$resolve_duplicates, {
        if (!resolve_after_start()) return(invisible(NULL))
        shiny::req(workflow_data$state_manager)
        shiny::showNotification("Resolving duplicate features...", id = working, duration = NULL)
        tryCatch({
            # KIT-075: the assay names from the state's metadata; the result carries the state's reference.
            assays <- omicQcHeadMetadata(workflow_data, config)$assay_names
            deliver <- function(result, figure) {
                stats <- omicsDuplicateResolutionStats(result$value, assays)
                resolved(stats)
                plot(omicQcDrawn(figure))
                output$resolution_results <- shiny::renderText(omicQcDuplicatesText(stats, binding$prefix))
                duplicates(NULL)
                removed <- sum(vapply(stats, function(stat) stat$removed, numeric(1)))
                logger::log_info(paste("Resolved duplicates: removed", removed, "rows"))
                shiny::removeNotification(working)
                shiny::showNotification(sprintf("Duplicates resolved: %d rows removed", removed), type = "message")
            }
            runWorkflowOmicsQcAction(binding, workflow_data, omicsDuplicateResolutionStep(binding), deliver,
                omicQcDraw(config, omic_type, "3_Duplicates_Resolved"),
                failed = function(error) shiny::removeNotification(working))
        }, error = function(e) {
            message <- paste("Error resolving duplicates:", e$message)
            logger::log_error(message)
            shiny::showNotification(message, type = "error", duration = 15)
            shiny::removeNotification(working)
        })
    })

    registerOmicQcRevertObserver(input, output, "revert_duplicates", workflow_data, "resolution_results", function() {
        resolved(NULL)
        duplicates(NULL)
        plot(NULL)
    }, "Reverted duplicate resolution to")

    registerOmicQcPlotOutput(output, plot)
    invisible(NULL)
}

# ------------------------------------------------------------------------------------------------
# The internal standards.
# ------------------------------------------------------------------------------------------------

#' The internal standards sub-tab's UI (mod_lipid_qc_itsd.R:37-80).
#' @param id Module id (`itsd`).
#' @param config An omic GUI configuration.
#' @noRd
mod_omic_qc_itsd_ui <- function(id, config) {
    ns <- shiny::NS(id)
    shiny::tabPanel(omicQcTabTitle(config, "omic_qc_itsd"), shiny::br(), shiny::fluidRow(
        shiny::column(4, shiny::wellPanel(
            shiny::h4("Internal Standard QC"),
            shiny::p(paste("Visualize internal standard (IS) performance across samples. Good IS should have low CV",
                "and consistent intensities.")),
            shiny::hr(),
            shiny::textInput(ns("is_pattern"), "Internal Standard Regex Pattern", value = "",
                placeholder = "e.g., ^IS_|_d[0-9]+$|ISTD"),
            shiny::helpText(config$qc$texts[["itsd_pattern_help"]]),
            shiny::hr(),
            shiny::actionButton(ns("analyze_is"), "Analyze Internal Standards", class = "btn-primary", width = "100%"),
            shiny::hr(),
            shiny::h5("IS Detection Summary"),
            shiny::uiOutput(ns("is_summary")))),
        shiny::column(8,
            shiny::verbatimTextOutput(ns("is_results")),
            shiny::br(),
            shiny::uiOutput(ns("is_viz_tabs")))))
}

#' The internal standards sub-tab's server (mod_lipid_qc_itsd.R:210-430).
#' @inheritParams mod_omic_qc_intensity_server
#' @noRd
mod_omic_qc_itsd_server <- function(id, workflow_data, omic_type, experiment_label, config) {
    shiny::moduleServer(id, function(input, output, session) {
        registerOmicQcItsd(input, output, session, workflow_data, config)
    })
}

#' The detection summary per assay with its median CV (mod_lipid_qc_itsd.R:82-122).
#' @noRd
omicQcItsdSummaryUi <- function(metrics) {
    if (is.null(metrics)) {
        return(shiny::p(shiny::icon("info-circle"), " Click 'Analyze' to detect internal standards.",
            style = "color: #666;"))
    }
    by_assay <- split(metrics, metrics$assay)
    items <- lapply(names(by_assay), function(assay) {
        median_cv <- stats::median(by_assay[[assay]]$cv, na.rm = TRUE)
        color <- if (median_cv <= 15) "green" else if (median_cv <= 30) "orange" else "red"
        icon <- if (median_cv <= 15) "check-circle" else if (median_cv <= 30) "exclamation-circle" else "times-circle"
        shiny::tags$li(shiny::icon(icon, style = paste("color:", color)),
            sprintf(" %s: %d IS (median CV: %.1f%%)", assay, nrow(by_assay[[assay]]), median_cv))
    })
    shiny::tagList(shiny::tags$ul(items, style = "list-style: none; padding-left: 0;"), shiny::hr(),
        shiny::tags$small(shiny::icon("info-circle"), " CV < 15%: Good | 15-30%: Acceptable | > 30%: Review",
            style = "color: #666;"))
}

#' The CV and intensity figures' tabs (mod_lipid_qc_itsd.R:124-153).
#' @noRd
omicQcItsdTabsUi <- function(metrics, ns) {
    if (is.null(metrics)) return(NULL)
    tabs <- list(
        shiny::tabPanel("CV Distribution", shiny::br(),
            shinyjqui::jqui_resizable(shiny::plotOutput(ns("cv_plot"), height = "500px"))),
        shiny::tabPanel("Intensity Trends", shiny::br(),
            shinyjqui::jqui_resizable(shiny::plotOutput(ns("intensity_plot"), height = "500px"))))
    do.call(shiny::tabsetPanel, c(list(id = ns("is_viz_tabset")), tabs))
}

#' The internal standards sub-tab inside its module session: Analyze (the core's analysis on the rows it can use),
#' the summary, the figures' tabs and the two figures.
#' @noRd
registerOmicQcItsd <- function(input, output, session, workflow_data, config) {
    binding <- config$binding$binding
    working <- config$qc$notification_ids[["itsd"]]
    metrics <- shiny::reactiveVal(NULL)
    long <- shiny::reactiveVal(NULL)

    # KIT-101 F-24: Analyze is display-only, and a click made before the start is not run either.
    analyze_after_start <- qcClickAfterStart(input, "analyze_is")
    shiny::observeEvent(input$analyze_is, {
        if (!analyze_after_start()) return(invisible(NULL))
        shiny::req(workflow_data$state_manager)
        shiny::showNotification("Analyzing internal standards...", id = working, duration = NULL)
        tryCatch({
            # KIT-082: only the rows the kernel can use are read; the tab keeps the metrics and the long table.
            found <- analyzeWorkflowOmicsItsd(binding, workflow_data$state_manager, input$is_pattern, config$qc$itsd)
            metrics(found$metrics)
            # With no intensities to read the figure is cleared (metab's handling; lipid's kept the previous one).
            long(found$long_data)
            output$is_results <- shiny::renderText(found$result_text)
            logger::log_info(paste("IS analysis complete:", found$total, "standards found"))
            shiny::removeNotification(working)
            shiny::showNotification(sprintf("Found %d internal standards", found$total), type = "message")
        }, error = function(e) {
            message <- paste("Error analyzing internal standards:", e$message)
            logger::log_error(message)
            shiny::showNotification(message, type = "error", duration = 10)
            shiny::removeNotification(working)
        })
    })

    output$is_summary <- shiny::renderUI(omicQcItsdSummaryUi(metrics()))
    output$is_viz_tabs <- shiny::renderUI(omicQcItsdTabsUi(metrics(), session$ns))
    output$cv_plot <- shiny::renderPlot({
        shown <- metrics()
        shiny::req(shown)
        omicsItsdCvPlot(shown)
    })
    output$intensity_plot <- shiny::renderPlot({
        shown <- long()
        shiny::req(shown)
        omicsItsdIntensityPlot(shown)
    })
    invisible(NULL)
}

# ------------------------------------------------------------------------------------------------
# The finalise.
# ------------------------------------------------------------------------------------------------

#' The finalise sub-tab's UI (mod_lipid_qc_s4.R:39-93): the per-assay statistics only with `finalize_per_assay`.
#' @param id Module id (`s4_finalize`).
#' @param config An omic GUI configuration.
#' @noRd
mod_omic_qc_finalize_ui <- function(id, config) {
    ns <- shiny::NS(id)
    shiny::tabPanel(omicQcTabTitle(config, "omic_qc_finalize"), shiny::br(),
        shiny::fluidRow(
            shiny::column(6, shiny::wellPanel(
                shiny::h4("QC Processing Summary"),
                shiny::p("Review the QC processing steps before finalizing. You can revert to any previous state if needed."),
                shiny::hr(),
                shiny::h5("Processing History"),
                shiny::uiOutput(ns("state_history")),
                shiny::hr(),
                shiny::div(shiny::actionButton(ns("finalize_qc"), "Finalize QC", class = "btn-success", width = "100%",
                    icon = shiny::icon("check"))))),
            shiny::column(6,
                shiny::wellPanel(shiny::h4("Current Data Summary"), shiny::hr(), shiny::uiOutput(ns("data_summary"))),
                if (isTRUE(config$qc$finalize_per_assay)) shiny::wellPanel(shiny::h4("Per-Assay Statistics"),
                    shiny::hr(), DT::DTOutput(ns("assay_stats_table"))))),
        shiny::fluidRow(shiny::column(12,
            shiny::verbatimTextOutput(ns("finalize_results")),
            shiny::br(),
            shinyjqui::jqui_resizable(shiny::plotOutput(ns("filter_plot"), height = "600px", width = "100%")))))
}

#' The finalise sub-tab's server (mod_lipid_qc_s4.R:597-656).
#' @inheritParams mod_omic_qc_intensity_server
#' @noRd
mod_omic_qc_finalize_server <- function(id, workflow_data, omic_type, experiment_label, config) {
    shiny::moduleServer(id, function(input, output, session) {
        registerOmicQcFinalize(input, output, workflow_data, omic_type, config)
    })
}

#' The processing history, the head marked (mod_lipid_qc_s4.R:95-123).
#' @noRd
omicQcStateHistoryUi <- function(history) {
    if (length(history) == 0) {
        return(shiny::p(shiny::icon("info-circle"), " No processing history available.", style = "color: #666;"))
    }
    items <- lapply(seq_along(history), function(i) {
        current <- i == length(history)
        shiny::tags$li(shiny::icon(if (current) "arrow-right" else "check",
            style = paste("color:", if (current) "blue" else "green")),
            shiny::tags$span(sprintf(" %d. %s", i, history[i]), style = if (current) "font-weight: bold;" else ""),
            if (current) shiny::tags$span(" (current)", style = "color: blue;") else NULL)
    })
    shiny::tags$ol(items, style = "list-style: none; padding-left: 0;")
}

#' The processing history's output (mod_lipid_qc_s4.R:125-145).
#' @noRd
registerOmicQcStateHistory <- function(output, workflow_data) {
    output$state_history <- shiny::renderUI({
        shiny::req(workflow_data$state_manager)
        history <- tryCatch(workflow_data$state_manager$getHistory(), error = function(e) character(0))
        omicQcStateHistoryUi(history)
    })
    invisible(output)
}

#' The features retained: unique ids per assay, counted from the store on a shell read (KIT-082).
#' @noRd
omicQcFeatureTotal <- function(shell, binding) {
    assays <- methods::slot(shell, binding$data_slot)
    id_column <- methods::slot(shell, binding$id_slot)
    sum(attr(shell, "feature_counts") %||% vapply(assays, function(assay) {
        if (id_column %in% names(assay)) length(unique(assay[[id_column]])) else nrow(assay)
    }, numeric(1)))
}

#' The data summary from the head's shell, with the per-assay counts (mod_lipid_qc_s4.R:147-213).
#' @noRd
omicQcDataSummaryUi <- function(shell, config) {
    binding <- config$binding$binding
    texts <- config$qc$texts
    if (is.null(shell) || !inherits(shell, binding$class)) {
        return(shiny::p(shiny::icon("exclamation-triangle", style = "color: orange;"),
            sprintf(" No %s object available.", binding$class)))
    }
    assays <- methods::slot(shell, binding$data_slot)
    samples <- if (length(assays) > 0) names(assays[[1]])[vapply(assays[[1]], is.numeric, logical(1))] else character(0)
    design <- shell@design_matrix
    groups <- if (nrow(design) > 0 && shell@group_id %in% names(design)) length(unique(design[[shell@group_id]])) else NA
    row <- function(label, value) shiny::tags$tr(shiny::tags$td(shiny::strong(label)), shiny::tags$td(value))
    shiny::tagList(shiny::tags$table(class = "table table-condensed", style = "margin-bottom: 0;", shiny::tags$tbody(
        row("Number of Assays:", length(assays)),
        row(texts[["finalize_total_label"]], format(omicQcFeatureTotal(shell, binding), big.mark = ",")),
        row("Number of Samples:", length(samples)),
        row("Experimental Groups:", if (!is.na(groups)) groups else "N/A"),
        row(texts[["finalize_id_label"]], methods::slot(shell, binding$id_slot)),
        row("Sample ID Column:", shell@sample_id))))
}

#' The one-table data summary from the head's metadata: the stored row count and the design; no table is read
#' (trx's `registerTrxQcSummaryOutput`, mod_trx_qc_filter.R:192-207 at d8f6dde4).
#' @noRd
omicQcMetadataSummaryUi <- function(meta, config) {
    binding <- config$binding$binding
    texts <- config$qc$texts
    if (is.null(meta) || !identical(meta$class, binding$class)) {
        return(shiny::p(shiny::icon("exclamation-triangle", style = "color: orange;"),
            sprintf(" No %s object available.", binding$class)))
    }
    row <- function(label, value) shiny::tags$tr(shiny::tags$td(shiny::strong(label)), shiny::tags$td(value))
    design <- meta$design_matrix
    shiny::tags$table(class = "table table-condensed", style = "margin-bottom: 0;", shiny::tags$tbody(
        row(texts[["finalize_total_label"]], format(sum(meta$assay_rows), big.mark = ",")),
        row("Number of Samples:", nrow(design)),
        row("Experimental Groups:", length(unique(design[[meta$slots$group_id]]))),
        row(texts[["finalize_id_label"]], meta$slots[[binding$id_slot]]),
        row("Sample ID Column:", meta$slots$sample_id)))
}

#' The per-assay statistics: features, samples and missingness, counted from the store (mod_lipid_qc_s4.R:238-286).
#' @return A data frame, or NULL with no state of the omic's class.
#' @noRd
omicQcAssayStatsFrame <- function(shell, config) {
    binding <- config$binding$binding
    if (is.null(shell) || !inherits(shell, binding$class)) return(NULL)
    assays <- methods::slot(shell, binding$data_slot)
    id_column <- methods::slot(shell, binding$id_slot)
    table <- data.frame(Assay = names(assays), Features = attr(shell, "feature_counts") %||% sapply(assays, function(a) {
        if (id_column %in% names(a)) length(unique(a[[id_column]])) else nrow(a)
    }), Samples = sapply(assays, function(a) sum(vapply(a, is.numeric, logical(1)))), stringsAsFactors = FALSE)
    names(table)[[2L]] <- config$qc$texts[["finalize_count_column"]]
    table$Missingness <- attr(shell, "missing_percent") %||% sapply(assays, function(a) {
        columns <- names(a)[vapply(a, is.numeric, logical(1))]
        if (length(columns) == 0) return(NA)
        missing <- sum(vapply(columns, function(column) sum(is.na(a[[column]]) | a[[column]] == 0), numeric(1)))
        round((missing / (nrow(a) * length(columns))) * 100, 1)
    })
    table
}

#' The finalise's data summary outputs: from the shell with the per-assay table (`finalize_per_assay`), or from the
#' metadata alone (mod_lipid_qc_s4.R:215-309; trx's mod_trx_qc_filter.R:192-207).
#' @param output The finalise sub-tab's output.
#' @param workflow_data The workflow bus.
#' @param config An omic GUI configuration.
#' @noRd
registerOmicQcFinalizeSummary <- function(output, workflow_data, config) {
    binding <- config$binding$binding
    if (!isTRUE(config$qc$finalize_per_assay)) {
        output$data_summary <- shiny::renderUI({
            shiny::req(workflow_data$state_manager)
            meta <- tryCatch(config$binding$read_metadata(workflow_data$state_manager), error = function(e) NULL)
            omicQcMetadataSummaryUi(meta, config)
        })
        return(invisible(output))
    }
    output$data_summary <- shiny::renderUI({
        shiny::req(workflow_data$state_manager)
        # KIT-082: the design, column names and id column; no assay table.
        shell <- tryCatch(readWorkflowOmicsShell(binding, workflow_data$state_manager, measures = "features"),
            error = function(e) NULL)
        omicQcDataSummaryUi(shell, config)
    })
    output$assay_stats_table <- DT::renderDT({
        shiny::req(workflow_data$state_manager)
        # KIT-082: counted a page at a time from the store.
        shell <- tryCatch(readWorkflowOmicsShell(binding, workflow_data$state_manager, measures = c("features", "missing")),
            error = function(e) NULL)
        table <- omicQcAssayStatsFrame(shell, config)
        if (is.null(table)) return(NULL)
        DT::datatable(table, options = list(dom = "t", paging = FALSE, ordering = FALSE), rownames = FALSE,
            class = "compact stripe")
    })
    invisible(output)
}

#' The finalise's results text (mod_lipid_qc_s4.R:311-337; trx's `trxQcFinalizeText` has no assay count).
#' @param total The features retained.
#' @param n_assays The number of assays, or NULL for an omic with one table.
#' @param history The owner's history after the commit.
#' @param config An omic GUI configuration.
#' @noRd
omicQcFinalizeText <- function(total, n_assays, history, config) {
    paste(c("QC Finalization Complete", "========================", "",
        sprintf(config$qc$texts[["finalize_total"]], total),
        if (!is.null(n_assays)) sprintf("Number of assays: %d", n_assays),
        sprintf("Processing steps completed: %d", length(history)), "", "Processing History:",
        paste(sprintf("  %d. %s", seq_along(history), history), collapse = "\n"), "",
        sprintf("State saved as: '%s_qc_complete'", config$binding$binding$prefix), "",
        "You can now proceed to the Normalization tab."), collapse = "\n")
}

#' The QC tab completed: `quality_control` complete, and with `finalize_opens_norm` a disabled or locked normalisation
#' tab pending (lipid's `completeLipidQcS4TabStatus`; metab sets QC only, D3).
#' @noRd
completeOmicQcTabStatus <- function(workflow_data, opens_norm) {
    status <- workflow_data$tab_status
    status[["quality_control"]] <- "complete"
    if (isTRUE(opens_norm) && "normalization" %in% names(status) && status[["normalization"]] %in% c("disabled", "locked")) {
        status[["normalization"]] <- "pending"
    }
    workflow_data$tab_status <- status
    invisible(status)
}

#' The QC progress figure at the finalise, from the state (in process) or the worker's metrics; NULL when it cannot
#' be drawn (mod_lipid_qc_s4.R:445-470).
#' @noRd
omicQcTrackingFigure <- function(state, omic_type, plot, config, metrics = NULL) {
    figure <- tryCatch(config$qc$filtering_figure(theObject = state, step_name = "4_QC_Complete", omics_type = omic_type,
        return_grid = TRUE, overwrite = TRUE, metrics = metrics), error = function(e) NULL)
    plot(figure)
    invisible(figure)
}

#' The finalise's success: the results text (in a frame that holds no state, KIT-075), the log and the notification.
#' @noRd
reportOmicQcFinalizeSuccess <- function(text, output, config) {
    output$finalize_results <- shiny::renderText(text)
    logger::log_info(paste(config$shell$label, "QC finalized successfully"))
    shiny::showNotification("QC complete! Proceed to Normalization.", type = "message", duration = 5)
    invisible(text)
}

#' The finalise's refusal or failure, logged and notified.
#' @noRd
reportOmicQcFinalizeError <- function(error) {
    message <- paste("Error finalizing QC:", conditionMessage(error))
    logger::log_error(message)
    shiny::showNotification(message, type = "error")
    invisible(message)
}

#' Finalise the QC: the head's shell checked (`validateWorkflowOmicsQcFinalize`: class, assays, design rows and the
#' omic's alignment), the store step `<prefix>_qc_finalize` committed under `<prefix>_qc_complete`, the progress
#' figure (the worker's metrics on a managed session), the tab completed and the results reported
#' (`runLipidQcS4FinalizeWorkflow`, mod_lipid_qc_s4.R:517-589). A refusal commits nothing.
#' @param workflow_data The workflow bus.
#' @param omic_type The QC tab's omic.
#' @param plot The sub-tab's figure `reactiveVal`.
#' @param output The sub-tab's output.
#' @param config An omic GUI configuration.
#' @noRd
runOmicQcFinalize <- function(workflow_data, omic_type, plot, output, config) {
    binding <- config$binding$binding
    qc <- config$qc
    manager <- workflow_data$state_manager
    tryCatch({
        # KIT-082: the checks read names, classes and row counts, and the text counts the id column. KIT-104 W6
        # (F104-04): an omic whose class cannot be rebuilt as the generic shell reads its own (`qc$finalize_shell`).
        read_shell <- qc$finalize_shell %||% function(binding, manager) {
            readWorkflowOmicsShell(binding, manager, measures = "features")
        }
        shell <- validateWorkflowOmicsQcFinalize(read_shell(binding, manager), qc, binding)
        assays <- methods::slot(shell, binding$data_slot)
        total <- if (isTRUE(qc$finalize_per_assay)) omicQcFeatureTotal(shell, binding) else
            sum(attr(shell, "assay_rows") %||% vapply(assays, nrow, integer(1)))
        n_assays <- if (isTRUE(qc$finalize_per_assay)) length(assays)
        rm(shell, assays)
        commitWorkflowOmicsStoreStep(manager, paste0(binding$prefix, "_qc_finalize"), paste0(binding$prefix, "_qc_complete"),
            workflow_data$config_list, "QC processing complete - ready for normalization")
        finish <- function() {
            completeOmicQcTabStatus(workflow_data, qc$finalize_opens_norm)
            history <- manager$getHistory()
            reportOmicQcFinalizeSuccess(omicQcFinalizeText(total, n_assays, history, config), output, config)
        }
        # KIT-082: the figure's metrics read the whole matrix, so a managed session's worker computes them.
        runner <- workflowOmicsManagedActions(binding$omic_type)
        if (!is.null(runner)) {
            runner$run(list(omicsManagedProduct(paste0(binding$prefix, "_filtering_metrics"), list(), function(result) {
                omicQcTrackingFigure(NULL, omic_type, plot, config,
                    metrics = readWorkflowOmicsFilteringMetrics(manager, result))
            })),
                # As in process, a figure's failure leaves the tab to complete with no figure; the state is committed.
                done = function() tryCatch(finish(), error = function(error) reportOmicQcFinalizeError(error)),
                failed = function(error) {
                    plot(NULL)
                    tryCatch(finish(), error = function(error) reportOmicQcFinalizeError(error))
                })
            return(invisible(NULL))
        }
        omicQcTrackingFigure(config$binding$read_state(manager), omic_type, plot, config)
        finish()
    }, error = function(e) reportOmicQcFinalizeError(e))
}

#' The finalise sub-tab inside its module session: the history, the data summary, Finalize QC and the figure.
#' @noRd
registerOmicQcFinalize <- function(input, output, workflow_data, omic_type, config) {
    plot <- shiny::reactiveVal(NULL)
    registerOmicQcStateHistory(output, workflow_data)
    registerOmicQcFinalizeSummary(output, workflow_data, config)
    finalize_after_start <- qcClickAfterStart(input, "finalize_qc") # KIT-101 F-24
    shiny::observeEvent(input$finalize_qc, {
        if (!finalize_after_start()) return(invisible(NULL))
        shiny::req(workflow_data$state_manager)
        runOmicQcFinalize(workflow_data, omic_type, plot, output, config)
    })
    registerOmicQcPlotOutput(output, plot)
    invisible(NULL)
}
