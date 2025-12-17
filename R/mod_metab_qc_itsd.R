# ============================================================================
# mod_metab_qc_itsd.R
# ============================================================================
# Purpose: Internal standard (ITSD) QC visualization Shiny module
#
# This module provides visualization of internal standard performance metrics
# including CV distribution and intensity trends across samples.
# ============================================================================

#' @title Internal Standard QC Visualization Module
#' @description A Shiny module for visualizing internal standard quality metrics.
#'              Displays IS detection, CV distributions, and intensity trends.
#' @name mod_metab_qc_itsd
NULL

#' @rdname mod_metab_qc_itsd
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 h5 p hr textInput actionButton verbatimTextOutput uiOutput plotOutput
#' @importFrom shinyjqui jqui_resizable
mod_metab_qc_itsd_ui <- function(id) {
    ns <- shiny::NS(id)
    
    shiny::tabPanel(
        "Internal Standards"
        , shiny::br()
        , shiny::fluidRow(
            shiny::column(4
                , shiny::wellPanel(
                    shiny::h4("Internal Standard QC")
                    , shiny::p("Visualize internal standard (IS) performance across samples. Good IS should have low CV and consistent intensities.")
                    , shiny::hr()
                    
                    # IS pattern input
                    , shiny::textInput(
                        ns("is_pattern")
                        , "Internal Standard Regex Pattern"
                        , value = ""
                        , placeholder = "e.g., ^IS_|_d[0-9]+$|ISTD"
                    )
                    , shiny::helpText("Regular expression to identify internal standards in metabolite IDs. Leave empty to use pattern from S4 object.")
                    
                    , shiny::hr()
                    , shiny::actionButton(
                        ns("analyze_is")
                        , "Analyze Internal Standards"
                        , class = "btn-primary"
                        , width = "100%"
                    )
                    
                    , shiny::hr()
                    , shiny::h5("IS Detection Summary")
                    , shiny::uiOutput(ns("is_summary"))
                )
            )
            , shiny::column(8
                , shiny::verbatimTextOutput(ns("is_results"))
                , shiny::br()
                # Per-assay IS visualization tabs
                , shiny::uiOutput(ns("is_viz_tabs"))
            )
        )
    )
}

#' @rdname mod_metab_qc_itsd
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification renderText renderUI renderPlot tabsetPanel tabPanel tags
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_point geom_line labs theme_minimal theme element_text coord_flip scale_fill_brewer facet_wrap geom_hline
#' @importFrom logger log_info log_error log_warn
mod_metab_qc_itsd_server <- function(id, workflow_data, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns
        
        is_metrics <- shiny::reactiveVal(NULL)
        is_data_long <- shiny::reactiveVal(NULL)
        
        # Analyze internal standards
        shiny::observeEvent(input$analyze_is, {
            shiny::req(workflow_data$state_manager)
            
            shiny::showNotification(
                "Analyzing internal standards..."
                , id = "is_analysis_working"
                , duration = NULL
            )
            
            tryCatch({
                current_s4 <- workflow_data$state_manager$getState()
                shiny::req(current_s4)
                
                if (!inherits(current_s4, "MetaboliteAssayData")) {
                    stop("Current state is not a MetaboliteAssayData object")
                }
                
                # Get IS pattern (from input or S4 object)
                is_pattern <- if (nzchar(input$is_pattern)) {
                    input$is_pattern
                } else if (!is.na(current_s4@internal_standard_regex) && nzchar(current_s4@internal_standard_regex)) {
                    current_s4@internal_standard_regex
                } else {
                    NULL
                }
                
                if (is.null(is_pattern) || !nzchar(is_pattern)) {
                    stop("No internal standard pattern provided. Please enter a regex pattern to identify IS.")
                }
                
                logger::log_info(paste("Analyzing internal standards with pattern:", is_pattern))
                
                metabolite_id_col <- current_s4@metabolite_id_column
                sample_id_col <- current_s4@sample_id
                assay_list <- current_s4@metabolite_data
                assay_names <- names(assay_list)
                
                # Calculate IS metrics per assay
                metrics_list <- list()
                long_data_list <- list()
                
                for (i in seq_along(assay_list)) {
                    assay_data <- assay_list[[i]]
                    assay_name <- if (!is.null(assay_names)) assay_names[i] else paste0("Assay_", i)
                    
                    # Get IS metrics using existing function
                    is_met <- getInternalStandardMetrics(
                        assay_data = assay_data
                        , is_pattern = is_pattern
                        , metabolite_id_col = metabolite_id_col
                        , sample_id_col = sample_id_col
                    )
                    
                    if (nrow(is_met) > 0) {
                        is_met$assay <- assay_name
                        metrics_list[[assay_name]] <- is_met
                        
                        # Also extract long-format IS data for intensity plots
                        quant_info <- getMetaboliteQuantData(assay_data)
                        sample_cols <- quant_info$sample_names
                        
                        if (length(sample_cols) > 0) {
                            # Filter to IS rows
                            is_ids <- is_met$is_id
                            is_rows <- assay_data[[metabolite_id_col]] %in% is_ids
                            
                            if (any(is_rows)) {
                                is_data <- assay_data[is_rows, c(metabolite_id_col, sample_cols), drop = FALSE]
                                
                                # Pivot to long format
                                is_long <- tidyr::pivot_longer(
                                    is_data
                                    , cols = dplyr::all_of(sample_cols)
                                    , names_to = "Sample"
                                    , values_to = "Intensity"
                                )
                                is_long$assay <- assay_name
                                names(is_long)[1] <- "IS_ID"
                                long_data_list[[assay_name]] <- is_long
                            }
                        }
                    }
                }
                
                if (length(metrics_list) == 0) {
                    stop(paste("No internal standards found matching pattern:", is_pattern))
                }
                
                # Combine results
                all_metrics <- do.call(rbind, metrics_list)
                is_metrics(all_metrics)
                
                if (length(long_data_list) > 0) {
                    all_long <- do.call(rbind, long_data_list)
                    is_data_long(all_long)
                }
                
                # Generate summary text
                n_is_total <- nrow(all_metrics)
                median_cv <- median(all_metrics$cv, na.rm = TRUE)
                max_cv <- max(all_metrics$cv, na.rm = TRUE)
                
                result_text <- paste(
                    "Internal Standard Analysis Complete"
                    , "===================================="
                    , sprintf("Pattern used: %s", is_pattern)
                    , sprintf("Total IS detected: %d", n_is_total)
                    , ""
                    , "CV Statistics:"
                    , sprintf("  Median CV: %.1f%%", median_cv)
                    , sprintf("  Max CV: %.1f%%", max_cv)
                    , sprintf("  IS with CV > 30%%: %d", sum(all_metrics$cv > 30, na.rm = TRUE))
                    , ""
                    , "Per-Assay IS Counts:"
                    , paste(sapply(names(metrics_list), function(name) {
                        sprintf("  %s: %d internal standards", name, nrow(metrics_list[[name]]))
                    }), collapse = "\n")
                    , sep = "\n"
                )
                
                output$is_results <- shiny::renderText(result_text)
                
                logger::log_info(paste("IS analysis complete:", n_is_total, "standards found"))
                
                shiny::removeNotification("is_analysis_working")
                shiny::showNotification(
                    sprintf("Found %d internal standards", n_is_total)
                    , type = "message"
                )
                
            }, error = function(e) {
                msg <- paste("Error analyzing internal standards:", e$message)
                logger::log_error(msg)
                shiny::showNotification(msg, type = "error", duration = 10)
                shiny::removeNotification("is_analysis_working")
            })
        })
        
        # Render IS summary
        output$is_summary <- shiny::renderUI({
            metrics <- is_metrics()
            
            if (is.null(metrics)) {
                return(shiny::p(
                    shiny::icon("info-circle")
                    , " Click 'Analyze' to detect internal standards."
                    , style = "color: #666;"
                ))
            }
            
            # Summarize by assay
            summary_by_assay <- split(metrics, metrics$assay)
            
            summary_items <- lapply(names(summary_by_assay), function(assay_name) {
                assay_metrics <- summary_by_assay[[assay_name]]
                n_is <- nrow(assay_metrics)
                median_cv <- median(assay_metrics$cv, na.rm = TRUE)
                
                # Color based on CV quality
                cv_color <- if (median_cv <= 15) "green" else if (median_cv <= 30) "orange" else "red"
                cv_icon <- if (median_cv <= 15) "check-circle" else if (median_cv <= 30) "exclamation-circle" else "times-circle"
                
                shiny::tags$li(
                    shiny::icon(cv_icon, style = paste("color:", cv_color))
                    , sprintf(" %s: %d IS (median CV: %.1f%%)", assay_name, n_is, median_cv)
                )
            })
            
            shiny::tagList(
                shiny::tags$ul(summary_items, style = "list-style: none; padding-left: 0;")
                , shiny::hr()
                , shiny::tags$small(
                    shiny::icon("info-circle")
                    , " CV < 15%: Good | 15-30%: Acceptable | > 30%: Review"
                    , style = "color: #666;"
                )
            )
        })
        
        # Render IS visualization tabs
        output$is_viz_tabs <- shiny::renderUI({
            metrics <- is_metrics()
            
            if (is.null(metrics)) {
                return(NULL)
            }
            
            assays <- unique(metrics$assay)
            
            tab_list <- list(
                # CV Distribution tab
                shiny::tabPanel(
                    "CV Distribution"
                    , shiny::br()
                    , shinyjqui::jqui_resizable(
                        shiny::plotOutput(ns("cv_plot"), height = "500px")
                    )
                )
                # Intensity Trends tab
                , shiny::tabPanel(
                    "Intensity Trends"
                    , shiny::br()
                    , shinyjqui::jqui_resizable(
                        shiny::plotOutput(ns("intensity_plot"), height = "500px")
                    )
                )
            )
            
            do.call(shiny::tabsetPanel, c(list(id = ns("is_viz_tabset")), tab_list))
        })
        
        # Render CV distribution plot
        output$cv_plot <- shiny::renderPlot({
            metrics <- is_metrics()
            shiny::req(metrics)
            
            ggplot2::ggplot(metrics, ggplot2::aes(x = is_id, y = cv, fill = assay)) +
                ggplot2::geom_bar(stat = "identity", position = "dodge") +
                ggplot2::geom_hline(yintercept = 15, linetype = "dashed", color = "green", linewidth = 1) +
                ggplot2::geom_hline(yintercept = 30, linetype = "dashed", color = "orange", linewidth = 1) +
                ggplot2::labs(
                    title = "Internal Standard CV Distribution"
                    , x = "Internal Standard ID"
                    , y = "Coefficient of Variation (%)"
                    , fill = "Assay"
                ) +
                ggplot2::theme_minimal() +
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8)
                    , legend.position = "bottom"
                ) +
                ggplot2::scale_fill_brewer(palette = "Set1") +
                ggplot2::coord_flip()
        })
        
        # Render intensity trend plot
        output$intensity_plot <- shiny::renderPlot({
            long_data <- is_data_long()
            shiny::req(long_data)
            
            ggplot2::ggplot(long_data, ggplot2::aes(x = Sample, y = Intensity, color = IS_ID, group = IS_ID)) +
                ggplot2::geom_line(alpha = 0.7) +
                ggplot2::geom_point(size = 1, alpha = 0.5) +
                ggplot2::facet_wrap(~assay, scales = "free_y") +
                ggplot2::labs(
                    title = "Internal Standard Intensity Across Samples"
                    , x = "Sample"
                    , y = "Intensity"
                    , color = "Internal Standard"
                ) +
                ggplot2::theme_minimal() +
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)
                    , legend.position = "bottom"
                )
        })
    })
}

