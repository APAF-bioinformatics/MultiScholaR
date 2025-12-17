# ============================================================================
# mod_metab_de.R
# ============================================================================
# Purpose: Metabolomics differential abundance analysis Shiny module
#
# This module provides per-assay differential abundance analysis using limma,
# with volcano plots and result tables.
# ============================================================================

#' @title Metabolomics Differential Analysis Module
#' @description A Shiny module for performing per-assay differential abundance
#'              analysis on metabolomics data using limma.
#' @name mod_metab_de
NULL

#' @rdname mod_metab_de
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 h4 h5 p hr br selectInput numericInput actionButton uiOutput verbatimTextOutput plotOutput tags downloadButton checkboxInput
#' @importFrom DT DTOutput
#' @importFrom shinyjqui jqui_resizable
mod_metab_de_ui <- function(id) {
    ns <- shiny::NS(id)
    
    shiny::tagList(
        shiny::fluidRow(
            shiny::column(12
                , shiny::wellPanel(
                    shiny::h3("Differential Abundance Analysis")
                    
                    , shiny::fluidRow(
                        # Left column: Parameters
                        shiny::column(4
                            , shiny::h4("Analysis Parameters")
                            
                            # Contrast selection
                            , shiny::selectInput(
                                ns("contrast_select")
                                , "Select Contrast"
                                , choices = NULL
                            )
                            
                            , shiny::hr()
                            
                            # Statistical thresholds
                            , shiny::numericInput(
                                ns("pval_threshold")
                                , "P-value threshold"
                                , value = 0.05
                                , min = 0.001
                                , max = 0.1
                                , step = 0.01
                            )
                            
                            , shiny::numericInput(
                                ns("fc_threshold")
                                , "Log2 fold change threshold"
                                , value = 1.0
                                , min = 0
                                , max = 5
                                , step = 0.25
                            )
                            
                            , shiny::checkboxInput(
                                ns("use_adj_pval")
                                , "Use adjusted p-values (BH)"
                                , value = TRUE
                            )
                            
                            , shiny::hr()
                            
                            , shiny::actionButton(
                                ns("run_de")
                                , "Run Differential Analysis"
                                , class = "btn-primary"
                                , width = "100%"
                                , icon = shiny::icon("chart-line")
                            )
                            
                            , shiny::br()
                            , shiny::br()
                            
                            # Results summary
                            , shiny::uiOutput(ns("de_summary"))
                            
                            , shiny::hr()
                            
                            # Download buttons
                            , shiny::downloadButton(
                                ns("download_results")
                                , "Download Results (CSV)"
                                , class = "btn-secondary"
                                , style = "width: 100%;"
                            )
                        )
                        
                        # Right column: Results
                        , shiny::column(8
                            # Assay tabs for results
                            , shiny::uiOutput(ns("results_tabs"))
                        )
                    )
                )
            )
        )
    )
}

#' @rdname mod_metab_de
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent req renderUI renderPlot renderText showNotification removeNotification updateSelectInput downloadHandler tags tabsetPanel tabPanel
#' @importFrom DT renderDT datatable formatRound formatStyle styleInterval
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline labs theme_minimal scale_color_manual facet_wrap
#' @importFrom logger log_info log_error log_warn
#' @importFrom limma makeContrasts lmFit contrasts.fit eBayes topTable
mod_metab_de_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns
        
        # Local reactive values
        local_data <- shiny::reactiveValues(
            de_results = NULL
            , de_results_list = NULL
            , significant_counts = NULL
        )
        
        # Update contrast dropdown when contrasts are available
        shiny::observe({
            contrasts <- workflow_data$contrasts
            
            if (!is.null(contrasts) && length(contrasts) > 0) {
                shiny::updateSelectInput(
                    session
                    , "contrast_select"
                    , choices = contrasts
                    , selected = contrasts[1]
                )
            }
        })
        
        # Run differential analysis
        shiny::observeEvent(input$run_de, {
            shiny::req(workflow_data$state_manager, input$contrast_select)
            
            shiny::showNotification(
                "Running differential analysis..."
                , id = "de_working"
                , duration = NULL
            )
            
            tryCatch({
                current_s4 <- workflow_data$state_manager$getState()
                shiny::req(current_s4)
                
                if (!inherits(current_s4, "MetaboliteAssayData")) {
                    stop("Current state is not a MetaboliteAssayData object")
                }
                
                assay_list <- current_s4@metabolite_data
                design_matrix <- current_s4@design_matrix
                sample_id_col <- current_s4@sample_id
                group_col <- current_s4@group_id
                metabolite_id_col <- current_s4@metabolite_id_column
                
                # Parse contrast
                contrast_str <- input$contrast_select
                contrast_parts <- strsplit(contrast_str, "-")[[1]]
                
                if (length(contrast_parts) != 2) {
                    stop("Invalid contrast format. Use: GroupA-GroupB")
                }
                
                group_a <- trimws(contrast_parts[1])
                group_b <- trimws(contrast_parts[2])
                
                # Run DE for each assay
                results_list <- list()
                sig_counts <- list()
                
                for (assay_name in names(assay_list)) {
                    logger::log_info(sprintf("Running DE for assay: %s", assay_name))
                    
                    assay_data <- assay_list[[assay_name]]
                    quant_info <- getMetaboliteQuantData(assay_data)
                    
                    # Get expression matrix (samples as columns)
                    expr_matrix <- as.matrix(quant_info$quant_data)
                    rownames(expr_matrix) <- assay_data[[metabolite_id_col]]
                    
                    # Match samples between expression and design
                    sample_cols <- colnames(expr_matrix)
                    dm_matched <- design_matrix[design_matrix[[sample_id_col]] %in% sample_cols, ]
                    
                    if (nrow(dm_matched) == 0) {
                        logger::log_warn(sprintf("No matching samples for assay: %s", assay_name))
                        next
                    }
                    
                    # Reorder expression matrix to match design
                    dm_matched <- dm_matched[order(match(dm_matched[[sample_id_col]], sample_cols)), ]
                    expr_matrix <- expr_matrix[, dm_matched[[sample_id_col]], drop = FALSE]
                    
                    # Create design matrix for limma
                    groups <- factor(dm_matched[[group_col]])
                    design <- model.matrix(~ 0 + groups)
                    colnames(design) <- levels(groups)
                    
                    # Make contrasts
                    if (!group_a %in% colnames(design) || !group_b %in% colnames(design)) {
                        logger::log_warn(sprintf(
                            "Groups '%s' or '%s' not found in design for assay: %s"
                            , group_a, group_b, assay_name
                        ))
                        next
                    }
                    
                    contrast_formula <- paste0(group_a, "-", group_b)
                    contrasts <- limma::makeContrasts(contrasts = contrast_formula, levels = design)
                    
                    # Fit limma model
                    fit <- limma::lmFit(expr_matrix, design)
                    fit2 <- limma::contrasts.fit(fit, contrasts)
                    fit2 <- limma::eBayes(fit2)
                    
                    # Extract results
                    results <- limma::topTable(
                        fit2
                        , coef = 1
                        , number = Inf
                        , adjust.method = "BH"
                        , sort.by = "P"
                    )
                    
                    results$metabolite_id <- rownames(results)
                    results$assay <- assay_name
                    results$contrast <- contrast_str
                    
                    # Classify significance
                    pval_col <- if (input$use_adj_pval) "adj.P.Val" else "P.Value"
                    results$significant <- ifelse(
                        results[[pval_col]] < input$pval_threshold & 
                            abs(results$logFC) >= input$fc_threshold
                        , ifelse(results$logFC > 0, "Up", "Down")
                        , "NS"
                    )
                    
                    results_list[[assay_name]] <- results
                    
                    # Count significant
                    sig_counts[[assay_name]] <- list(
                        up = sum(results$significant == "Up")
                        , down = sum(results$significant == "Down")
                        , ns = sum(results$significant == "NS")
                    )
                }
                
                if (length(results_list) == 0) {
                    stop("No results generated. Check your contrast and sample matching.")
                }
                
                # Combine results
                all_results <- do.call(rbind, results_list)
                rownames(all_results) <- NULL
                
                local_data$de_results <- all_results
                local_data$de_results_list <- results_list
                local_data$significant_counts <- sig_counts
                
                # Update processing log
                workflow_data$processing_log$differential_analysis <- list(
                    timestamp = Sys.time()
                    , contrast = contrast_str
                    , pval_threshold = input$pval_threshold
                    , fc_threshold = input$fc_threshold
                    , significant_counts = sig_counts
                )
                
                workflow_data$tab_status$differential_analysis <- "complete"
                
                logger::log_info(sprintf(
                    "DE analysis complete: %d results across %d assays"
                    , nrow(all_results)
                    , length(results_list)
                ))
                
                shiny::removeNotification("de_working")
                shiny::showNotification("Differential analysis complete!", type = "message")
                
            }, error = function(e) {
                logger::log_error(paste("DE analysis error:", e$message))
                shiny::removeNotification("de_working")
                shiny::showNotification(paste("Error:", e$message), type = "error", duration = 10)
            })
        })
        
        # DE summary
        output$de_summary <- shiny::renderUI({
            sig_counts <- local_data$significant_counts
            
            if (is.null(sig_counts)) {
                return(shiny::p(shiny::icon("info-circle"), " Run analysis to see results.", style = "color: #666;"))
            }
            
            summary_items <- lapply(names(sig_counts), function(assay_name) {
                counts <- sig_counts[[assay_name]]
                shiny::tags$li(
                    shiny::tags$strong(assay_name, ": ")
                    , shiny::tags$span(sprintf("%d ", counts$up), style = "color: #e74c3c;")
                    , shiny::icon("arrow-up", style = "color: #e74c3c;")
                    , " | "
                    , shiny::tags$span(sprintf("%d ", counts$down), style = "color: #3498db;")
                    , shiny::icon("arrow-down", style = "color: #3498db;")
                )
            })
            
            shiny::tagList(
                shiny::h5("Significant Metabolites")
                , shiny::tags$ul(summary_items, style = "list-style: none; padding-left: 0;")
            )
        })
        
        # Results tabs
        output$results_tabs <- shiny::renderUI({
            results_list <- local_data$de_results_list
            
            if (is.null(results_list) || length(results_list) == 0) {
                return(shiny::div(
                    class = "alert alert-info"
                    , shiny::icon("info-circle")
                    , " Run differential analysis to see results."
                ))
            }
            
            tab_list <- list(
                # Combined volcano plot
                shiny::tabPanel(
                    "Volcano Plot (All)"
                    , shiny::br()
                    , shinyjqui::jqui_resizable(
                        shiny::plotOutput(ns("volcano_combined"), height = "500px")
                    )
                )
                # Combined results table
                , shiny::tabPanel(
                    "Results Table"
                    , shiny::br()
                    , DT::DTOutput(ns("results_table"))
                )
            )
            
            # Add per-assay tabs
            for (assay_name in names(results_list)) {
                tab_list <- c(tab_list, list(
                    shiny::tabPanel(
                        assay_name
                        , shiny::br()
                        , shiny::fluidRow(
                            shiny::column(6
                                , shinyjqui::jqui_resizable(
                                    shiny::plotOutput(ns(paste0("volcano_", gsub("[^a-zA-Z0-9]", "_", assay_name))), height = "400px")
                                )
                            )
                            , shiny::column(6
                                , DT::DTOutput(ns(paste0("table_", gsub("[^a-zA-Z0-9]", "_", assay_name))))
                            )
                        )
                    )
                ))
            }
            
            do.call(shiny::tabsetPanel, c(list(id = ns("de_results_tabs")), tab_list))
        })
        
        # Combined volcano plot
        output$volcano_combined <- shiny::renderPlot({
            results <- local_data$de_results
            shiny::req(results)
            
            pval_col <- if (input$use_adj_pval) "adj.P.Val" else "P.Value"
            
            ggplot2::ggplot(results, ggplot2::aes(
                x = logFC
                , y = -log10(get(pval_col))
                , color = significant
            )) +
                ggplot2::geom_point(alpha = 0.6, size = 1.5) +
                ggplot2::geom_hline(
                    yintercept = -log10(input$pval_threshold)
                    , linetype = "dashed"
                    , color = "gray50"
                ) +
                ggplot2::geom_vline(
                    xintercept = c(-input$fc_threshold, input$fc_threshold)
                    , linetype = "dashed"
                    , color = "gray50"
                ) +
                ggplot2::facet_wrap(~assay, scales = "free") +
                ggplot2::scale_color_manual(
                    values = c("Up" = "#e74c3c", "Down" = "#3498db", "NS" = "gray70")
                ) +
                ggplot2::labs(
                    title = paste("Volcano Plot:", input$contrast_select)
                    , x = "Log2 Fold Change"
                    , y = paste("-log10(", pval_col, ")")
                    , color = "Significance"
                ) +
                ggplot2::theme_minimal() +
                ggplot2::theme(legend.position = "bottom")
        })
        
        # Results table
        output$results_table <- DT::renderDT({
            results <- local_data$de_results
            shiny::req(results)
            
            # Select columns to display
            display_cols <- c("metabolite_id", "assay", "logFC", "AveExpr", "P.Value", "adj.P.Val", "significant")
            display_df <- results[, display_cols]
            
            DT::datatable(
                display_df
                , options = list(
                    pageLength = 25
                    , scrollX = TRUE
                    , order = list(list(4, 'asc'))
                )
                , rownames = FALSE
                , filter = "top"
            ) %>%
                DT::formatRound(columns = c("logFC", "AveExpr"), digits = 3) %>%
                DT::formatRound(columns = c("P.Value", "adj.P.Val"), digits = 6) %>%
                DT::formatStyle(
                    "significant"
                    , backgroundColor = DT::styleEqual(
                        c("Up", "Down", "NS")
                        , c("#ffcccc", "#cce5ff", "white")
                    )
                )
        })
        
        # Create per-assay outputs dynamically
        shiny::observe({
            results_list <- local_data$de_results_list
            shiny::req(results_list)
            
            for (assay_name in names(results_list)) {
                local({
                    an <- assay_name
                    safe_name <- gsub("[^a-zA-Z0-9]", "_", an)
                    
                    # Per-assay volcano
                    output[[paste0("volcano_", safe_name)]] <- shiny::renderPlot({
                        results <- results_list[[an]]
                        shiny::req(results)
                        
                        pval_col <- if (input$use_adj_pval) "adj.P.Val" else "P.Value"
                        
                        ggplot2::ggplot(results, ggplot2::aes(
                            x = logFC
                            , y = -log10(get(pval_col))
                            , color = significant
                        )) +
                            ggplot2::geom_point(alpha = 0.6, size = 2) +
                            ggplot2::geom_hline(
                                yintercept = -log10(input$pval_threshold)
                                , linetype = "dashed"
                                , color = "gray50"
                            ) +
                            ggplot2::geom_vline(
                                xintercept = c(-input$fc_threshold, input$fc_threshold)
                                , linetype = "dashed"
                                , color = "gray50"
                            ) +
                            ggplot2::scale_color_manual(
                                values = c("Up" = "#e74c3c", "Down" = "#3498db", "NS" = "gray70")
                            ) +
                            ggplot2::labs(
                                title = paste("Volcano Plot:", an)
                                , x = "Log2 Fold Change"
                                , y = paste("-log10(", pval_col, ")")
                            ) +
                            ggplot2::theme_minimal()
                    })
                    
                    # Per-assay table
                    output[[paste0("table_", safe_name)]] <- DT::renderDT({
                        results <- results_list[[an]]
                        shiny::req(results)
                        
                        # Top significant
                        sig_results <- results[results$significant != "NS", ]
                        sig_results <- sig_results[order(sig_results$P.Value), ]
                        
                        display_cols <- c("metabolite_id", "logFC", "P.Value", "adj.P.Val", "significant")
                        
                        DT::datatable(
                            sig_results[, display_cols]
                            , options = list(
                                pageLength = 10
                                , dom = 'ftp'
                            )
                            , rownames = FALSE
                        ) %>%
                            DT::formatRound(columns = c("logFC"), digits = 3) %>%
                            DT::formatRound(columns = c("P.Value", "adj.P.Val"), digits = 6)
                    })
                })
            }
        })
        
        # Download handler
        output$download_results <- shiny::downloadHandler(
            filename = function() {
                paste0("metabolomics_de_results_", Sys.Date(), ".csv")
            }
            , content = function(file) {
                results <- local_data$de_results
                shiny::req(results)
                write.csv(results, file, row.names = FALSE)
            }
        )
    })
}

