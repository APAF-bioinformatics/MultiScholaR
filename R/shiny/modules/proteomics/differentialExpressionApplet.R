#' @title differentialExpressionAppletModule
#'
#' @description A Shiny module for the Differential Expression step of the proteomics
#' workflow. Handles limma-based statistical analysis and visualization of results
#' across user-defined contrasts.
#'
#' @name differentialExpressionAppletModule
NULL

#' @rdname differentialExpressionAppletModule
#' @export
#' @import shiny
#' @import shinydashboard
differentialExpressionAppletUI <- function(id) {
  ns <- NS(id)
  
  shiny::fluidRow(
    shiny::column(3,
      shiny::wellPanel(
        shiny::h4("DE Analysis Settings"),
        
        # Formula input (will be populated from S4 @args)
        shiny::textAreaInput(
          ns("formula_string"),
          "Model Formula:",
          value = "~ 0 + group",
          height = "60px"
        ),
        shiny::helpText("Formula fetched from S4 object @args"),
        
        shiny::hr(),
        
        # Analysis parameters
        shiny::h5("Analysis Parameters"),
        shiny::numericInput(
          ns("de_q_val_thresh"),
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
        
        shiny::hr(),
        
        # Available contrasts display
        shiny::h5("Available Contrasts"),
        shiny::verbatimTextOutput(ns("contrasts_display")),
        shiny::helpText("Contrasts defined in design matrix"),
        
        shiny::br(),
        
        # Main action button
        shiny::actionButton(
          ns("run_de_analysis"),
          "Run DE Analysis",
          class = "btn-primary",
          width = "100%",
          icon = shiny::icon("play")
        ),
        
        shiny::br(),
        shiny::br(),
        
        # Download button
        shiny::downloadButton(
          ns("download_de_results"),
          "Download All Results",
          class = "btn-success",
          width = "100%"
        ),
        
        shiny::br(),
        shiny::br(),
        
        # Status display
        shiny::h5("Analysis Status"),
        shiny::verbatimTextOutput(ns("de_status"))
      )
    ),
    
    # Main results panel with 3 tabs
    shiny::column(9,
      shiny::tabsetPanel(
        id = ns("de_results_tabs"),
        
        # Tab 1: Volcano Plot (Interactive with Glimma)
        shiny::tabPanel(
          "Volcano Plot",
          icon = shiny::icon("chart-scatter"),
          shiny::br(),
          
          # Contrast selector for volcano plot
          shiny::fluidRow(
            shiny::column(6,
              shiny::selectInput(
                ns("volcano_contrast"),
                "Select Contrast:",
                choices = NULL,
                width = "100%"
              )
            ),
            shiny::column(6,
              shiny::checkboxInput(
                ns("volcano_interactive"),
                "Interactive Plot (Glimma)",
                value = TRUE
              )
            )
          ),
          
          shiny::hr(),
          
          # Volcano plot output (will switch between static and interactive)
          shiny::conditionalPanel(
            condition = "input.volcano_interactive == true",
            ns = ns,
            shiny::div(
              id = ns("volcano_glimma_container"),
              style = "height: 600px;",
              shiny::htmlOutput(ns("volcano_glimma"))
            )
          ),
          
          shiny::conditionalPanel(
            condition = "input.volcano_interactive == false", 
            ns = ns,
            shinyjqui::jqui_resizable(
              shiny::plotOutput(ns("volcano_static"), height = "600px")
            )
          )
        ),
        
        # Tab 2: Heatmap (Interactive)
        shiny::tabPanel(
          "Heatmap",
          icon = shiny::icon("th"),
          shiny::br(),
          
          # Contrast selector and heatmap options
          shiny::fluidRow(
            shiny::column(3,
              shiny::selectInput(
                ns("heatmap_contrast"),
                "Select Contrast:",
                choices = NULL,
                width = "100%"
              )
            ),
            shiny::column(3,
              shiny::numericInput(
                ns("heatmap_top_n"),
                "Top N Genes:",
                value = 50,
                min = 10,
                max = 500,
                step = 10
              )
            ),
            shiny::column(3,
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
            shiny::column(3,
              shiny::selectInput(
                ns("heatmap_scaling"),
                "Data Scaling:",
                choices = list(
                  "Row (gene) scaling" = "row",
                  "Column (sample) scaling" = "column",
                  "Both" = "both",
                  "None" = "none"
                ),
                selected = "row"
              )
            )
          ),
          
          # Advanced clustering controls (collapsible)
          shiny::conditionalPanel(
            condition = "input.heatmap_clustering != 'none'",
            ns = ns,
            shiny::wellPanel(
              shiny::h5("Advanced Clustering Options"),
              
              shiny::fluidRow(
                shiny::column(3,
                  shiny::selectInput(
                    ns("heatmap_cluster_method"),
                    "Clustering Method:",
                    choices = list(
                      "Ward (minimum variance)" = "ward.D2",
                      "Ward (original)" = "ward.D",
                      "Complete linkage" = "complete",
                      "Single linkage" = "single", 
                      "Average linkage" = "average",
                      "McQuitty (WPGMA)" = "mcquitty",
                      "Median (WPGMC)" = "median",
                      "Centroid (UPGMC)" = "centroid"
                    ),
                    selected = "ward.D2"
                  )
                ),
                shiny::column(3,
                  shiny::selectInput(
                    ns("heatmap_distance_method"),
                    "Distance Metric:",
                    choices = list(
                      "Euclidean" = "euclidean",
                      "Manhattan" = "manhattan",
                      "Pearson correlation" = "pearson",
                      "Spearman correlation" = "spearman",
                      "Maximum" = "maximum",
                      "Canberra" = "canberra",
                      "Binary" = "binary",
                      "Minkowski" = "minkowski"
                    ),
                    selected = "euclidean"
                  )
                ),
                shiny::column(3,
                  shiny::selectInput(
                    ns("heatmap_tree_cut_method"),
                    "Tree Cutting:",
                    choices = list(
                      "Number of clusters" = "k_clusters",
                      "Height cutoff" = "height_cutoff",
                      "Dynamic tree cutting" = "dynamic",
                      "No cutting" = "none"
                    ),
                    selected = "k_clusters"
                  )
                ),
                shiny::column(3,
                  # Dynamic UI for tree cutting parameters
                  shiny::conditionalPanel(
                    condition = "input.heatmap_tree_cut_method == 'k_clusters'",
                    ns = ns,
                    shiny::numericInput(
                      ns("heatmap_n_clusters"),
                      "Number of Clusters:",
                      value = 4,
                      min = 2,
                      max = 20,
                      step = 1
                    )
                  ),
                  shiny::conditionalPanel(
                    condition = "input.heatmap_tree_cut_method == 'height_cutoff'",
                    ns = ns,
                    shiny::numericInput(
                      ns("heatmap_cut_height"),
                      "Cut Height:",
                      value = 0.5,
                      min = 0.1,
                      max = 2.0,
                      step = 0.1
                    )
                  ),
                  shiny::conditionalPanel(
                    condition = "input.heatmap_tree_cut_method == 'dynamic'",
                    ns = ns,
                    shiny::numericInput(
                      ns("heatmap_min_cluster_size"),
                      "Min Cluster Size:",
                      value = 3,
                      min = 2,
                      max = 20,
                      step = 1
                    )
                  )
                )
              ),
              
              # Additional heatmap display options
              shiny::fluidRow(
                shiny::column(4,
                  shiny::checkboxInput(
                    ns("heatmap_show_dendro"),
                    "Show Dendrogram",
                    value = TRUE
                  )
                ),
                shiny::column(4,
                  shiny::checkboxInput(
                    ns("heatmap_show_labels"),
                    "Show Gene Labels",
                    value = FALSE
                  )
                ),
                shiny::column(4,
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
              )
            )
          ),
          
          shiny::hr(),
          
          # Interactive heatmap output
          shiny::div(
            id = ns("heatmap_container"),
            style = "height: 650px;",
            shinyjqui::jqui_resizable(
              shiny::plotOutput(ns("heatmap_plot"), height = "600px")
            )
          ),
          
          # Cluster summary (if tree cutting is applied)
          shiny::conditionalPanel(
            condition = "input.heatmap_tree_cut_method != 'none'",
            ns = ns,
            shiny::wellPanel(
              shiny::h5("Cluster Information"),
              shiny::verbatimTextOutput(ns("cluster_summary"))
            )
          )
        ),
        
        # Tab 3: DE Results Table (Interactive)
        shiny::tabPanel(
          "DE Results Table",
          icon = shiny::icon("table"),
          shiny::br(),
          
          # Contrast selector and table options
          shiny::fluidRow(
            shiny::column(4,
              shiny::selectInput(
                ns("table_contrast"),
                "Select Contrast:",
                choices = NULL,
                width = "100%"
              )
            ),
            shiny::column(4,
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
            shiny::column(4,
              shiny::numericInput(
                ns("table_max_rows"),
                "Max Rows:",
                value = 1000,
                min = 100,
                max = 10000,
                step = 100
              )
            )
          ),
          
          shiny::hr(),
          
          # Summary statistics
          shiny::fluidRow(
            shiny::column(12,
              shiny::wellPanel(
                shiny::h5("Summary Statistics"),
                shiny::verbatimTextOutput(ns("de_summary_stats"))
              )
            )
          ),
          
          # Interactive results table
          DT::DTOutput(ns("de_results_table"))
        )
      )
    )
  )
}

#' @rdname differentialExpressionAppletModule 
#' @export
differentialExpressionAppletServer <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    message("=== DIFFERENTIAL EXPRESSION MODULE SERVER STARTED ===")
    message(sprintf("Module ID: %s", id))
    
    # Initialize reactive values for DE state
    de_data <- reactiveValues(
      de_results_list = NULL,
      contrasts_available = NULL,
      analysis_complete = FALSE,
      current_s4_object = NULL,
      formula_from_s4 = NULL,
      current_row_clusters = NULL,
      current_col_clusters = NULL
    )
    
    # Fetch formula and contrasts from S4 object when module loads
    observe({
      message("=== DE TAB: Checking for S4 object and contrasts ===")
      
      # Get current S4 object from state manager
      if (!is.null(workflow_data$state_manager)) {
        current_state <- workflow_data$state_manager$current_state
        message(sprintf("   DE TAB: Current state = %s", current_state))
        
        # Should be getting the correlation-filtered protein object (final state before DE)
        if (current_state == "correlation_filtered") {
          message("   DE TAB: State is valid for DE analysis (correlation-filtered state found)")
          current_s4 <- workflow_data$state_manager$getState(current_state)
          
          if (!is.null(current_s4)) {
            message(sprintf("   DE TAB: S4 object retrieved, class = %s", class(current_s4)))
            de_data$current_s4_object <- current_s4
            
            # Extract formula from S4 @args
            if ("deAnalysisParameters" %in% names(current_s4@args)) {
              if ("formula_string" %in% names(current_s4@args$deAnalysisParameters)) {
                formula_from_s4 <- current_s4@args$deAnalysisParameters$formula_string
                de_data$formula_from_s4 <- formula_from_s4
                message(sprintf("   DE TAB: Formula from S4 = %s", formula_from_s4))
                
                # Update UI with formula from S4
                shiny::updateTextAreaInput(
                  session,
                  "formula_string",
                  value = formula_from_s4
                )
              } else {
                message("   DE TAB: No formula_string in deAnalysisParameters")
              }
            } else {
              message("   DE TAB: No deAnalysisParameters in S4 @args")
            }
            
            # Get contrasts from design matrix or S4 args
            if (!is.null(current_s4@design_matrix)) {
              message("   DE TAB: Design matrix found in S4 object")
              message(sprintf("   DE TAB: Design matrix dims = %d rows, %d cols", nrow(current_s4@design_matrix), ncol(current_s4@design_matrix)))
              
              # Check for contrasts_tbl in global environment
              if (exists("contrasts_tbl", envir = .GlobalEnv)) {
                contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
                                 message("   DE TAB: Found contrasts_tbl in global environment")
                 message("   DE TAB: contrasts_tbl structure:")
                 utils::str(contrasts_tbl)
                
                if ("comparison" %in% names(contrasts_tbl)) {
                  de_data$contrasts_available <- contrasts_tbl$comparison
                  message(sprintf("   DE TAB: Set contrasts from comparison column: %s", paste(de_data$contrasts_available, collapse = ", ")))
                } else if ("contrasts" %in% names(contrasts_tbl)) {
                  de_data$contrasts_available <- contrasts_tbl$contrasts
                  message(sprintf("   DE TAB: Set contrasts from contrasts column: %s", paste(de_data$contrasts_available, collapse = ", ")))
                } else {
                  message("   DE TAB: contrasts_tbl found but no recognized column names")
                  # Try first column
                  de_data$contrasts_available <- contrasts_tbl[[1]]
                  message(sprintf("   DE TAB: Using first column: %s", paste(de_data$contrasts_available, collapse = ", ")))
                }
              } else {
                message("   DE TAB: No contrasts_tbl in global environment, creating basic contrasts")
                # Create basic contrasts from groups
                groups <- unique(current_s4@design_matrix$group)
                message(sprintf("   DE TAB: Found groups in design matrix: %s", paste(groups, collapse = ", ")))
                
                if (length(groups) >= 2) {
                  # Create pairwise contrasts (placeholder)
                  basic_contrasts <- combn(groups, 2, function(x) paste(x[1], "-", x[2]), simplify = TRUE)
                  de_data$contrasts_available <- basic_contrasts
                  message(sprintf("   DE TAB: Created basic contrasts: %s", paste(basic_contrasts, collapse = ", ")))
                } else {
                  message("   DE TAB: Not enough groups to create contrasts")
                }
              }
            } else {
              message("   DE TAB: No design matrix found in S4 object")
            }
          } else {
            message("   DE TAB: S4 object is NULL")
          }
        } else {
          message(sprintf("   DE TAB: State '%s' not valid for DE analysis (expecting: correlation_filtered)", current_state))
        }
      } else {
        message("   DE TAB: workflow_data$state_manager is NULL")
      }
      
      message("=== DE TAB: Contrast detection complete ===")
    })
    
    # Update contrast dropdowns when contrasts become available
    observe({
      if (!is.null(de_data$contrasts_available)) {
        contrast_choices <- setNames(de_data$contrasts_available, de_data$contrasts_available)
        
        # Update all contrast selectors
        shiny::updateSelectInput(session, "volcano_contrast", choices = contrast_choices)
        shiny::updateSelectInput(session, "heatmap_contrast", choices = contrast_choices)
        shiny::updateSelectInput(session, "table_contrast", choices = contrast_choices)
      }
    })
    
    # Display available contrasts
    output$contrasts_display <- shiny::renderText({
      if (!is.null(de_data$contrasts_available)) {
        paste(de_data$contrasts_available, collapse = "\n")
      } else {
        "No contrasts available.\nComplete normalization and\ncorrelation filtering first."
      }
    })
    
    # Display analysis status
    output$de_status <- shiny::renderText({
      if (de_data$analysis_complete) {
        paste(
          "✅ Analysis Complete\n",
          sprintf("Contrasts analyzed: %d\n", length(de_data$contrasts_available)),
          "Results available in all tabs"
        )
      } else {
        "⏳ Waiting for analysis...\nClick 'Run DE Analysis' to start"
      }
    })
    
    # Main DE Analysis Button Logic
    observeEvent(input$run_de_analysis, {
      message("=== STARTING DIFFERENTIAL EXPRESSION ANALYSIS ===")
      
      shiny::req(de_data$current_s4_object, de_data$contrasts_available)
      
      shiny::showNotification("Running differential expression analysis...", id = "de_working", duration = NULL)
      
      tryCatch({
        shiny::withProgress(message = "Running DE analysis...", value = 0, {
          
          # Step 1: Prepare contrasts table for analysis
          shiny::incProgress(0.2, detail = "Preparing contrasts for analysis...")
          
          # Get contrasts from global environment or create basic ones
          contrasts_tbl <- NULL
          if (exists("contrasts_tbl", envir = .GlobalEnv)) {
            contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
          } else {
            # Create basic contrasts from available contrasts
            contrasts_tbl <- data.frame(
              contrasts = de_data$contrasts_available,
              stringsAsFactors = FALSE
            )
          }
          
          # Step 2: Run actual differential expression analysis
          shiny::incProgress(0.4, detail = "Running limma analysis...")
          
          # Use the new modular DE analysis function
          de_results_list <- differentialExpressionAnalysis(
            theObject = de_data$current_s4_object,
            contrasts_tbl = contrasts_tbl,
            formula_string = input$formula_string,
            de_q_val_thresh = input$de_q_val_thresh,
            treat_lfc_cutoff = input$treat_lfc_cutoff,
            qvalue_column = "fdr_qvalue",
            raw_pvalue_column = "raw_pvalue"
          )
          
          # Step 3: Store results
          shiny::incProgress(0.8, detail = "Processing results...")
          
          de_data$de_results_list <- de_results_list
          de_data$analysis_complete <- TRUE
          
          # Update workflow data
          workflow_data$de_analysis_results_list <- de_results_list
          workflow_data$tab_status$differential_expression <- "complete"
          workflow_data$tab_status$enrichment_analysis <- "pending"
          
          shiny::incProgress(1.0, detail = "Complete!")
        })
        
        shiny::showNotification(
          "Differential expression analysis completed successfully!",
          type = "success",
          duration = 5
        )
        
        message("=== DIFFERENTIAL EXPRESSION ANALYSIS COMPLETED ===")
        
      }, error = function(e) {
        message(sprintf("*** ERROR in DE analysis: %s", e$message))
        shiny::showNotification(
          sprintf("Error in DE analysis: %s", e$message),
          type = "error",
          duration = 10
        )
      })
      
      shiny::removeNotification("de_working")
    })
    
    # Volcano Plot Rendering
    output$volcano_glimma <- shiny::renderUI({
      shiny::req(input$volcano_contrast, de_data$de_results_list)
      
      tryCatch({
        # Get UniProt annotations if available
        uniprot_tbl <- NULL
        if (exists("uniprot_dat_cln", envir = .GlobalEnv)) {
          uniprot_tbl <- get("uniprot_dat_cln", envir = .GlobalEnv)
        }
        
        # Generate interactive Glimma volcano plot using new function
        glimma_widget <- generateVolcanoPlotGlimma(
          de_results_list = de_data$de_results_list,
          selected_contrast = input$volcano_contrast,
          uniprot_tbl = uniprot_tbl,
          de_q_val_thresh = input$de_q_val_thresh
        )
        
        if (!is.null(glimma_widget)) {
          # Return the Glimma widget
          glimma_widget
        } else {
          shiny::div(
            style = "text-align: center; padding: 20px;",
            shiny::h5("No data available for selected contrast"),
            shiny::p("Please check contrast selection and ensure analysis has completed.")
          )
        }
        
      }, error = function(e) {
        message(sprintf("*** ERROR in volcano plot generation: %s", e$message))
        shiny::div(
          style = "text-align: center; padding: 20px;",
          shiny::h5("Error generating volcano plot"),
          shiny::p(sprintf("Error: %s", e$message))
        )
      })
    })
    
    output$volcano_static <- shiny::renderPlot({
      shiny::req(input$volcano_contrast, de_data$de_results_list)
      
      # Placeholder for static volcano plot
      plot(1:10, 1:10, main = paste("Static Volcano Plot -", input$volcano_contrast),
           xlab = "Log2 Fold Change", ylab = "-Log10 P-value")
    })
    
    # Heatmap Rendering
    output$heatmap_plot <- shiny::renderPlot({
      shiny::req(input$heatmap_contrast, de_data$de_results_list)
      
      tryCatch({
        # Generate heatmap using new function
        heatmap_plot <- generateDEHeatmap(
          de_results_list = de_data$de_results_list,
          selected_contrast = input$heatmap_contrast,
          top_n_genes = input$heatmap_top_n,
          clustering_method = input$heatmap_cluster_method,
          distance_method = input$heatmap_distance_method,
          cluster_rows = input$heatmap_clustering %in% c("both", "row"),
          cluster_cols = input$heatmap_clustering %in% c("both", "column"),
          scale_data = input$heatmap_scaling,
          color_scheme = input$heatmap_color_scheme,
          show_gene_names = input$heatmap_show_labels,
          de_q_val_thresh = input$de_q_val_thresh
        )
        
        if (!is.null(heatmap_plot)) {
          if (is.function(heatmap_plot)) {
            # For base R heatmap functions
            heatmap_plot()
          } else {
            # For ComplexHeatmap objects
            heatmap_plot
          }
        } else {
          # No significant genes found
          plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
               main = paste("No significant genes found for contrast:", input$heatmap_contrast))
          text(1, 1, "Adjust significance thresholds\nor select different contrast", cex = 1.2)
        }
        
      }, error = function(e) {
        message(sprintf("*** ERROR in heatmap generation: %s", e$message))
        plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
             main = "Error generating heatmap")
        text(1, 1, sprintf("Error: %s", e$message), cex = 1.2)
      })
    })
    
    # Cluster Summary Output
    output$cluster_summary <- shiny::renderText({
      if (!is.null(de_data$current_row_clusters) && !is.null(de_data$current_col_clusters)) {
        row_cluster_summary <- table(de_data$current_row_clusters)
        col_cluster_summary <- table(de_data$current_col_clusters)
        
        paste0(
          "Gene Clusters:\n",
          paste(sprintf("Cluster %s: %d genes", names(row_cluster_summary), row_cluster_summary), collapse = "\n"),
          "\n\nSample Clusters:\n",
          paste(sprintf("Cluster %s: %d samples", names(col_cluster_summary), col_cluster_summary), collapse = "\n"),
          "\n\nClustering Method: ", input$heatmap_cluster_method,
          "\nDistance Metric: ", input$heatmap_distance_method,
          "\nTree Cutting: ", input$heatmap_tree_cut_method
        )
      } else {
        "No cluster information available.\nEnable tree cutting to see cluster details."
      }
    })
    
    # DE Results Table Rendering
    output$de_results_table <- DT::renderDT({
      shiny::req(input$table_contrast, de_data$de_results_list)
      
      tryCatch({
        # Get DE results from the new format
        if (!is.null(de_data$de_results_list$de_proteins_long)) {
          # Filter for selected contrast
          current_results <- de_data$de_results_list$de_proteins_long |>
            dplyr::filter(comparison == input$table_contrast)
          
          if (nrow(current_results) > 0) {
            # Filter based on significance selection
            if (input$table_significance == "significant") {
              current_results <- current_results[current_results$fdr_qvalue < input$de_q_val_thresh, ]
            } else if (input$table_significance == "up") {
              current_results <- current_results[current_results$fdr_qvalue < input$de_q_val_thresh & 
                                               current_results$log2FC > input$treat_lfc_cutoff, ]
            } else if (input$table_significance == "down") {
              current_results <- current_results[current_results$fdr_qvalue < input$de_q_val_thresh & 
                                               current_results$log2FC < -input$treat_lfc_cutoff, ]
            }
            
            # Limit rows
            if (nrow(current_results) > input$table_max_rows) {
              current_results <- current_results[1:input$table_max_rows, ]
            }
            
            # Select relevant columns for display
            display_columns <- c("uniprot_acc", "log2FC", "raw_pvalue", "fdr_qvalue")
            current_results <- current_results |>
              dplyr::select(any_of(display_columns))
            
            DT::datatable(
              current_results,
              options = list(
                pageLength = 25,
                scrollX = TRUE,
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel')
              ),
              extensions = 'Buttons'
            ) |>
            DT::formatRound(columns = c("log2FC", "raw_pvalue", "fdr_qvalue"), digits = 4)
          } else {
            # No results for this contrast
            DT::datatable(data.frame(Message = "No results available for selected contrast"))
          }
        } else {
          # No DE results available
          DT::datatable(data.frame(Message = "No DE analysis results available"))
        }
        
      }, error = function(e) {
        message(sprintf("*** ERROR in DE results table: %s", e$message))
        DT::datatable(data.frame(Message = sprintf("Error: %s", e$message)))
      })
    })
    
    # Summary statistics
    output$de_summary_stats <- shiny::renderText({
      shiny::req(input$table_contrast, de_data$de_results_list)
      
      tryCatch({
        if (!is.null(de_data$de_results_list$de_proteins_long)) {
          # Filter for selected contrast
          current_results <- de_data$de_results_list$de_proteins_long |>
            dplyr::filter(comparison == input$table_contrast)
          
          if (nrow(current_results) > 0) {
            total_genes <- nrow(current_results)
            significant <- sum(current_results$fdr_qvalue < input$de_q_val_thresh, na.rm = TRUE)
            up_reg <- sum(current_results$fdr_qvalue < input$de_q_val_thresh & 
                         current_results$log2FC > input$treat_lfc_cutoff, na.rm = TRUE)
            down_reg <- sum(current_results$fdr_qvalue < input$de_q_val_thresh & 
                           current_results$log2FC < -input$treat_lfc_cutoff, na.rm = TRUE)
            
            paste(
              sprintf("Total genes: %d", total_genes),
              sprintf("Significant (q < %.3f): %d", input$de_q_val_thresh, significant),
              sprintf("Up-regulated: %d", up_reg),
              sprintf("Down-regulated: %d", down_reg),
              sprintf("Fold-change cutoff: %.2f", input$treat_lfc_cutoff),
              sep = "\n"
            )
          } else {
            "No results available for selected contrast"
          }
        } else {
          "No DE analysis results available"
        }
        
      }, error = function(e) {
        sprintf("Error calculating statistics: %s", e$message)
      })
    })
    
    # Download handler for results
    output$download_de_results <- shiny::downloadHandler(
      filename = function() {
        paste0("DE_results_", Sys.Date(), ".zip")
      },
      content = function(file) {
        # Placeholder for creating downloadable zip file
        writeLines("Differential expression results would be packaged here", file)
      }
    )
    
    # Return DE data for potential use by parent module
    return(de_data)
  })
} 